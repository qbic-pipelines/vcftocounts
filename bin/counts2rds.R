#!/usr/bin/env Rscript
## counts2rds.R -- merge vcf2counts per-sample CSVs into the objects the compareDiversity
## analysis reads, i.e. what its own `prepare` step writes
## (compareDiversity R/steps.R::step_prepare, R/hill_calculations.R::prepare_null):
##
##   null.rds      prepared null: integer 0/1/2 call matrix plus the variant/sample index
##   variants.rds  id, chr, pos of every variant, in the null's variant order (plan, regions)
##   info.rds      samples, code-to-description maps, pipelines, row and variant counts
##   q2.rds        call-pattern codes for the analysis' --fastQ2 path (q = 2 only). Written
##                 unless the number of call patterns exceeds 1e8, which is the same limit
##                 step_prepare applies; above it --fastQ2 cannot be used at all.
##
## Every file is written as <prefix>.<name> when --prefix is given, so several sample sets
## can share one directory. The analysis reads the bare names from --prepDir, so hand it a
## directory holding exactly one set, with the prefix stripped.
##
## Genotype coding: 0 = hom ref, 1 = het, 2 = hom alt. An uncalled genotype (variant absent
## from the sample's CSV, or the caller left it uncalled/multi-allelic) becomes 0: the Hill
## numbers count absence, not missing data. chrX/chrY are dropped unless --include-sex, and
## samples without a population or superpopulation label are dropped, because they cannot
## be grouped - the same rules the analysis' own count-directory reader applies.
##
## Usage:
##   counts2rds.R --counts counts/ --metadata 1000g.tsv --population populations.tsv \
##       --out-dir out --threads 8

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

## ===========================================================================
## Metadata
## ===========================================================================

## Source column -> analysis column, first match wins. This reads both an analysis-shaped
## table and the 1000G sequence index without a flag. Only the five columns read_demo
## produces are kept: the caller columns are derived as "everything else", so a stray
## character column would be taken for a caller.
ALIASES <- list(
  sample_name = c("Run", "SAMPLE_NAME", "sample_name", "Sample", "sample"),
  population_code = c("population", "POPULATION", "population_code", "Population", "Population Code"),
  super_population_code = c("Super_Population_Code", "super_population_code",
                            "SUPER_POPULATION_CODE", "super_population", "Super Population"),
  population_description = c("Population_Description", "population_description",
                             "POPULATION_DESCRIPTION", "Population Description"),
  super_population_description = c("Super_Population_Description", "super_population_description",
                                   "SUPER_POPULATION_DESCRIPTION", "Super Population Description")
)
BLANKS <- c("", "na", "nan", "none", "null", ".")     # strings that mean "no value"
## populations.tsv carries no superpopulation description, only the code; labelling the
## five known codes keeps the report readable. Unknown codes stay NA and are shown as-is.
SUPER_DESC <- c(AFR = "African", AMR = "Admixed American", EAS = "East Asian",
                EUR = "European", SAS = "South Asian")

read_table <- function(path) {
  ## The separator is sniffed from the header: the sequence index is tab-separated,
  ## the analysis-shaped tables and populations.tsv are not.
  fread(path, sep = if (grepl("\t", readLines(path, n = 1L), fixed = TRUE)) "\t" else "auto",
        showProgress = FALSE, colClasses = "character")
}

read_demo <- function(path, lookup = NULL) {
  ## Same five columns, and the same filling rules, as compareDiversity R/data_io.R::read_demo.
  raw <- read_table(path)
  ## Five columns of the file's full length: a missing source column becomes all-NA, never a
  ## short vector that a later assignment cannot fit.
  demo <- as.data.table(lapply(names(ALIASES), function(col) {
    hit <- intersect(ALIASES[[col]], names(raw))[1]
    if (is.na(hit)) rep(NA_character_, nrow(raw)) else as.character(raw[[hit]])
  }))
  setnames(demo, names(ALIASES))
  demo <- demo[!is.na(sample_name) & !sample_name %in% BLANKS]
  demo <- unique(demo, by = "sample_name")        # the sequence index has many rows per sample
  for (col in names(demo)) {
    v <- demo[[col]]
    v[v %in% BLANKS] <- NA_character_
    set(demo, j = col, value = v)
  }
  if (!is.null(lookup)) {
    lk <- read_table(lookup)
    code <- intersect(ALIASES$population_code, names(lk))[1]
    sup <- intersect(ALIASES$super_population_code, names(lk))[1]
    if (!is.na(code) && !is.na(sup)) {
      ## Columns are pulled by [[, never selected by name inside i: a header like
      ## "Population Code" is parsed as an expression there and comes back as a logical.
      ## One row per code (the file lists a code per phase), matched, not joined.
      ok <- !is.na(lk[[code]]) & !lk[[code]] %in% BLANKS & !is.na(lk[[sup]])
      ref <- data.table(population_code = as.character(lk[[code]][ok]), super = as.character(lk[[sup]][ok]))
      ref <- ref[!duplicated(population_code)]
      miss <- is.na(demo$super_population_code)
      if (any(miss) && nrow(ref) > 0L) {
        idx <- match(demo$population_code, ref$population_code)
        set(demo, which(miss), "super_population_code", ref$super[idx[miss]])
      }
    }
  }
  ## A population's superpopulation and descriptions are only carried on some rows of the
  ## demographics file; fill them from the siblings that have them.
  fill <- function(target) {
    have <- !is.na(demo[[target]])
    if (!any(!have)) return(invisible(NULL))
    ref <- demo[have][!duplicated(population_code)]
    if (nrow(ref) == 0L) return(invisible(NULL))
    idx <- match(demo$population_code, ref$population_code)
    set(demo, which(!have), target, ref[[target]][idx[!have]])
  }
  for (col in c("super_population_code", "super_population_description", "population_description"))
    fill(col)
  ## Labels for the five superpopulation codes, for reports. Only fills what is still missing.
  known <- SUPER_DESC[demo$super_population_code]
  fill_known <- is.na(demo$super_population_description) & !is.na(known)
  if (any(fill_known)) set(demo, which(fill_known), "super_population_description", known[fill_known])
  ## Ungroupable samples are dropped here, as read_demo does, so the null never holds them.
  demo[!is.na(population_code) & !is.na(super_population_code)]
}

## ===========================================================================
## Count tables
## ===========================================================================

read_sample <- function(path, exclude, include_sex) {
  x <- fread(path, showProgress = FALSE, colClasses = list(character = 1))
  sample_name <- sub("\\.csv$", "", basename(path))
  setnames(x, 1L, "variant_called")
  ## vcf2counts may prefix the caller columns with "<sample>_"; the analysis wants bare names.
  pre <- paste0("^", sample_name, "_")
  setnames(x, names(x)[-1L], sub(pre, "", names(x)[-1L]))
  callers <- setdiff(setdiff(names(x), "variant_called"), exclude)
  x <- x[, c("variant_called", callers), with = FALSE]
  n_filled <- 0L
  for (p in callers) {
    v <- suppressWarnings(as.integer(x[[p]]))
    n_filled <- n_filled + sum(is.na(v) | v < 0L)
    v[is.na(v) | v < 0L] <- 0L        # uncalled becomes 0: absence, not missing data
    set(x, j = p, value = v)
  }
  if (!include_sex) x <- x[!grepl("^chr[XY]:", variant_called)]
  x <- x[rowSums(x[, callers, with = FALSE]) > 0L]
  x$sample_name <- sample_name
  attr(x, "n_filled") <- n_filled
  x
}

read_counts <- function(input_dir, exclude, include_sex, threads, demo) {
  files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  tables <- if (threads > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(files, read_sample, exclude, include_sex, mc.cores = threads) else
    lapply(files, read_sample, exclude, include_sex)
  ## Every sample arrives with the same callers; the caller set is read off the tables.
  callers <- setdiff(names(tables[[1L]]), c("variant_called", "sample_name"))
  n_filled <- sum(vapply(tables, attr, numeric(1), "n_filled"))
  d <- rbindlist(tables, use.names = TRUE)
  rm(tables); invisible(gc())
  ## Caller columns have to be integer before prepare_null indexes or maxes them.
  for (p in callers) set(d, j = p, value = as.integer(d[[p]]))
  ## One variant can appear under several rsIDs; the strongest call wins. Calls are already
  ## 0-filled, so there is no NA case to decide here.
  d <- d[, lapply(.SD, max), by = c("sample_name", "variant_called"), .SDcols = callers]
  ## Ungroupable samples cannot enter the null.
  d <- d[sample_name %in% demo$sample_name]
  list(data = d, callers = callers, n_filled = n_filled)
}

## ===========================================================================
## Prepare (ported from compareDiversity R/hill_calculations.R, kept verbatim where
## possible so the objects stay identical to the analysis' own prepare step)
## ===========================================================================

## Code-to-label map for the report. Both the description and the code are character here
## (read_demo stores them that way), which is what setNames needs.
create_code_mapping <- function(samples, code_col, desc_col) {
  mapping <- setNames(samples[[desc_col]], samples[[code_col]])
  mapping[!duplicated(names(mapping))]
}

prepare_null <- function(null_data, demo) {
  vlev <- levels(null_data$variant_called)
  chr <- sub(":.*", "", vlev)
  pos <- suppressWarnings(as.numeric(sub("^[^:]*:([0-9]+).*", "\\1", vlev)))
  if (anyNA(pos))
    stop("Cannot parse a position from variant_called (expected chr:pos_ref/alt), e.g. ",
         vlev[which(is.na(pos))[1]], call. = FALSE)
  ord <- order(chr, pos, vlev, method = "radix")   # radix: no locale collation of ~24M strings
  rank <- integer(length(vlev)); rank[ord] <- seq_along(ord)
  code <- rank[as.integer(null_data$variant_called)]
  o <- order(code)
  row_var <- code[o]; rm(code)

  pipes <- setdiff(names(null_data), c("sample_name", "variant_called"))
  M <- matrix(0L, nrow(null_data), length(pipes), dimnames = list(NULL, pipes))
  for (j in seq_along(pipes)) {                    # one pipeline at a time: no full double copy
    v <- null_data[[pipes[j]]][o]
    v[is.na(v) | v < 0L] <- 0L
    M[, j] <- v
  }
  row_sample <- null_data$sample_name[o]
  rm(o)
  samples <- demo[match(levels(row_sample), demo$sample_name), ]
  var_len <- tabulate(row_var, nbins = length(vlev))
  nd <- list(
    M = M, row_var = row_var, var_start = cumsum(c(1L, var_len[-length(var_len)])),
    var_len = var_len, variant_chr = chr[ord], row_sample = row_sample,
    sample_pop = factor(samples$population_code),           # indexed by as.integer(row_sample)
    sample_super = factor(samples$super_population_code),
    samples = samples[, c("sample_name", "population_code", "super_population_code")],
    pop_desc = create_code_mapping(samples, "population_code", "population_description"),
    super_desc = create_code_mapping(samples, "super_population_code", "super_population_description")
  )
  list(nd = nd, variant_ids = vlev[ord])
}

q2_records <- function(nd) {
  pipes <- colnames(nd$M); k <- length(pipes)
  ## Same guard as step_prepare, and in the same order: the product is computed in double
  ## before any as.integer(), which would overflow to NA for more than a few pipelines.
  base <- max(nd$M) + 1L; npop <- nlevels(nd$sample_pop)
  if (as.numeric(base)^k * npop > 1e8) return(NULL)
  npat <- as.integer(base^k)
  pop <- as.integer(nd$sample_pop)[as.integer(nd$row_sample)]
  code <- (pop - 1L) * npat + 1L
  for (j in seq_len(k)) code <- code + nd$M[, j] * as.integer(base^(j - 1))
  ## Within each variant, rows are ordered by superpopulation then population, so each
  ## (variant, group) is one contiguous run for the by-variant sums.
  pops <- levels(nd$sample_pop)
  pop_super <- nd$samples$super_population_code[match(pops, nd$samples$population_code)]
  code <- code[order(nd$row_var, match(pop_super, sort(unique(pop_super)))[pop], code,
                     method = "radix")]
  nv <- length(nd$var_len)
  pats <- as.matrix(expand.grid(rep(list(0:(base - 1L)), k)))   # row p = pattern p, first pipeline fastest
  colnames(pats) <- pipes
  list(code = code, var_len = nd$var_len, var_start = cumsum(c(1L, nd$var_len[-nv])),
       variant_chr = nd$variant_chr, pats = pats, npat = npat, pops = pops, pop_super = pop_super)
}

## ===========================================================================
## main
## ===========================================================================

main <- function() {
  opt <- parse_args(OptionParser(option_list = list(
    make_option(c("-i", "--counts"), help = "Directory with one vcf2counts CSV per sample"),
    make_option(c("-m", "--metadata"), help = "Sample metadata: demographics table or 1000G sequence index"),
    make_option(c("-p", "--population"), default = NULL,
                help = "Table keyed on population code giving the superpopulation, e.g. populations.tsv"),
    make_option(c("-o", "--out-dir"), default = ".",
                help = "Where the files go, no subdirectory is made [default: the working dir]"),
    make_option(c("-x", "--prefix"), default = NULL,
                help = "Prefix for every output file; omit to write unprefixed names"),
    make_option(c("-e", "--exclude-callers"), default = NULL,
                help = "Comma-separated caller columns to drop, e.g. DeepVariant,GATK"),
    make_option(c("-s", "--include-sex"), action = "store_true", default = FALSE,
                help = "Keep chrX/chrY variants (dropped by default)"),
    make_option(c("-t", "--threads"), type = "integer", default = 1L,
                help = "Samples read in parallel [default %default]")
  )), positional_arguments = FALSE, convert_hyphens_to_underscores = TRUE)

  ## Everything else is assumed to arrive cleaned and checked from the pipeline.
  for (required in c("counts", "metadata"))
    if (is.null(opt[[required]]))
      stop("Missing required option: --", gsub("_", "-", required), call. = FALSE)
  exclude <- if (is.null(opt$exclude_callers)) character(0) else
    trimws(strsplit(opt$exclude_callers, ",", fixed = TRUE)[[1]])

  demo <- read_demo(opt$metadata, opt$population)
  loaded <- read_counts(opt$counts, exclude, opt$include_sex, max(1L, opt$threads), demo)
  ## prepare_null works off factor codes; factor levels are ids only, the order comes from
  ## (chr, pos, id). sample_name has to be a factor too: q2_records indexes sample_pop by it.
  d <- loaded$data
  ## Same shapes as the analysis' own tidy_null output, because prepare_null depends on them:
  ## levels give the ids, as.integer gives the codes, and row_sample's levels fix which
  ## population each row belongs to. tidy_null hashes, i.e. levels in order of appearance;
  ## base factor() sorts, which reorders the samples but not the variant rank (prepare_null
  ## sorts variants by chr, pos, id itself).
  d[, `:=`(variant_called = factor(variant_called), sample_name = factor(sample_name))]

  message("Preparing the null (sorting by variant)...")
  p <- prepare_null(d, demo)
  rm(d); invisible(gc())

  info <- list(samples = p$nd$samples, pop_desc = p$nd$pop_desc, super_desc = p$nd$super_desc,
               pipelines = colnames(p$nd$M), null_rows = nrow(p$nd$M),
               n_variants = length(p$variant_ids), nullCSV = NA_character_,
               countDirs = normalizePath(opt$counts), excludeCallers = exclude,
               includeSex = opt$include_sex, demoInfo = normalizePath(opt$metadata))
  variants <- data.table(id = p$variant_ids, chr = p$nd$variant_chr,
                         pos = as.numeric(sub("^[^:]*:([0-9]+).*", "\\1", p$variant_ids)))
  ## With --prefix every file is named <prefix>.<name>, so several sample sets can share
  ## one directory. The analysis expects the bare names in --prepDir.
  path <- function(f) file.path(opt$out_dir, if (is.null(opt$prefix)) f else paste0(opt$prefix, ".", f))
  ## A file that this run does not produce must not survive from an earlier one: the
  ## analysis prefers q2.rds when it exists, so a stale copy would silently win.
  drop_stale <- function(f) if (file.exists(f)) file.remove(f)
  saveRDS(p$nd, path("null.rds"), compress = FALSE)      # shard workers read it whole
  saveRDS(variants, path("variants.rds"), compress = FALSE)
  saveRDS(info, path("info.rds"))
  ## Always written: --fastQ2 is the cheap shard path, and the analysis would rebuild it from
  ## the null anyway. Only the pattern ceiling can prevent it, and then --fastQ2 is unusable
  ## for this dataset regardless of what is written here.
  q2 <- q2_records(p$nd)
  if (is.null(q2)) {
    drop_stale(path("q2.rds"))
    warning("Too many call patterns for q2.rds; --fastQ2 is unavailable for this dataset",
            call. = FALSE)
  } else {
    saveRDS(q2, path("q2.rds"), compress = FALSE)
  }

  ## What went in.
  fwrite(data.table(metric = c("n_samples", "n_callers", "callers", "n_rows", "n_variants",
                               "n_na_filled", "n_populations", "n_superpopulations",
                               "include_sex", "metadata_file", "population_file"),
                    value = c(nrow(p$nd$samples), length(loaded$callers),
                              paste(loaded$callers, collapse = "|"), nrow(p$nd$M),
                              length(p$variant_ids), loaded$n_filled, nlevels(p$nd$sample_pop),
                              nlevels(p$nd$sample_super), opt$include_sex, opt$metadata,
                              if (is.null(opt$population)) "none" else opt$population)),
         path("summary.csv"))

  message(sprintf("Prepared null: %d row(s), %d variant(s), %d sample(s), %d pipeline(s) -> %s",
                  nrow(p$nd$M), length(p$variant_ids), nrow(p$nd$samples), ncol(p$nd$M),
                  normalizePath(opt$out_dir)))
  if (loaded$n_filled > 0L)
    message(sprintf("Uncalled genotypes filled with 0: %d cell(s)", loaded$n_filled))
  invisible(NULL)
}

if (sys.nframe() == 0L) main()
