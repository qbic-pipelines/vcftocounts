#!/usr/bin/env Rscript
'VCF to count matrix converter

Usage:
    vcf2counts.R --help
    vcf2counts.R --output=<ofile> VCF

Options:
    -h, --help                  help screen
    -o, --output=<ofile>        output file name [default: mat.csv]

Arguments:
    VCF                     input vcf file
'->doc

suppressMessages(library(VariantAnnotation, warn.conflicts = FALSE, quietly=TRUE))
suppressMessages(library(docopt, warn.conflicts = FALSE, quietly=TRUE))
suppressMessages(library(Matrix, warn.conflicts = FALSE, quietly=TRUE))

generateMatrixfromVCF <- function(VCF, ofile) {
    # Read in VCF file
    vcfobj <- readVcf(VCF)
    # Convert genotype to SNP matrix
    genomat <- geno(vcfobj)$GT

    variantmat <- apply(genomat, c(1, 2), function(x) {
        xstrip <- gsub("[[:punct:]]", "", x)
        if (xstrip == "11") {
            return (2)
        } else if (xstrip %in% c("01", "10")) {
            return (1)
        } else if (xstrip %in% c("00", "")) {
            return (0)
        } else {
            return (NA)
        }
    })
    write.csv(variantmat, file = ofile)
}

opt <- docopt(doc)

generateMatrixfromVCF(opt$VCF, opt[["--output"]])
