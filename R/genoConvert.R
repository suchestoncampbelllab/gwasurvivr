# Convert a VCF genotype field to an ALT-allele dosage matrix (SNP x sample).
#
# gwasurvivr models additive ALT-allele dosage. Imputation servers emit several
# encodings; we support:
#   DS  - ALT dosage already (Minimac/Michigan, Sanger). Used as-is.
#   GT  - hard genotype calls ("0|1", "1/1", ...). Converted to allele count.
#   HDS - per-haplotype dosages ("0.02,0.98"). Summed to a diploid dosage.
# The conversions are fully vectorized (a lookup table for GT) to preserve the
# package's speed on genome-scale data.

genotypeFieldToDosage <- function(mat, field, verbose = TRUE) {
    field <- match.arg(field, c("DS", "GT", "HDS"))
    if (field == "DS") {
        storage.mode(mat) <- "double"
        return(mat)
    }
    if (field == "GT") return(gtToDosage(mat, verbose = verbose))
    if (field == "HDS") return(hdsToDosage(mat))
}

# Hard genotype calls -> ALT allele count (0/1/2). Missing ("." / "./." / ".|.")
# and any non-biallelic call become NA (and trigger a one-time message).
gtToDosage <- function(gt, verbose = TRUE) {
    lut <- c("0/0" = 0, "0|0" = 0,
             "0/1" = 1, "1/0" = 1, "0|1" = 1, "1|0" = 1,
             "1/1" = 2, "1|1" = 2)
    v <- as.vector(gt)
    d <- unname(lut[v])
    missing.call <- v %in% c(".", "./.", ".|.", "") | is.na(v)
    unknown <- is.na(d) & !missing.call
    if (verbose && any(unknown)) {
        message(sum(unknown), " genotype call(s) were not biallelic 0/1 ",
                "(e.g. multiallelic or malformed) and were set to missing.")
    }
    matrix(d, nrow(gt), ncol(gt), dimnames = dimnames(gt))
}

# Per-haplotype dosages -> diploid ALT dosage. VariantAnnotation returns HDS
# either as a character matrix of "a,b" strings or as a 3-D array
# (snp x sample x 2); handle both.
hdsToDosage <- function(hds) {
    if (is.array(hds) && length(dim(hds)) == 3L) {
        d <- hds[, , 1L, drop = TRUE] + hds[, , 2L, drop = TRUE]
        return(matrix(as.numeric(d), nrow(hds), ncol(hds),
                      dimnames = dimnames(hds)[1:2]))
    }
    v <- as.vector(hds)
    parts <- strsplit(v, ",", fixed = TRUE)
    d <- vapply(parts, function(p) {
        if (length(p) == 0L || any(p %in% c(".", ""))) return(NA_real_)
        sum(as.numeric(p))
    }, numeric(1))
    matrix(d, nrow(hds), ncol(hds), dimnames = dimnames(hds))
}
