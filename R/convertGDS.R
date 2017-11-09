convertImputeGds <- function(chunk.impute, chunk.sample, chr, outfile.name){
        gdsfile <- paste0(outfile.name, ".gds")
        snpfile <- paste0(outfile.name, ".snp.rdata")
        scanfile <- paste0(outfile.name, ".scan.rdata")
        imputedDosageFile(input.files=chunk.impute, chunk.sample),
                          filename=gdsfile,
                          chromosome=as.numeric(chr),
                          input.type="IMPUTE2",
                          input.dosage=F,
                          file.type="gds",
                          snp.annot.filename = snpfile,
                          scan.annot.filename = scanfile)
}
