subl#' Convert IMPUTE2 data to GDS format
#'
#' Convert flat file in genotype format (.gen or .impute) to compressed GDS format. This should be done as the first step of performing survival analysis with IMPUTE2 files. 
#'
#' @param chunk.impute character(1) file name of impute genotype file
#'
#' @param chunk.sample character(1) file name of impute sample file
#'
#' @param chr integer(1) defining the chromosome that is being compressed 
#'
#' @param outfile.name character(1) string defining output file names WITHOUT extension
#' 
#' @return
#' Three files: 1) GDS file, 2) SNP RData file, 3) Scan RData file
#' 
#' @example
#' chunk.impute <- "chr21.25000005-25500000.impute"
#' chunk.sample <- "chr21.25000005-25500000.impute.sample"
#' chr <- 21
#' outfile <- "chr21.25000005-25500000"
#' 
#' convertImputeGds(chunk.impute, chunk.sample, chr, outfile)
#' 
#' @export


convertImputeGds <- function(chunk.impute, chunk.sample, chr, outfile.name){
        gdsfile <- paste0(outfile.name, ".gds")
        snpfile <- paste0(outfile.name, ".snp.rdata")
        scanfile <- paste0(outfile.name, ".scan.rdata")
        imputedDosageFile(input.files=c(chunk.impute, chunk.sample),
                          filename=gdsfile,
                          chromosome=as.numeric(chr),
                          input.type="IMPUTE2",
                          input.dosage=F,
                          file.type="gds",
                          snp.annot.filename = snpfile,
                          scan.annot.filename = scanfile)
}
