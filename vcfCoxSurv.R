library(VariantAnnotation)
vcf <- VcfFile("chr21.dose.vcf.gz")
vcf <- VcfFile("chr21.dose.vcf.gz", yieldSize=5000)
scanVcfHeader(vcf)
yieldTabix(vcf)
scanVcfHeader(vcf)

info(scanVcfHeader(vcf))
vcf.file <- readVcf(vcf)

open(vcf.file)


## write function

vcfCoxSurv <- function(vcf.file, chunk.size){
        open(vcf.file)
        repeat{
                vcf <- VcfFile(vcf.file, yieldSize=chunk.size)
                data <- readVcf(vcf)
                if(nrow(data==0)){
                        break
                }
                
                # grab data so that it is in SE
                
                # write a survival function
                
                # save output in chunks        
                
        }
        close(vcf.file)

}




repeat {
        data <- readVcf(vcf)
        if(nrow(data==0)) break
        ## do whatever with chunk
        # right here you'd fit 100 survival models and accumulate results
        # save results and move onto next chunk
        # next time through loop assign next chunk to same variable
        ## so the previous memory is free to be garbage collected
}




newVcf <- readVcf(vcf,param=ScanVcfParam(geno="DS"))

newVcf <- readVcf(vcf,param=ScanVcfParam(geno="DS", fixed=c("REF", "ALT"), info=c("MAF", "R2")))

# just the DS matrix
vcf.ds <- readGeno(vcf, "DS")