library(SummarizedExperiment)
library(GWASTools)
library(gdsfmt)
#library(Biobase)
#library(snpStats)
library(snow)
library(survival)
library(dplyr)
library(data.table)

readImpute <- function(input.files, info.file, chromosome=NULL, input.type, input.dosage, select.snps=NULL, ...){
        
        # creates temp files on disk
        gdsfile <- tempfile()
        snpfile <- tempfile()
        scanfile <- tempfile()
        
        # if no chromosome is given, assign it to NA
        if(is.null(chromosome)){
                chromosome <- NA
        }
        
        # GWASTools function (cite GWASTools)
        imputedDosageFile(input.files=input.files,
                          filename=gdsfile,
                          chromosome=chromosome,
                          input.type=input.type,
                          input.dosage=input.dosage,
                          file.type="gds",
                          snp.annot.filename = snpfile,
                          scan.annot.filename = scanfile)
        
        gds <- GdsGenotypeReader(gdsfile)
        scanAnnot <- getobj(scanfile)
        snpAnnot <- getobj(snpfile)
        genoData <- GenotypeData(gds,
                                 scanAnnot=scanAnnot,
                                 snpAnnot=snpAnnot)
        
        
        # store genotype, sample info and, and snp info
        dat <- list()
        dat$genotypes <- getGenotype(genoData)

        # grab map data and remove the row number column
        dat$map <- genoData@snpAnnot@data[,-1]
        # grab sample file data
        dat$sample <- genoData@scanAnnot@data
        
        
        # read in info table
        infoFile <- read.table(infoFile, header=T)
        # only grab columns of interest
        infoFile <- infoFile[,c("snp_id", "rs_id", "exp_freq_a1", "info", "certainty")]
        # need to be able to generalize the following renaming in function
        infoFile <- infoFile %>% rename(rsID=rs_id, snp=snp_id)
        
        # merge 'map' file with info file
        dat$map <- dat$map %>% left_join(infoFile)
        
        # create 'end' position -- since these are SNPs they are the same 
        dat$map$end <- dat$map$position 
        
        # rename position to start so we can convert into GRanges object
        dat$map <- dat$map %>% rename(start=position, chr=chromosome) 
        
        # create GRanges object because needed for rowRanges in se
        dat.gr <- makeGRangesFromDataFrame(dat$map, keep.extra.columns = T)
        
        
        # parse sample file ... but might not do this
        dat$sample <- dat$sample[,c("ID_1", "ID_2", "missing", "sex")]
        dat$sample <- dat$sample %>% rename(sex_fromSampleFile="sex")
        #dat$sample <- sapply(strsplit(dat$sample, " "), "[", 2)
        dat$sample <- DataFrame(dat$sample)
        
        # assign names (so summarizedexperiment object has colnames and rownames filled out)
        dimnames(dat$genotypes) <- list(dat.gr$rsID, dat$sample$ID_2)
        
        
        # if(input.type=="BEAGLE"){
        #         snpid <- dat$map$snp
        # }else if(input.type=="IMPUTE2"){
        #         snpid <- dat$map$rsID
        # }else if(input.type=="MaCH"){
        #         snpid <- dat$map$snp
        # }
        
        # create SummarizedExperiment object
        
        se <- SummarizedExperiment(assays=list(input.files=dat$genotypes),
                                        colData=dat$sample,
                                        rowRanges=dat.gr)
        
        
        if(!is.null(select.snps)){
                se <- se[select.snps,]
        }
        close(gds)
        unlink(c(gdsfile, snpfile, scanfile))
        return(se)
}


# test


library(microbenchmark)

microbenchmark(se.1 <- readImpute(input.files, infoFile, chromosome, input.type, input.dosage=F), times=1)
