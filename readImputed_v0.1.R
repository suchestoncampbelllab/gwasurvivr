library(SummarizedExperiment)
library(GWASTools)
library(gdsfmt)
library(Biobase)
library(snpStats)
library(snow)
library(survival)
library(dplyr)
library(data.table)


## import imputed dosage files
path <- "~/GoogleDrive/OSU_PHD/survivR/data"
setwd(path)

gdsfile <- "chr7-55000000-60000000.v2.gds"
snpfile <- "chr7-55000000-60000000.snp.RData"
scanfile <- "chr7-55000000-60000000.scan.RData"
logfile <- "chr7-55000000-60000000.log"

input.files <- c(paste0(path, "/", "BMT093013_forImpute.chr7-55000000-60000000.impute2"),
                   paste0(path, "/", "BMT093013_forImpute.chr7-55000000-60000000.impute2_samples"))

chromosome <- 7

input.type <- "IMPUTE2"
input.dosage <- F

readImpute <- function(input.files, info.file, chromosome, input.type, input.dosage, select.snps=NULL, ...){

# gdsfile <- tempfile()
# snpfile <- tempfile()
# scanfile <- tempfile()
# 
# if(is.null(chromosome)){
#         chromosome <- NA
# }

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
dim(dat$genotypes)
dimnames(dat$genotypes) <- list(dat.gr$rsID, dat$fam$sampleID)




# grab map data and remove the row number column
dat$map <- genoData@snpAnnot@data[,-1]
# grab fam data
dat$fam <- genoData@scanAnnot@data


# read in info table
infoFile <- read.table("BMT093013_forImpute.chr7-55000000-60000000.impute2_info", header=T)
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
dat$fam <- dat$fam[,"sampleID"]
dat$fam <- sapply(strsplit(dat$fam, " "), "[", 2)
dat$fam <- DataFrame(sampleID=dat$fam)


# if(input.type=="BEAGLE"){
#         snpid <- dat$map$snp
# }else if(input.type=="IMPUTE2"){
#         snpid <- dat$map$rsID
# }else if(input.type=="MaCH"){
#         snpid <- dat$map$snp
# }

# create SummarizedExperiment object

se <- SummarizedExperiment(assays=list("chr7-55000000-60000000.impute2"=dat$genotypes),
                                colData=dat$fam,
                                rowRanges=dat.gr)


if(!is.null(select.snps)){
        se <- se[select.snps,]
}
close(gds)
unlink(c(gdsfile, snpfile, scanfile))
return(se)
}

# library(survival)



