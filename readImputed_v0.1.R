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
input.dosage=F

#read.imputed <- function(input.files, chromosome, input.type, input.dosage, select.snps=NULL, ...){
gdsfile <- tempfile()
snpfile <- tempfile()
scanfile <- tempfile()

if(is.null(chromosome)){
        chromosome <- NA
}

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




# grab map data and remove the row number column
dat$map <- genoData@snpAnnot@data[,-1]
# grab fam data
dat$fam <- genoData@scanAnnot@data



# read in info table
infoFile <- read.table("BMT093013_forImpute.chr7-55000000-60000000.impute2_info", header=T)
# only grab columns of interest
infoFile <- infoFile[,c("rs_id", "exp_freq_a1", "info", "certainty")]
# need to be able to generalize the following renaming in function
infoFile <- infoFile %>% rename(rsID=rs_id)
# merge 'map' file with info file
dat$map <- dat$map %>% left_join(infoFile)
dat$map$end <- dat$map$position 
dat$map <- dat$map %>% rename(start=position, chr=chromosome) 

# create GRanges object because needed for rowRanges in se
dat.gr <- makeGRangesFromDataFrame(dat$map, keep.extra.columns = T)



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
# rownames(dat$genotypes) <- rownames(dat$map) <- snpid
# colnames(dat$genotypes) <- rownames(dat$fam) <- dat$fam$sampleID
# 


#####################        
# subset for testing
#####################

dat.sub <- list()
dat.sub$genotypes <- dat$genotypes[1:1000,1:6805]
dat.sub$rowData <- dat.gr[1:1000,]
dat.sub$colData <- dat$fam
dimnames(dat.sub$genotypes) <- list(dat.sub$rowData$rsID, dat.sub$colData$sampleID)
head(dat$fam)


se <- SummarizedExperiment(assays=list("chr7-55000000-60000000.impute2"=dat.sub$genotypes),
                      colData=dat.sub$colData,
                      rowRanges=dat.sub$rowData)

############################
## load phenotype data
############################

pheno <- fread("~/GoogleDrive/OSU_PHD/mhcproject/BMT_PHENOTYPE_FILE_QC_allclinical_5_5_17.txt")

#######################
## clean phenotype data
#######################

# remove SI.rsxxxx I.rsxxxx ST, T.rsxxxxxx columns
pheno <- pheno %>% 
        select(-c(grep("I.r", colnames(pheno)),
                  grep("SI.", colnames(pheno)),
                  grep("ST.", colnames(pheno)),
                  grep("T.", colnames(pheno))))


pheno <- pheno %>% 
        select(-PC1.AA.c1:-PC10.HIS.c1)

pheno$IID.c1.c2 <- NA 
pheno$IID.c1.c2[!is.na(pheno$IID.EA.c1)] <- pheno$IID.EA.c1[!is.na(pheno$IID.EA.c1)]
pheno$IID.c1.c2[!is.na(pheno$IID.EA.c2)] <- pheno$IID.EA.c2[!is.na(pheno$IID.EA.c2)]

head(pheno)

# change those 3 patients IDs that are swapped 
pheno$IID.c1.c2 <- pheno$IID.c1.c2 %>% dplyr::recode(`G_BMT_OmniExp_BMT-125_0085017851`="G_BMT_OmniExp_BMT-125_0085017823",
                                              `G_BMT_OmniExp_BMT-152_0085017639`="G_BMT_OmniExp_BMT-152_1023993633",
                                              `G_BMT_OmniExp_BMT-147_0105388831`="G_BMT_OmniExp_BMT-119_0074380675") 

# create dummy columns for trm covariates

trm.cov <- pheno %>%
        select(IID.c1.c2, intxsurv_1Y, trm_overall, age, bmi_cat, graftype) %>%
        mutate(bmiOBS=ifelse(bmi_cat=="obese", 1, 0),
               bmiOVWT=ifelse(bmi_cat=="overweight", 1, 0),
               PBlood=ifelse(graftype=="Peripheral blood", 1, 0)) %>%
        select(-bmi_cat, -graftype) %>%
        rename(sampleID=IID.c1.c2) %>%
        data.table %>%
        na.omit

length(colData(se)$sampleID)

se.sub <- se[,colnames(se) %in% trm.cov$sampleID]



data.frame(colData(se)) %>% left_join(trm.cov) %>% na.omit %>% nrow


# if(!is.null(select.snps)){
#         se <- se[select.snps,]
# }
# close(gds)
# unlink(c(gdsfile, snpfile, scanfile))
# return(se)
#}

# library(survival)



