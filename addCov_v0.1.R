

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

addCov <- function(se, covFile){
        se <- se[,colnames(se) %in% covFile[[1]]]
        colData(se) <- merge(colData(se), covFile)
        return(se)
}


se.cov <- addCov(se.1, trm.cov)



