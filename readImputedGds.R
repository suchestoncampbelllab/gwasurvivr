readImputedGds <- function(gdsfile, scanfile, snpfile, infofile){
        library(GWASTools)
        library(magrittr)
        library(dplyr)
        library(data.table)
        library(SummarizedExperiment)
        
        # read genotype
        gds <- GdsGenotypeReader(gdsfile)
        # read in snp data
        snpAnnot <- getobj(snpfile)
        # read scan
        scanAnnot <- getobj(scanfile)
        # put into GenotypeData coding 
        genoData <- GenotypeData(gds,
                                 snpAnnot=snpAnnot,
                                 scanAnnot=scanAnnot)

        # store genotype, sample info and, and snp info
        genotypes <- getGenotype(genoData)
        
        # grab map data and remove the row number column
        snp <- getAnnotation(getSnpAnnotation(genoData)) 
        # grab sample file data
        scan <- getAnnotation(getScanAnnotation(genoData))
        
        # read in info table
        infofile <- read.table(infofile, header=T, stringsAsFactors = F)
        colnames(infofile) <- c("snp_id",
                                "rs_id",
                                "position",
                                "exp_freq_a1",
                                "info",
                                "certainty",
                                "type",
                                "info_type0", 
                                "concord_type0",
                                "r2_type0")

        infofile <- infofile %>% dplyr::select(snp_id, rs_id, exp_freq_a1, info, certainty) %>%  # only grab columns of interest
        dplyr::rename(rsID=rs_id, snp=snp_id) # need to be able to generalize the following renaming in function
        
        # merge snp file with info file
        snp <- snp %>% dplyr::left_join(infofile) %>%
                mutate(end=position) %>% # create 'end' position -- since these are SNPs they are the same 
                dplyr::rename(start=position, chr=chromosome) %>% # rename position to start so we can convert into GRanges object
                makeGRangesFromDataFrame(keep.extra.columns=T) # create GRanges object because needed for rowRanges in se
        
        # parse sample file ... but might not do this
        scan <- scan %>% select(ID_1, ID_2, missing, sex) %>% 
                dplyr::rename(sex_fromSampleFile="sex") %>%
                DataFrame()
        
        # assign names (so summarizedexperiment object has colnames and rownames filled out)
        dimnames(genotypes) <- list(snp$rsID, scan$ID_2)
        
        # put into summarizedexperiment 
        se <- SummarizedExperiment(assays=list(input.files=genotypes),
                                   colData=scan,
                                   rowRanges=snp)
        
        # close gds files so you can reopen them
        unlink(gdsfile, scanfile, snpfile)
        close(gds)
        return(se)
}

