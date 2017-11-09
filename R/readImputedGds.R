readImputedGds <- function(gdsfile, scanfile, snpfile, infofile){
        library(GWASTools)
        library(magrittr)
        library(dplyr)
        library(data.table)
        library(SummarizedExperiment)
        
        # read genotype
        gds <- GdsGenotypeReader(gdsfile)
        # close gds file on exit of the function
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
        snp <- getAnnotation(getSnpAnnotation(genoData)) %>%
                dplyr::rename(snp.index=snpID,
                              rsid=rsID,
                              snpid=snp)
        # grab sample file data
        scan <- getAnnotation(getScanAnnotation(genoData))
        
        # read in info table
        infofile <- read.table(infofile, header=T, stringsAsFactors = F)
        colnames(infofile) <- c("snpid",
                                "rsid",
                                "position",
                                "exp_freq_a1",
                                "info",
                                "certainty",
                                "type",
                                "info_type0", 
                                "concord_type0",
                                "r2_type0")

        infofile <- infofile %>%
                dplyr::select(snpid, rsid, exp_freq_a1, info, certainty)
        
        # merge snp file with info file
        snp <- snp %>%
                dplyr::left_join(infofile) %>%
                mutate(end=position) %>% # create 'end' position -- since these are SNPs they are the same 
                dplyr::rename(start=position,
                              chr=chromosome, 
                              allele0=alleleA, 
                              allele1=alleleB) %>% # rename position to start so we can convert into GRanges object
                makeGRangesFromDataFrame(keep.extra.columns=T) # create GRanges object because needed for rowRanges in se
        
        # parse sample file ... but might not do this
        scan <- scan %>%
                dplyr::select(ID_1, ID_2, missing, sex) %>% 
                dplyr::rename(sex_fromSampleFile="sex") %>%
                DataFrame()
        
        # assign names (so summarizedexperiment object has colnames and rownames filled out)
        dimnames(genotypes) <- list(snp$rsid, scan$ID_2)
        
        chunk.name <- as.character(gsub(".gds", "", gdsfile))
        
        # put into summarizedexperiment 
        se <- SummarizedExperiment(assays=setNames(list(genotypes), chunk.name),
                                   colData=scan,
                                   rowRanges=snp)
        
        close(gds)
        
        return(se)
}

