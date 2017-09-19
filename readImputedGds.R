readImputedGds <- function(gdsfile, snpfile, scanfile, infofile){
        # read genotype
        gds <- GdsGenotypeReader(gdsfile)
        
        # read in snp data
        snpAnnot <- getobj(snpfile)
        
        # read scan 
        scanAnnot <- getobj(scanfile)
        
        # put into GenotypeData coding 
        genoData <- GenotypeData(gds,
                                 scanAnnot=scanAnnot,
                                 snpAnnot=snpAnnot)
        
        # store genotype, sample info and, and snp info
        genotypes <- getGenotype(genoData)
        
        # grab map data and remove the row number column
        snp <- getAnnotation(getSnpAnnotation(genoData)) 
        # grab sample file data
        scan <- getAnnotation(getScanAnnotation(genoData))
        
        # read in info table
        infofile <- read.table(infofile, header=F, stringsAsFactors = F)
        colnames(infoFile) <- c("snp_id",
                                "rs_id",
                                "position",
                                "exp_freq_a1",
                                "info",
                                "certainty",
                                "type",
                                "info_type0", 
                                "concord_type0",
                                "r2_type0")
        
        # only grab columns of interest
        infofile <- infofile[,c("snp_id", "rs_id", "exp_freq_a1", "info", "certainty")]
        
        # need to be able to generalize the following renaming in function
        infofile <- infofile %>% dplyr::rename(rsID=rs_id, snp=snp_id)
        
        # merge 'map' file with info file
        snp <- snp %>% dplyr::left_join(infofile)

        # create 'end' position -- since these are SNPs they are the same 
        snp$end <- snp$position 
        
        # rename position to start so we can convert into GRanges object
        snp <- snp %>% dplyr::rename(start=position, chr=chromosome) 
        
        # create GRanges object because needed for rowRanges in se
        dat.gr <- makeGRangesFromDataFrame(snp, keep.extra.columns = T)
        
        # parse sample file ... but might not do this
        scan <- scan[,c("ID_1", "ID_2", "missing", "sex")]
        scan <- scan %>% dplyr::rename(sex_fromSampleFile="sex")
        
        # convert to DataFrame so it can be put in colData
        scan <- DataFrame(scan)
        
        # assign names (so summarizedexperiment object has colnames and rownames filled out)
        dimnames(genotypes) <- list(dat.gr$rsID, scan$ID_2)
        
        # put into summarizedexperiment 
        se <- SummarizedExperiment(assays=list(input.files=genotypes),
                                   colData=scan,
                                   rowRanges=dat.gr)
        
        return(se)
}