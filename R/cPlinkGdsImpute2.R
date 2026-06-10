loadProcessWrite.PlinkGdsImpute2CoxSurv <-
    function(x, cl, cox.params) {
        writeFileHeadings(x, cox.params = cox.params)
        
        genoData <- getGenoData(x)
        
        results <- runOnChunks(x, genoData, cox.params, cl)
        
        return(list(
            snps_removed = results$snp.drop.n,
            snps_analyzed = results$snp.n
        ))
    }

writeFileHeadings <- function(x, cox.params) {
    # set up columns for output
    write.table(
        t(x$columnHeadings),
        paste0(x$out.file, ".snps_removed"),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE,
        append = FALSE
    )
    
    colnames(x$snp.df) <- x$columnHeadings
    rownames(x$snp.df) <- NULL
    
    cox.out <- getSnpSpikeCoxOut(x$inter.term,
                                 snp.spike = createSnpSpike(x, cox.params),
                                 cox.params,
                                 x$print.covs)
    
    res.cols <- colnames(coxExtract(cox.out,
                                    x$snp.df,
                                    print.covs = x$print.covs))
    
    write.table(
        t(res.cols),
        paste0(x$out.file, ".coxph"),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE,
        append = FALSE
    )
    
}


runOnChunks <- function(x, genoData, cox.params, cl) {
    # number of snps in segment
    snp.start <- 1
    snp.end <- nsnp(genoData)
    # number of dropped snps
    snp.drop.n <- 0
    snp.n <- 0
    
    # get genotypes for certain chunk size
    nsnp.seg <- snp.end - snp.start + 1
    nchunks <- ceiling(nsnp.seg / x$chunk.size)
    
    for (i in seq_len(nchunks)) {
        if (x$verbose)
            message("Analyzing part ", i, "/", nchunks, "...")
        
        # set up chunks
        next.chunk <- (i - 1) * x$chunk.size
        next.chunk.start <- snp.start + next.chunk
        snp.chunk <- ifelse(
            next.chunk.start + x$chunk.size > snp.end,
            snp.end - next.chunk.start + 1,
            x$chunk.size
        )
        chunk.idx <- (next.chunk + 1):(next.chunk + snp.chunk)
        
        # get genotypes for chunk
        genotypes <- getGenotype(
            genoData,
            snp = c(next.chunk.start, snp.chunk),
            scan = c(1, -1),
            drop = FALSE
        )
        
        # get the snp info file
        snp <- getAnnotation(getSnpAnnotation(genoData))[chunk.idx, ]
        
        
        colnames(snp) <- x$snp.cols
        snp <- snp[, x$snp.ord]
        
        # grab sample file data
        scanAnn <- getAnnotation(getScanAnnotation(genoData))
        
        
        listSNPGenotype <- processSNPGenotypes(
            x,
            snp = snp,
            genotypes = genotypes,
            scanAnn = scanAnn,
            exclude.snps = x$exclude.snps,
            cox.params = cox.params,
            verbose = x$verbose
        )
        
        snp <- listSNPGenotype$snp
        genotypes <- listSNPGenotype$genotypes
        
        if (x$flip.dosage)
            genotypes <- 2 - genotypes
        ########################################################################
        
        ###############################################################
        ##### SNP filtering ###########################################
        #
        # Compute allele frequency + sample MAF for EVERY SNP in the chunk
        # first, then build a single keep/drop mask. This fixes two bugs:
        #  (1) the old code subset `snp`/`genotypes` BEFORE deriving snp.drop,
        #      so dropped rows were taken from the already-subset frame with a
        #      wrong-length mask (garbage / lost SNPs);
        #  (2) it stamped bogus exp_freq_A1=1 / SAMP_MAF=0 onto SD-dropped SNPs
        #      and rbound mismatched schemas, corrupting .snps_removed (#36).
        # Doing it once also removes the redundant double-subsetting (speed).

        # allele frequency of A1 and sample minor allele frequency
        snp$exp_freq_A1 <-
            round(rowMeans2(genotypes, na.rm = TRUE) * 0.5, 4)
        snp$SAMP_MAF <- ifelse(snp$exp_freq_A1 > 0.5,
                               1 - snp$exp_freq_A1,
                               snp$exp_freq_A1)

        # Drop near-monomorphic SNPs (SD < 1e-4). For perspective: a sample of
        # 100,000,000 with a single carrier has SD = 1e-4.
        #   x <- c(rep(0, 1e8), 1); sd(x)
        ok.sd <- rowSds(genotypes, na.rm = TRUE) > 1e-4

        # Drop SNPs failing the user MAF threshold (if any). SAMP_MAF is already
        # min(af, 1-af), so a single one-sided test is correct.
        ok.maf <- if (!is.null(x$maf.filter)) {
            snp$SAMP_MAF > x$maf.filter
        } else {
            rep(TRUE, nrow(snp))
        }

        keep <- ok.sd & ok.maf
        # snp.drop carries the SAME columns as kept SNPs (incl. real freq/MAF).
        snp.drop <- snp[!keep, , drop = FALSE]
        snp <- snp[keep, , drop = FALSE]
        genotypes <- genotypes[keep, , drop = FALSE]

        if (nrow(snp.drop) > 0) {
            write.table(
                snp.drop,
                paste0(x$out.file, ".snps_removed"),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE,
                append = TRUE
            )
            snp.drop.n <- snp.drop.n + nrow(snp.drop)
        }
        
        if (nrow(genotypes) > 0) {
            # fit models in parallel
            cox.out <-
                getGenotypesCoxOut(x$inter.term,
                                   genotypes,
                                   cl,
                                   cox.params,
                                   x$print.covs)
            
            res <- coxExtract(cox.out,
                              snp,
                              x$print.covs)
            
            write.table(
                res,
                file = paste0(x$out.file, ".coxph"),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                append = TRUE
            )
            
            snp.n <- nrow(genotypes) + snp.n
            
        }
        
        
    }
    
    return(list(snp.drop.n = snp.drop.n,
                snp.n = snp.n))
    
}

replaceFileExt <- function(file.path, ext) {
    sub("\\.[^.]*?$", ext, file.path)
}