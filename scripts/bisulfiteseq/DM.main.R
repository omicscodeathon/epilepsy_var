#!/usr/bin/env Rscript
# DM.R for DMRichR
# Author: Ben Laufer
# Contributors: Hyeyeon Hwang and Charles Mordaunt

## corrected, shoreten DM.R script for run with block and using imperviously, and no annotation part
## and corrected for PCA and heatmap 
# Initialize --------------------------------------------------------------
# Initialize --------------------------------------------------------------

cat("\n[DMRichR] Initializing \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")


if(length(grep("genomecenter.ucdavis.edu", .libPaths())) > 0){
 .libPaths("/share/lasallelab/programs/DMRichR/R_4.0")
 AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.0")
 ExperimentHub::setExperimentHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.0")
}

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
if(suppressPackageStartupMessages(!requireNamespace("DMRichR", quietly = TRUE))){
  Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)
  BiocManager::install("ben-laufer/DMRichR")
}

# Check if optparse library is installed
if (!requireNamespace("optparse", quietly = TRUE)) {
  # If not installed, install it using remotes
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  # Install optparse
  remotes::install_github("trevorld/optparse")
}


suppressPackageStartupMessages(library(DMRichR))
library(DMRichR)
library(optparse)
library(getopt)
library(openxlsx)

#processBismark

processBismark <- function(files = list.files(path = getwd(), pattern = "*.CpG_report.txt.gz"),
                      meta = openxlsx::read.xlsx("sample_info.xlsx", colNames = TRUE) %>% dplyr::mutate_if(is.character, as.factor),
                      testCovariate = testCovariate,
                      adjustCovariate = NULL,
                      matchCovariate = NULL,
                      coverage = coverage,
                      cores = cores,
                      perGroup = perGroup,
                      sexCheck = FALSE){
  
  cat("\n[DMRichR] Processing Bismark cytosine reports \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  start_time <- Sys.time()
  print(glue::glue("Selecting files..."))
  files.idx <- pmatch(meta$Name, files)
# files.idx <- na.omit(files.idx)
 files <- files[files.idx]
   #names <- as.data.frame(gsub( "_.*$","", files[files.idx])) # For colData, but jumbles file order with parallel processing
  #colnames(names) <- "Name"
  #rownames(names) <- names[,1]
  #names[,1] <- NULL
  
  # glue::glue("Determining parallelization...") # Does not work on some clusters due to use of BiocParallel, but speeds up desktops 
  # if(cores >= 4){
  #  BPPARAM <- BiocParallel::MulticoreParam(workers = floor(cores/4), progressbar = TRUE)
  #  nThread <- as.integer(floor(cores/floor(cores/4)))
  #  glue::glue("Parallel processing will be used with {floor(cores/4)} cores consisting of {nThread} threads each")
  # }else if(cores < 4){
  #  BPPARAM <- BiocParallel::MulticoreParam(workers = 1, progressbar = TRUE)
  #  nThread <- as.integer(1)
  #  glue::glue("Parallel processing will not be used")
  # }
  
  print(glue::glue("Reading cytosine reports..."))
  bs <- bsseq::read.bismark(files = files,
                            #colData = names,
                            rmZeroCov = FALSE,
                            strandCollapse = TRUE,
                            verbose = TRUE,
     #                       BPPARAM = BiocParallel::MulticoreParam(workers = cores, progressbar = FALSE), # BPPARAM # bpparam() # MulticoreParam(workers = cores, progressbar = TRUE)
                            nThread = 1) # 1L # nThread
  
  print(glue::glue("Assigning sample metadata with {testCovariate} as factor of interest..."))
  sampleNames(bs) <- gsub( "_.*$","", sampleNames(bs))
  meta <- meta[order(match(meta[,1],sampleNames(bs))),]
  stopifnot(sampleNames(bs) == as.character(meta$Name))
  pData(bs) <- cbind(pData(bs), meta[2:length(meta)])
  print(pData(bs))
  bs <- GenomeInfoDb::keepStandardChromosomes(bs, pruning.mode = "coarse")
  GenomeInfoDb::seqlevelsStyle(bs) <- "UCSC"
  
  if (sexCheck == TRUE) {
    
    # Check sex of samples using k-means clustering
    print(glue::glue("Checking sex of samples..."))
    bs.chrX <- bs[seqnames(bs) == 'chrX']
    bs.chrY <- bs[seqnames(bs) == 'chrY']
    
    coverageChrX <- bsseq::getCoverage(bs.chrX) %>%
      DelayedMatrixStats::colSums2()
    coverageChrY <- bsseq::getCoverage(bs.chrY) %>%
      DelayedMatrixStats::colSums2()
    sexCluster <- kmeans(coverageChrY / coverageChrX, centers = 2)
    allSameFlag <- "No"
    # If the value of one center is greater than 2x the value of the other
    if (max(sexCluster$centers) / min(sexCluster$centers) > 2) {
      maleIdx <- which(sexCluster$centers == max(sexCluster$centers))
      predictedSex <- character()
      for (idx in sexCluster$cluster) {
        if (idx == maleIdx) {
          predictedSex <- c(predictedSex, "M")
        } else {
          predictedSex <- c(predictedSex, "F")
        }
      }
    } else {
      allSameFlag <- "Yes"
      # Samples are either all male or all female
      predictedSex <- rep("all Male or all Female", length(sexCluster$cluster))
    }
    # Check for mismatch between predicted sex and sample info sex
    
    sampleInfo <- bs %>% pData()
    sampleInfo$Sex <- sampleInfo$Sex %>% as.character()
    
    sexMismatch <- character()
    mismatchSamples <- character()
    for (i in 1:length(sexCluster$cluster)) {
      if (sampleInfo$Sex[i] %in% c("Male", "male", "M", "m")) {
        sampleInfo$Sex[i] = "M"
      } else if (sampleInfo$Sex[i] %in% c("Female", "female", "F", "f")) {
        sampleInfo$Sex[i] = "F"
      }
      
      if (allSameFlag == "No") {
        if (predictedSex[i] == sampleInfo$Sex[i]) {
          sexMismatch <- sexMismatch %>% append(".")
        } else {
          sexMismatch <- sexMismatch %>% append("Mismatch")
          mismatchSamples <- mismatchSamples %>% append(sampleInfo %>% rownames() %>% .[i])
        }
      }
    }
    
    if (allSameFlag == "No") {
      if (length(mismatchSamples) == 0) {
        print(glue::glue("Sex of all samples matched correctly. Sex choromosomes will now be dropped"))
        bs <- GenomeInfoDb::dropSeqlevels(bs,
                                          c("chrX", "chrY"),
                                          pruning.mode = "coarse")
      } else {
        stop("Sex mismatched for the following ", toString(length(mismatchSamples)), " sample(s): ", toString(mismatchSamples), ". Rerun after correcting sample info file.")
      }
    } else {
      # allSameFlag == "Yes"
      if (length(unique(sampleInfo$Sex)) == 1) {
        print(glue::glue("Sex of samples match correctly as all male or all female."))
      } else {
        stop("Sex of samples predicted to be all male or all female. Sample info file is inconsistent with prediction. Rerun after correcting sample info file.")
      }
    }
    
    #    sexCheckResult <- data.frame(
    #      "Sample_name" = sampleInfo %>% rownames(),
    #      "ChrX_coverage" = coverageChrX,
    #      "ChrY_coverage" = coverageChrY,
    #      "ChrY_ChrX_ratio" = (coverageChrY / coverageChrX),
    #      "ChrY_ChrX_percent" = (coverageChrY / coverageChrX) * 100,
    #      "Predicted_sex" = predictedSex,
    #      "Sample_info_sex" = bs %>% pData() %>% .$Sex,
    #      "Sex_mismatch" = sexMismatch
    #    )
    #    save(sexCheckResult, file = "sexCheckResult.RData")
  }
  
  
  print(glue::glue("Filtering CpGs for {testCovariate}..."))
  pData(bs)[[testCovariate]] <- as.factor(pData(bs)[[testCovariate]])
  loci.cov <- bsseq::getCoverage(bs, type = "Cov")
  
  if(!is.null(adjustCovariate)){
    excludeCovar <- NULL
    for(i in 1:length(adjustCovariate)){
      if(is.numeric(pData(bs)[, adjustCovariate[i]])){
        print(glue::glue("Assuming adjustment covariate {adjustCovariate[i]} is continuous and excluding it from filtering..."))
        excludeCovar <- c(excludeCovar, adjustCovariate[i])
        
      }else{
        print(glue::glue("Assuming adjustment covariate {adjustCovariate[i]} is discrete and including it for filtering..."))
      }
    }
    adjustCovariate <- adjustCovariate[!adjustCovariate %in% excludeCovar]
  }
  
  if(!is.null(matchCovariate)){
    if(length(matchCovariate) > 1){
      stop(print(glue::glue("Only one matching covariate can be used")))
      
    }else if(is.numeric(pData(bs)[, matchCovariate])){
      stop(print(glue::glue("Matching covariate {matchCovariate} must be discrete")))
      
    }else{
      print(glue::glue("Assuming matching covariate {matchCovariate} is discrete and including it for filtering..."))
    }
  }
  
  covar.groups <- apply(pData(bs)[, as.character(c(testCovariate, adjustCovariate, matchCovariate))] %>% as.data.frame(), 
                        MARGIN = 1, FUN = paste, collapse = "_") %>% 
    as.factor() # Covariate combination groups
  
  group.samples <- split(t(loci.cov >= coverage) %>% as.data.frame(), f = covar.groups) %>% 
    parallel::mclapply(FUN = as.matrix, mc.cores = cores) %>% 
    parallel::mclapply(FUN = DelayedMatrixStats::colSums2, mc.cores = cores) %>%
    simplify2array() %>%
    as.data.frame() # Samples in each cov.group meeting coverage threshold by CpG (slow)
  
  print(glue::glue("Making coverage filter table..."))
  perGroup.seq <- seq(0,1,0.05)
  covFilter <- NULL
  for(i in 1:length(perGroup.seq)){
    groups.n <- (table(covar.groups) * perGroup.seq[i]) %>% ceiling(.) %>% as.integer()
    perGroup.seq.test <- mapply(function(x, y){x >= y}, 
                                x = group.samples, 
                                y = (table(covar.groups) * perGroup.seq[i]) %>%
                                  ceiling(.) %>%
                                  as.integer()) # Test if enough samples are in each group by CpG
    CpGs <- sum(DelayedMatrixStats::rowSums2(perGroup.seq.test) >= length(unique(covar.groups))) # Total CpGs meeting coverage threshold in at least perGroup of all covariate combos
    temp <- c(perGroup.seq[i] * 100, groups.n, CpGs, round(CpGs * 100 / length(bs), 2))
    covFilter <- rbind(covFilter, temp)
  }
  covFilter <- as.data.frame(covFilter, row.names = 1:nrow(covFilter))
  colnames(covFilter) <- c("perGroup", paste("n", levels(covar.groups), sep = "_"), "nCpG", "perCpG")
  print(covFilter)
  
  if(perGroup == 1){
    print(glue::glue("Filtering for {coverage}x coverage in all samples"))
    sample.idx <- which(pData(bs)[[testCovariate]] %in% levels(pData(bs)[[testCovariate]]))
    loci.idx <- which(DelayedMatrixStats::rowSums2(bsseq::getCoverage(bs, type = "Cov") >= coverage) >= length(sample.idx))
    bs.filtered <- bs[loci.idx, sample.idx]
    
  }else if(perGroup < 1){
    print(glue::glue("Filtering for {coverage}x coverage in at least {perGroup*100}% of samples for \\
                                 all combinations of covariates..."))
    sample.idx <- which(pData(bs)[[testCovariate]] %in% levels(pData(bs)[[testCovariate]]))
    perGroup.test <- mapply(function(x, y){x >= y}, 
                            x = group.samples, 
                            y = (table(covar.groups) * perGroup) %>% ceiling(.) %>% as.integer()) # Test if enough samples are in each group by CpG
    loci.idx <- which(DelayedMatrixStats::rowSums2(perGroup.test) >= length(unique(covar.groups))) # Which CpGs meet coverage threshold in at least perGroup of all covariate combos
    bs.filtered <- bs[loci.idx, sample.idx]
    
  }else if(perGroup > 1){
    stop(print(glue::glue("perGroup is {perGroup} and cannot be greater than 1, which is 100% of samples")))
    
  }else{
    stop(print(glue::glue("processBismark arguments")))
  } 
  
  glue::glue("Assigning colors for plotting...")
  pData <- pData(bs.filtered)
  if(length(levels(pData[,testCovariate])) == 2){
    pData$col <- NULL
    pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "mediumblue"
    pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <- "firebrick3"
    pData(bs.filtered) <- pData
  }
  
  print(glue::glue("processBismark timing..."))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  print(glue::glue("Before filtering for {coverage}x coverage there were {nrow(bs)} CpGs, \\
                         after filtering there are {nrow(bs.filtered)} CpGs, \\
                         which is {round(nrow(bs.filtered)/nrow(bs)*100,1)}% of all CpGs."))
  
  return(bs.filtered)
}

bispro <- get("processBismark", envir = asNamespace('DMRichR'))
environment(processBismark) <- environment(bispro)
attributes(processBismark) <- attributes(bispro)
assignInNamespace("processBismark", processBismark, ns="DMRichR")

#windows 

windows <- function(bs.filtered.bsseq = bs.filtered.bsseq,
                    size = 5e4,
                    goi = goi){
  print(glue::glue("Obtaining {size/1000} Kb window individual smoothed methylation values from the {BSgenome::commonName(goi)} genome"))
  goi %>%
    GenomeInfoDb::seqlengths() %>%
    GenomicRanges::tileGenome(tilewidth = size,
                              cut.last.tile.in.chrom = TRUE) %>%
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
    bsseq::getMeth(BSseq = bs.filtered.bsseq,
                   regions = .,
                   type = "smooth",
                   what = "perRegion") %>% 
    na.omit() %>%
    return()
}

win <- get("windows", envir = asNamespace('DMRichR'))
environment(windows) <- environment(win)
attributes(windows) <- attributes(win)
assignInNamespace("windows", windows, ns="DMRichR")

#heatmap
smoothPheatmap <- function(bs.filtered.bsseq = bs.filtered.bsseq,
                           sigRegions = sigRegions,
                           testCovariate = testCovariate,
                           ...){
  cat("\n[DMRichR] DMR heatmap \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  bsseq::getMeth(BSseq = bs.filtered.bsseq,
                 regions = sigRegions,
                 type = "smooth",
                 what = "perRegion") %>%
    na.omit() %>%
    as.matrix() %>%
    pheatmap::pheatmap(.,
                       scale = "row",
                       annotation_col =  bs.filtered.bsseq %>%
                         pData %>%
                         as.data.frame() %>%
                         dplyr::select_if(~ nlevels(.) > 1),
                       color = RColorBrewer::brewer.pal(11,
                                                        name = "RdBu") %>%
                         rev(),
                       show_colnames = FALSE,
                       #angle_col = 45,
                       border_color = "grey",
                       main = glue::glue("Z-Scores of {length(sigRegions)} Differentially Methylated Regions"),
                       fontsize = 10,
                       filename = "./DMRs/heatmap.pdf",
                       width = 10,
                       height = 14,
                       annotation_colors = DMRichR::gg_color_hue(2) %>%
                         setNames(bs.filtered.bsseq %>%
                                    pData() %>%
                                    as.data.frame() %>%
                                    purrr::pluck(testCovariate) %>%
                                    unique() %>%
                                    sort() %>%
                                    rev()) %>%
                         list(testCovariate = .) %>%
                         setNames(testCovariate),
                       ...
    ) %>%
    return()
}

htmap <- get("smoothPheatmap", envir = asNamespace('DMRichR'))
environment(smoothPheatmap) <- environment(htmap)
attributes(smoothPheatmap) <- attributes(htmap)
assignInNamespace("smoothPheatmap", smoothPheatmap, ns="DMRichR")

#PCA
PCA <- function(matrix = matrix,
                testCovariate = testCovariate,
                bs.filtered.bsseq = bs.filtered.bsseq){
  
  print(glue::glue("PCA of {length(matrix)} sites"))
  
  matrix %>%
    PCAtools::pca(scale = TRUE,
                  metadata = pData(bs.filtered.bsseq),
                  removeVar = 0.1) %>%
    PCAtools::biplot(colby = testCovariate,
                     colkey = DMRichR::gg_color_hue(2) %>%
                       setNames(bs.filtered.bsseq %>%
                                  pData() %>%
                                  as.data.frame() %>%
                                  purrr::pluck(testCovariate) %>%
                                  unique() %>%
                                  sort() %>%
                                  rev()),
                     pointSize = 6,
                     labSize = 4,
                     legendPosition = 'top',
                     legendLabSize = 16,
                     legendIconSize = 8.0)
}

pc <- get("PCA", envir = asNamespace('DMRichR'))
environment(PCA) <- environment(pc)
attributes(PCA) <- attributes(pc)
assignInNamespace("PCA", PCA, ns="DMRichR")




#DM.R--------------------

DM.R <- function(genome = c("hg38", "hg19", "mm10", "mm9", "rheMac10",
			                                "rheMac8", "rn6", "danRer11", "galGal6",
							                            "bosTau9", "panTro6", "dm6", "susScr11",
							                            "canFam3", "TAIR10", "TAIR9"),
		                  coverage = 1,
				                   perGroup = 0.05,
				                   minCpGs = 5,
						                    maxPerms = 8,
						                    maxBlockPerms = 8,
								                     cutoff = 0.025,
								                     testCovariate = testCovariate,
										                      adjustCovariate = NULL,
										                      matchCovariate = NULL,
												                       cores = 20,
												                       GOfuncR = TRUE,
														                        sexCheck = TRUE,
														                        EnsDb = FALSE,
																	                 cellComposition = FALSE){
	  
	  
	  # Check dmrseq version 
	  if(Biobase::package.version("dmrseq") %>%
	          stringr::str_remove("1.") %>%
		       as.numeric() < 7.3){
		      warning(paste("Your version of dmrseq is out of date and contains a bug.",
				                      "This bug won't affect the DMRichR run but could affect your custom follow up analyses.",
						                        "See the install section of the DMRichR README for the code to manually update.",
						                        "Read more about the issue: https://github.com/kdkorthauer/dmrseq/issues/37"))
	    }
  
  # Set options
  options(scipen = 999)
    options(readr.num_columns = 0)
    
    # Check for requirements
    stopifnot(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10",
			                              "rheMac8", "rn6", "danRer11", "galGal6",
						                                "bosTau9", "panTro6", "dm6", "susScr11",
						                                "canFam3", "TAIR10", "TAIR9"))
      stopifnot(!is.null(testCovariate))
      stopifnot(coverage >= 1)
        
        # Check for more permutations than samples
        nSamples <- openxlsx::read.xlsx("sample_info.xlsx", colNames = TRUE) %>%
		    nrow()
	      
	      if(nSamples < maxPerms){
		          print(glue::glue("Warning: You have requested {maxPerms} permutations for the DMR analysis, \\
					                      which is more than the {nSamples} samples you have. \\
							                         maxPerms will now be changed to {nSamples}."))
							          maxPerms <- nSamples
								    }
	        
	        if(nSamples < maxBlockPerms){
			    print(glue::glue("Warning: You have requested {maxBlockPerms} permutations for the block analysis, \\
					                        which is more than the {nSamples} samples you have. \\
								                   maxBlockPerms will now be changed to {nSamples}."))
								    maxBlockPerms <- nSamples
								      }
	        
	        rm(nSamples)
		  
		  # Print
		  print(glue::glue("genome = {genome}"))
		  print(glue::glue("coverage = {coverage}"))
		    print(glue::glue("perGroup = {perGroup}"))
		    print(glue::glue("minCpGs = {minCpGs}"))
		      print(glue::glue("maxPerms = {maxPerms}"))
		      print(glue::glue("maxBlockPerms = {maxBlockPerms}"))
		        print(glue::glue("cutoff = {cutoff}"))
		        print(glue::glue("testCovariate = {testCovariate}"))
			  print(glue::glue("adjustCovariate = {adjustCovariate}"))
			  print(glue::glue("matchCovariate = {matchCovariate}"))
			    print(glue::glue("cores = {cores}"))
			    print(glue::glue("cellComposition = {cellComposition}"))
			      print(glue::glue("sexCheck = {sexCheck}"))
			      print(glue::glue("EnsDb = {EnsDb}"))
			        print(glue::glue("GOfuncR = {GOfuncR}"))
			        
			        # Setup annotation databases ----------------------------------------------
			        
			        cat("\n[DMRichR] Selecting annotation databases \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
				  
				  DMRichR::annotationDatabases(genome = genome,
						                                      EnsDb = EnsDb)
				  
				  print(glue::glue("Saving Rdata..."))
				    dir.create("RData")
				   settings_env <- ls(all = TRUE)
				     save(list = settings_env, file = "RData/settings.RData")
				   #   load("RData/settings.RData")
				      
				      # Load and process samples ------------------------------------------------
				      
				     bs.filtered <- DMRichR::processBismark(files = list.files(path = getwd(),
				                                                               pattern = "*.CpG_report.txt.gz"),
				                                           meta = openxlsx::read.xlsx("sample_info.xlsx", colNames = TRUE) %>%
														                                                 dplyr::mutate_if(is.character, as.factor),
																		                                          testCovariate = testCovariate,
																		                                          adjustCovariate = adjustCovariate,
															                                          matchCovariate = matchCovariate,
																						                                           coverage = coverage,
																						                                           cores = cores,
																													                                            perGroup = perGroup,
																												                                            sexCheck = sexCheck)
				 
				      
				      print(glue::glue("Saving Rdata..."))
				      save(bs.filtered, file = "RData/bismark.RData")
			#		 load("RData/bismark.RData")
				  
				  print(glue::glue("Building annotations for plotting..."))
				 if(is(TxDb, "TxDb")){
					      annoTrack <- dmrseq::getAnnot(genome)
				   }else if(is(TxDb, "EnsDb")){
					        annoTrack <- GenomicRanges::GRangesList(CpGs = DMRichR::getCpGs(genome),
												                                            Exons = DMRichR::getExons(TxDb),
																	                                                compress = FALSE)
					      }
					    
				#	    # Background --------------------------------------------------------------
					    
					    cat("\n[DMRichR] Getting bsseq background regions \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
					    dir.create("Extra")
					      
					      DMRichR::getBackground(bs.filtered,
								                              minNumRegion = minCpGs,
										                               maxGap = 1000) %>% 
					        write.table(file = "Extra/bsseq_background.csv",
							                    sep = ",",
									                    quote = FALSE,
									                    row.names = FALSE)
						  
						  
					      
					      # Blocks ------------------------------------------------------------------
					      
					      cat("\n[DMRichR] Testing for blocks with dmrseq \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
					      start_time <- Sys.time()
					     
					     tryCatch({
					       blocks <- dmrseq::dmrseq(bs = bs.filtered,
					                                 cutoff = cutoff,
					                                 maxPerms = maxBlockPerms,
					                                 testCovariate = testCovariate,
					                                 adjustCovariate = adjustCovariate,
					                                 matchCovariate = matchCovariate,
					                                 block = TRUE,
					                                 minInSpan = 500,
					                                 bpSpan = 5e4,
					                                 maxGapSmooth = 1e6,
					                                 maxGap = 5e3,
					                                 minNumRegion = (minCpGs*2)
					                     #            BPPARAM = BiocParallel::MulticoreParam(workers = cores)
					        )
					        
					        print(glue::glue("Selecting significant blocks..."))
					        
					        if(length(blocks) != 0){
					          blocks <- blocks %>% 
					            plyranges::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
					                                                           stat < 0 ~ "Hypomethylated"),
					                              difference = round(beta/pi *100))
					        }
					    
						if(sum(blocks$qval < 0.05, na.rm = TRUE) == 0 & sum(blocks$pval < 0.05, na.rm = TRUE) != 0) {
					          sigBlocks <- blocks %>%
					            plyranges::filter(pval < 0.05)
					        }else if(sum(blocks$qval < 0.05) >= 1){
					          sigBlocks <- blocks %>%
					            plyranges::filter(qval < 0.05)
					        }else if(sum(blocks$pval < 0.05) == 0 & length(blocks) != 0){
					          glue::glue("No significant blocks detected in {length(blocks)} background blocks")
					        }else if(length(blocks) == 0){
					          glue::glue("No background blocks detected")
					        }
					        
					        if(length(blocks) != 0){
					          print(glue::glue("Exporting block and background information..."))
					          
					          dir.create("Blocks")
					         gr2bed(blocks, "Blocks/backgroundBlocks.bed")
					          if(sum(blocks$pval < 0.05) > 0){
					            print(glue::glue("{length(sigBlocks)} significant blocks of differential methylation in {length(blocks)} background blocks"))
					            gr2bed(sigBlocks, "Blocks/blocks.bed")
					            
					            print(glue::glue("Annotating and plotting blocks..."))
					            pdf("Blocks/Blocks.pdf", height = 7.50, width = 11.50)
					            dmrseq::plotDMRs(bs.filtered,
					                             regions = blocks,
					                             testCovariate = testCovariate,
					                             annoTrack = annoTrack,
					                             regionCol = "#FF00001A",
					                             qval = FALSE,
					                             stat = FALSE)
					            dev.off()
					          }
					        }
					        
					        print(glue::glue("Blocks timing..."))
					        end_time <- Sys.time()
					   	end_time - start_time
					        
					        print(glue::glue("Saving RData..."))
					        save(blocks, file = "RData/Blocks.RData")
					        load("RData/Blocks.RData")
					        
					      },
					      error = function(error_condition) {
					        print(glue::glue("Warning: Block analysis has produced an error"))
					      })
						
						  
					      # DMRs --------------------------------------------------------------------
					      
					      cat("\n[DMRichR] Testing for DMRs with dmrseq \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
					      start_time <- Sys.time()
					      
					      regions <- dmrseq::dmrseq(bs = bs.filtered,
					                                cutoff = cutoff,
					                                minNumRegion = minCpGs,
					                                maxPerms = maxPerms,
					                                testCovariate = testCovariate,
					                                adjustCovariate = adjustCovariate,
					                                matchCovariate = matchCovariate,
					                                BPPARAM = BiocParallel::MulticoreParam(workers = cores)
					      )
					      
					      print(glue::glue("Selecting significant DMRs..."))
					      
					      regions <- regions %>% 
					        plyranges::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
					                                                       stat < 0 ~ "Hypomethylated"),
					                          difference = round(beta/pi *100))
					      
					      if(sum(regions$qval < 0.05, na.rm = TRUE) < 100 & sum(regions$pval < 0.05, na.rm = TRUE) != 0){
					        sigRegions <- regions %>%
					          plyranges::filter(pval < 0.05)
		          		      }else if(sum(regions$qval < 0.05) >= 100){
					        sigRegions <- regions %>%
					          plyranges::filter(qval < 0.05)
				      }else if(sum(regions$pval < 0.05) == 0){
					        stop(glue::glue("No significant DMRs detected in {length(regions)} background regions"))
					      }
					      
					      
				      print(glue::glue("Exporting DMR and background region information..."))
				      dir.create("DMRs")
				      gr2bed(sigRegions, "DMRs/DMRs.bed")
				      gr2bed(regions, "DMRs/backgroundRegions.bed")
					      
					      if(sum(sigRegions$stat > 0) > 0 & sum(sigRegions$stat < 0) > 0){
					        
					        print(glue::glue("Summary: There are {tidySigRegions} DMRs \\
						                  ({tidyHyper}% hypermethylated, {tidyHypo}% hypomethylated) \\
												               from {tidyRegions} background regions consisting of {tidyCpGs} CpGs \\
													               assayed at {coverage}x coverage", 
				                         tidySigRegions = length(sigRegions),
				                         tidyHyper = round(sum(sigRegions$stat > 0) / length(sigRegions), digits = 2)*100,
				                         tidyHypo = round(sum(sigRegions$stat < 0) / length(sigRegions), digits = 2)*100,
					                         tidyRegions = length(regions),
					                         tidyCpGs = nrow(bs.filtered)))
					      }
					      
					      print(glue::glue("DMR timing..."))
					      end_time <- Sys.time()
					      end_time - start_time
					      
					      print(glue::glue("Saving Rdata..."))
					      save(regions, sigRegions, file = "RData/DMRs.RData")
					      #load("RData/DMRs.RData")
#					      
					      print(glue::glue("Annotating DMRs and plotting..."))
				      
					      pdf("DMRs/DMRs.pdf", height = 4, width = 8)
					      tryCatch({
					        DMRichR::plotDMRs2(bs.filtered,
					                           regions = sigRegions,
					                           testCovariate = testCovariate,
					                           extend = (end(sigRegions) - start(sigRegions) + 1)*2,
					                           addRegions = sigRegions,
					                           annoTrack = annoTrack,
					                           regionCol = "#FF00001A",
					                           lwd = 2,
					                           qval = FALSE,
					                           stat = FALSE,
					                           horizLegend = FALSE)
					      },
					      error = function(error_condition) {
					        print(glue::glue("Warning: One (or more) of your DMRs can't be plotted, \\
														                         try again later by manually loading R Data and subsetting sigRegions"))
					      })
					      dev.off()
					      
					      # Annotate DMRs with gene symbols -----------------------------------------
					      
						    cat("\n[DMRichR] Annotating DMRs with gene symbols \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  sigRegions %>%
    DMRichR::annotateRegions(TxDb = TxDb,
                             annoDb = annoDb) %T>%
    DMRichR::DMReport(regions = regions,
                      bs.filtered = bs.filtered,
                      coverage = coverage,
                      name = "DMReport") %>% 
    openxlsx::write.xlsx(file = "DMRs/DMRs_annotated.xlsx")
  
  print(glue::glue("Annotating background regions with gene symbols..."))
  regions %>%
    DMRichR::annotateRegions(TxDb = TxDb,
                             annoDb = annoDb) %>% 
    openxlsx::write.xlsx(file = "DMRs/background_annotated.xlsx")
  
	
					  # Individual smoothed values ----------------------------------------------
											  
											  cat("\n[DMRichR] Smoothing individual methylation values \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
											    start_time <- Sys.time()
											    
											    bs.filtered.bsseq <- bsseq::BSmooth(bs.filtered,
																                                      BPPARAM = BiocParallel::MulticoreParam(workers = cores,
																									                                                                                  progressbar = TRUE))
											    # Drop chrY in Rat only due to poor quality (some CpGs in females map to Y)
											    if(genome == "rn6"){
											      bs.filtered.bsseq <- GenomeInfoDb::dropSeqlevels(bs.filtered.bsseq,
											                                                       "chrY",
											                                                       pruning.mode = "coarse")
											      GenomeInfoDb::seqlevels(bs.filtered.bsseq)
											    }
											    
											    bs.filtered.bsseq
											        
											        print(glue::glue("Extracting individual smoothed methylation values of DMRs..."))
											        bs.filtered.bsseq %>%
													    DMRichR::smooth2txt(regions = sigRegions,
																                        txt = "DMRs/DMR_individual_smoothed_methylation.txt")
												  
												  print(glue::glue("Extracting individual smoothed methylation values of background regions..."))
												    bs.filtered.bsseq %>%
													        DMRichR::smooth2txt(regions = regions,
																                            txt = "DMRs/background_region_individual_smoothed_methylation.txt")
												    
												    print(glue::glue("Individual smoothing timing..."))
												      end_time <- Sys.time()
												      end_time - start_time
												        
												        print(glue::glue("Saving Rdata..."))
												        save(bs.filtered.bsseq,
													            file = "RData/bsseq.RData")
													  #load("RData/bsseq.RData")
													  
													  # ChromHMM and Reference Epigenomes ---------------------------------------
													  
													  if(length(grep("genomecenter.ucdavis.edu", .libPaths())) > 0 & genome == "hg38"){
														      
														      dir.create("LOLA")
													      setwd("LOLA")
													          
													          dmrList <- sigRegions %>% 
															        DMRichR::dmrList()
															    
															    LOLA <- function(x){
																          
																          dir.create(names(dmrList)[x])
															          setwd(names(dmrList)[x])
																        
																        dmrList[x] %>%
																		        DMRichR::chromHMM(regions = regions,
																					                            cores = floor(cores/3)) %>% 
																          DMRichR::chromHMM_heatmap()
																        
																        dmrList[x] %>%
																		        DMRichR::roadmap(regions = regions,
																					                          cores = floor(cores/3)) %>% 
																	        DMRichR::roadmap_heatmap()
																	      
																	      if(file.exists("Rplots.pdf")){file.remove("Rplots.pdf")}
																	          }
															        
															        parallel::mclapply(seq_along(dmrList),
																		                          LOLA,
																					                         mc.cores = 3,
																					                         mc.silent = TRUE)
															        
															        setwd("..")
																  }
													  
													  # HOMER -------------------------------------------------------------------
													  
													  sigRegions %>% 
														      DMRichR::prepareHOMER(regions = regions)
													        
													        DMRichR::HOMER(genome = genome,
															                        cores = cores)
														  
														  # Smoothed global, chromosomal, and CGi methylation statistics ------------
														  
														  dir.create("Global")
														  
														  bs.filtered.bsseq %>%
															      DMRichR::globalStats(genome = genome,
																		                            testCovariate = testCovariate,
																					                             adjustCovariate = adjustCovariate,
																					                             matchCovariate = matchCovariate) %>%
														      openxlsx::write.xlsx("Global/smoothed_globalStats.xlsx") 
													        
													        # Global plots ------------------------------------------------------------
														  
														  # Heatmap -----------------------------------------------------------------
														  
														  sigRegions %>%
														    DMRichR::smoothPheatmap(bs.filtered.bsseq = bs.filtered.bsseq,
														                            testCovariate = testCovariate)
														  # Global plots ------------------------------------------------------------
													        
													        windows <- bs.filtered.bsseq %>%
															    DMRichR::windows(goi = goi)
														      
														      CpGs <- bs.filtered.bsseq %>%
															          DMRichR::CpGs()
															    
															    plots <- c("windows", "CpGs")
															      
															      if(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6",
																	                          "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")){
																          
																          CGi <- bs.filtered.bsseq %>% 
																		        DMRichR::CGi(genome = genome)
																		    
																		    plots <- c("windows", "CpGs", "CGi")
																		      }
															      
															      purrr::walk(plots,
																	                function(plotMatrix,
																				                        group =  bs.filtered.bsseq %>%
																								                         pData() %>%
																											                          dplyr::as_tibble() %>%
																														                           dplyr::pull(!!testCovariate) %>%
																																	                            forcats::fct_rev()){
																				                
																				                title <- dplyr::case_when(plotMatrix == "windows" ~ "20Kb Windows",
																									                                            plotMatrix == "CpGs" ~ "Single CpG",
																														                                              plotMatrix == "CGi" ~ "CpG Island")
																			                
																			                plotMatrix %>%
																						                  get() %>% 
																								                    DMRichR::PCA(testCovariate = testCovariate,
																												                                bs.filtered.bsseq = bs.filtered.bsseq) %>%
																					                  ggplot2::ggsave(glue::glue("Global/{title} PCA.pdf"),
																									                                    plot = .,
																													                                      device = NULL,
																													                                      width = 11,
																																	                                        height = 8.5)
																							                  
																							                  plotMatrix %>%
																										                    get() %>% 
																												                      DMRichR::densityPlot(group = group) %>% 
																														                        ggplot2::ggsave(glue::glue("Global/{title} Density Plot.pdf"),
																																			                                  plot = .,
																																							                                    device = NULL,
																																							                                    width = 11,
																																											                                      height = 4)
																									                  
																									                  Glimma::glMDSPlot(plotMatrix %>%
																													                                        get(),
																																	                                  groups = cbind(bsseq::sampleNames(bs.filtered.bsseq),
																																							                                                  pData(bs.filtered.bsseq)) %>%
																																		                                    dplyr::as_tibble() %>% 
																																						                                        dplyr::select(-col) %>%
																																											                                    dplyr::rename(Name = bsseq..sampleNames.bs.filtered.bsseq.),
																																														                                      path = getwd(),
																																						                                      folder = "interactiveMDS",
																																										                                        html = glue::glue("{title} MDS plot"),
																																										                                        launch = FALSE)
																											                })
															     
																  
																  # CpG and genic enrichment testing ----------------------------------------
																  
																  cat("\n[DMRichR] Performing DMRichments \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
																    
																    DMRich <- function(x){
																	        
																	        if(genome %in% c("fHypTra1", "hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6", "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")){
																			      print(glue::glue("Running CpG annotation enrichments for {names(dmrList)[x]}"))
																        dmrList[x] %>% 
																		        DMRichR::DMRichCpG(regions = regions,
																					                              genome = genome) %T>%
																	        openxlsx::write.xlsx(file = glue::glue("DMRichments/{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
																		        DMRichR::DMRichPlot(type = "CpG") %>% 
																			        ggplot2::ggsave(glue::glue("DMRichments/{names(dmrList)[x]}_CpG_enrichments.pdf"),
																						                        plot = ., 
																									                        width = 4,
																									                        height = 3)
																		    }
																      
																      print(glue::glue("Running gene region annotation enrichments for {names(dmrList)[x]}"))
																          dmrList[x] %>% 
																		        DMRichR::DMRichGenic(regions = regions,
																					                                TxDb = TxDb,
																									                           annoDb = annoDb) %T>%
																            openxlsx::write.xlsx(file = glue::glue("DMRichments/{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
																	          DMRichR::DMRichPlot(type = "genic") %>% 
																		        ggplot2::ggsave(glue::glue("DMRichments/{names(dmrList)[x]}_genic_enrichments.pdf"),
																					                      plot = ., 
																							                            width = 4,
																							                            height = 4)
																	      }
																    
																    dmrList <- sigRegions %>% 
																	        DMRichR::dmrList()
																	  
																	  dir.create("DMRichments")
																	    
																	    purrr::walk(seq_along(dmrList),
																			              DMRich)
																	    
																	    purrr::walk(dplyr::case_when(genome %in% c("fHypTra1", "hg38", "hg19", "mm10", "mm9", "rn6") ~ c("CpG", "genic"),
																					                                TRUE ~ "genic") %>%
																	                    unique(),
																		                  function(type){
																					                  
																					                  print(glue::glue("Creating DMRichMultiPlots for {type} annotations"))
																			                    
																			                    DMRichR::DMparseR(direction =  c("All DMRs",
																									                                                      "Hypermethylated DMRs",
																															                                                       "Hypomethylated DMRs"),
																							                                        type = type) %>%
																					                      DMRichR::DMRichPlot(type = type,
																										                                        multi = TRUE) %>% 
																							                        ggplot2::ggsave(glue::glue("DMRichments/{type}_multi_plot.pdf"),
																												                                  plot = .,
																																                                    device = NULL,
																																                                    height = dplyr::case_when(type == "genic" ~ 5,
																																							                                                                  type == "CpG" ~ 3.5),
																												                                  width = 7)
																										              })
																	      
																	      # Overlap with human imprinted genes --------------------------------------
																	      
																	      cat("\n[DMRichR] Testing for imprinted gene enrichment \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
																	      
																	      dmrList <- sigRegions %>% 
																		          DMRichR::dmrList()
																		    
																		    sink("DMRs/human_imprinted_gene_overlaps.txt")
																		      
																		      purrr::walk(seq_along(dmrList),
																				                function(x){
																							                print(glue::glue("Analyzing {names(dmrList)[x]}"))
																						                
																						                dmrList[x] %>%
																									                  DMRichR::imprintOverlap(regions = regions,
																														                                            TxDb = TxDb,
																																			                                              annoDb = annoDb)
																								              })
																		      
																		      sink()
																		        
																		        # Manhattan plot ----------------------------------------------------------
																		        
																		        regions %>%
																				    DMRichR::annotateRegions(TxDb = TxDb,
																							                                  annoDb = annoDb) %>% 
																		          DMRichR::Manhattan()
																		    
																		    # Gene Ontology analyses --------------------------------------------------
																		    
																		    cat("\n[DMRichR] Performing gene ontology analyses \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
																		      
																		      dir.create("Ontologies")
																		      
																		      if(genome %in% c("fHypTra1", "hg38", "hg19", "mm10", "mm9")){
																			          
																			          print(glue::glue("Running GREAT"))
																		          GREATjob <- sigRegions %>%
																				        dplyr::as_tibble() %>%
																					      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
																					            rGREAT::submitGreatJob(bg = regions,
																									                                species = genome,
																													                             rule = "oneClosest",
																													                             request_interval = 1,
																																                                  version = "4.0.4")
																			      
																			      print(glue::glue("Saving and plotting GREAT results"))
																			          GREATjob %>%
																					        rGREAT::getEnrichmentTables(category = "GO") %T>% #%>%
																						      #purrr::map(~ dplyr::filter(., Hyper_Adjp_BH < 0.05)) %T>%
																						      openxlsx::write.xlsx(file = glue::glue("Ontologies/GREAT_results.xlsx")) %>%
																						            DMRichR::slimGO(tool = "rGREAT",
																									                          annoDb = annoDb,
																												                        plots = FALSE) %T>%
																			            openxlsx::write.xlsx(file = glue::glue("Ontologies/GREAT_slimmed_results.xlsx")) %>%
																				          DMRichR::GOplot() %>%
																					        ggplot2::ggsave(glue::glue("Ontologies/GREAT_plot.pdf"),
																								                      plot = .,
																										                            device = NULL,
																										                            height = 8.5,
																													                          width = 10)
																				        
																				        # pdf(glue::glue("Ontologies/GREAT_gene_associations_graph.pdf"),
																				        #     height = 8.5,
																				        #     width = 11)
																				        # par(mfrow = c(1, 3))
																				        # res <- rGREAT::plotRegionGeneAssociationGraphs(GREATjob)
																				        # dev.off()
																				        # write.csv(as.data.frame(res),
																				        #           file = glue::glue("Ontologies/GREATannotations.csv"),
																				        #           row.names = FALSE)
																				      }
																		        
																		        if(GOfuncR == TRUE){
																				    print(glue::glue("Running GOfuncR"))
																			    sigRegions %>% 
																				          DMRichR::GOfuncR(regions = regions,
																							                          n_randsets = 1000,
																										                         upstream = 5000,
																										                         downstream = 1000,
																													                        annoDb = annoDb,
																													                        TxDb = TxDb) %T>%
																			          openxlsx::write.xlsx(glue::glue("Ontologies/GOfuncR.xlsx")) %>% 
																				        DMRichR::slimGO(tool = "GOfuncR",
																							                      annoDb = annoDb,
																									                            plots = FALSE) %T>%
																				        openxlsx::write.xlsx(file = glue::glue("Ontologies/GOfuncR_slimmed_results.xlsx")) %>% 
																					      DMRichR::GOplot() %>% 
																					            ggplot2::ggsave(glue::glue("Ontologies/GOfuncR_plot.pdf"),
																								                          plot = .,
																											                        device = NULL,
																											                        height = 8.5,
																														                      width = 10)
																					  }
																		        
																		        if(genome != "TAIR10" & genome != "TAIR9"){
																				    tryCatch({
																					          print(glue::glue("Running enrichR"))
																						        
																						        enrichR:::.onAttach() # Needed or else "EnrichR website not responding"
																						        #dbs <- enrichR::listEnrichrDbs()
																						        dbs <- c("GO_Biological_Process_2018",
																								                "GO_Cellular_Component_2018",
																										               "GO_Molecular_Function_2018",
																										               "KEGG_2019_Human",
																											                      "Panther_2016",
																											                      "Reactome_2016",
																													                     "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")
																							      
																							      if(genome %in% c("mm10", "mm9", "rn6")){
																								              dbs %>%
																										                gsub(pattern = "Human", replacement = "Mouse")
																											      }else if(genome %in% c("danRer11", "dm6")){
																												              if(genome == "danRer11"){
																														                enrichR::setEnrichrSite("FishEnrichr")
																											              }else if(genome == "dm6"){
																													                enrichR::setEnrichrSite("FlyEnrichr")}
																											              dbs <- c("GO_Biological_Process_2018",
																													                        "GO_Cellular_Component_2018",
																																                 "GO_Molecular_Function_2018",
																																                 "KEGG_2019")
																												            }
																							      
																							      sigRegions %>%
																								              DMRichR::annotateRegions(TxDb = TxDb,
																												                                        annoDb = annoDb) %>%  
																							              dplyr::select(geneSymbol) %>%
																								              purrr::flatten() %>%
																									              enrichR::enrichr(dbs) %>% 
																										              purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %T>% # %>% 
																											              #purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05)) %T>%
																											              openxlsx::write.xlsx(file = glue::glue("Ontologies/enrichr.xlsx")) %>%
																												              DMRichR::slimGO(tool = "enrichR",
																															                              annoDb = annoDb,
																																		                              plots = FALSE) %T>%
																								              openxlsx::write.xlsx(file = glue::glue("Ontologies/enrichr_slimmed_results.xlsx")) %>% 
																									              DMRichR::GOplot() %>% 
																										              ggplot2::ggsave(glue::glue("Ontologies/enrichr_plot.pdf"),
																													                              plot = .,
																																                              device = NULL,
																																                              height = 8.5,
																																			                              width = 10)
																									            
																									          },
																										      error = function(error_condition) {
																											            print(glue::glue("Warning: enrichR did not finish. \\
																														                           The website may be down or there are internet connection issues."))
																																	       })
																										    }
																			  
																			  # Machine learning --------------------------------------------------------
																			  tryCatch({
																				      methylLearnOutput <- DMRichR::methylLearn(bs.filtered.bsseq = bs.filtered.bsseq,
																										                                              sigRegions = sigRegions,
																															                                                    testCovariate = testCovariate,
																															                                                    TxDb = TxDb,
																																					                                                  annoDb = annoDb,
																																					                                                  topPercent = 1,
																																											                                                output = "all",
																																											                                                saveHtmlReport = TRUE)
																				          
																				          if(!dir.exists("./Machine_learning")) {
																						        dir.create("./Machine_learning")
																					      } 
																				          
																				          if(length(methylLearnOutput) == 1) {
																						        openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput), 
																									                                file = "./Machine_learning/Machine_learning_output_one.xlsx") 
																					      } else {
																						            openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput$`Annotated common DMRs`,
																										                                      RF_Ranking_All_DMRs = methylLearnOutput$`RF ranking`,
																														                                      SVM_Ranking_All_DMRs = methylLearnOutput$`SVM ranking`),
																										                            file = "./Machine_learning/Machine_learning_output_all.xlsx") 
																					          }
																					      
																					      print(glue::glue("Saving RData..."))
																					      save(methylLearnOutput, file = "RData/machineLearning.RData")
																					          #load("RData/machineLearing.RData")
																					        },
																					        error = function(error_condition) {
																							    print(glue::glue("Warning: methylLearn did not finish. \\
																									                           You may have not had enough top DMRs across algrothims."))
																												     })
																						  
																						  # Cell composition --------------------------------------------------------
																						  
																						  if(cellComposition == TRUE & genome %in% c("hg38", "hg19")){
																							      
																							      bs.filtered.bsseq.cc <- bs.filtered.bsseq %>%
																								            DMRichR::prepareCC(genome = genome)
																								        
																								        rm(bs.filtered.bsseq)
																									    
																									    # Houseman method
																									    
																									    HousemanCC <- bs.filtered.bsseq.cc %>%
																										          DMRichR::Houseman()
																										      
																										      # methylCC method
																										      
																										      ccDMRs <- DMRichR::find_dmrs2(mset_train_flow_sort = "FlowSorted.Blood.EPIC",
																														                                      include_cpgs = FALSE,
																																		                                        include_dmrs = TRUE)
																										          
																										          methylCC <- bs.filtered.bsseq.cc %>%
																												        methylCC::estimatecc(include_cpgs = FALSE,
																															                                include_dmrs = TRUE,
																																			                           find_dmrs_object = ccDMRs)
																											      
																											      dir.create("Cell Composition")
																											          save(HousemanCC, methylCC, ccDMRs, file = "RData/cellComposition.RData")
																											          
																											          purrr::walk(c("Houseman", "methylCC"),
																													                      function(method = method){
																																                        if(method == "Houseman"){
																																				                    CC <- HousemanCC
																															                        }else if(method == "methylCC"){
																																			                    CC <- methylCC %>%
																																						                          methylCC::cell_counts()
																																								                    }
																															                        
																															                        CC %>% 
																																			                    "*"(100) %>%
																																					                        "/"(rowSums(.)) %>% 
																																								                    as.data.frame() %>% 
																																										                        DMRichR::CCstats(bs.filtered.bsseq.cc = bs.filtered.bsseq.cc,
																																															                                      testCovariate = testCovariate,
																																																			                                           adjustCovariate = adjustCovariate,
																																																			                                           matchCovariate = matchCovariate
																																																								                       ) %T>%
																																		                    openxlsx::write.xlsx(glue::glue("Cell Composition/{method}_stats.xlsx")) %>%
																																				                        DMRichR::CCplot(testCovariate = testCovariate,
																																									                                    adjustCovariate = adjustCovariate,
																																													                                        matchCovariate = matchCovariate
																																													                        ) %>% 
																																				                        ggplot2::ggsave(glue::glue("Cell Composition/{method}_plot.pdf"),
																																									                                    plot = .,
																																													                                        device = NULL,
																																													                                        height = 6,
																																																		                                    width = 6)
																																							                })
																												    }
																						    
																						    # End ---------------------------------------------------------------------
																						    
																						    cat("\n[DMRichR] Summary \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
																						    
																						    print(glue::glue("Summary: There were {dmrLength} DMRs that covered {sigRegionPercent} of the genome. \\
																								                        The DMRs were identified from {backgroundLength } background regions that covered {regionPercent} of the genome.
																											                   {tidyHyper} of the DMRs were hypermethylated, and {tidyHypo} were hypomethylated. \\
																											                   The methylomes consisted of {tidyCpGs} CpGs.", 
																													                      dmrLength = sigRegions %>%
																																                           length() %>%
																																			                        formatC(format = "d", big.mark = ","),
																																					                   backgroundLength = regions %>%
																																								                        length() %>%
																																											                     formatC(format = "d", big.mark = ","),
																																												                        tidyHyper = (sum(sigRegions$stat > 0) / length(sigRegions)) %>%
																																																                     scales::percent(),
																																																	                        tidyHypo = (sum(sigRegions$stat < 0) / length(sigRegions)) %>%
																																																					                     scales::percent(),
																																																						                        tidyCpGs = nrow(bs.filtered) %>%
																																																										                     formatC(format = "d", big.mark = ","),
																																																											                        genomeSize = goi %>%
																																																															                     seqinfo() %>%
																																																																	                          GenomeInfoDb::keepStandardChromosomes() %>%
																																																																				                       as.data.frame() %>%
																																																																						                            purrr::pluck("seqlengths") %>%
																																																																									                         sum(),
																																																																											                    dmrSize = sigRegions %>%
																																																																														                         dplyr::as_tibble() %>%
																																																																																	                      purrr::pluck("width") %>%
																																																																																			                           sum(),
																																																																																					                      backgroundSize = regions  %>%
																																																																																								                           dplyr::as_tibble() %>%
																																																																																											                        purrr::pluck("width") %>%
																																																																																														                     sum(),
																																																																																															                        sigRegionPercent = (dmrSize/genomeSize) %>%
																																																																																																			                     scales::percent(accuracy = 0.01),
																																																																																																				                        regionPercent = (backgroundSize/genomeSize) %>%
																																																																																																								                     scales::percent(accuracy = 0.01)
																																																																																																									       ))
																						      
																						      try(if(sum(blocks$pval < 0.05) > 0 & length(blocks) != 0){
																								      print(glue::glue("{length(sigBlocks)} significant blocks of differential methylation \\
																										                  in {length(blocks)} background blocks"))
																												    }, silent = TRUE)
																								        
																								        writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
																									  if(file.exists("Rplots.pdf")){file.remove("Rplots.pdf")}
																									  
																									  print(glue::glue("Done..."))
}





dm <- get("DM.R", envir = asNamespace('DMRichR'))
environment(DM.R) <- environment(dm)
attributes(DM.R) <- attributes(dm)
assignInNamespace("DM.R", DM.R, ns="DMRichR")

# Global variables --------------------------------------------------------

cat("\n[DMRichR] Processing arguments from script \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

option_list <- list(
  optparse::make_option(c("-g", "--genome"), type = "character", default = NULL,
                        help = "Choose a genome (fHypTra1, hg38, hg19, mm10, mm9, rheMac10, rheMac8, rn6, danRer11, galGal6, bosTau9, panTro6, dm6, susScr11, canFam3, TAIR10, or TAIR9) [required]"),
  optparse::make_option(c("-x", "--coverage"), type = "integer", default = 1,
                        help = "Choose a CpG coverage cutoff [default = %default]"),
  optparse::make_option(c("-s", "--perGroup"), type = "double", default = 1,
                        help = "Choose the percent [values from 0 to 1] of samples in all combinations of covariates meeting CpG coverage cutoff [default = %default]"),
  optparse::make_option(c("-n", "--minCpGs"), type = "integer", default = 5,
                        help = "Choose the minimum number of CpGs for a DMR [default = %default]"),
  optparse::make_option(c("-p", "--maxPerms"), type = "integer", default = 10,
                        help = "Choose the number of permutations for the DMR analysis [default = %default]"),
  optparse::make_option(c("-b", "--maxBlockPerms"), type = "integer", default = 10,
                        help = "Choose the number of permutations for the block analysis [default = %default]"),
  optparse::make_option(c("-o", "--cutoff"), type = "double", default = 0.05,
                        help = "Choose the cutoff value [from 0 to 1] for the single CpG coefficient utilized to discover testable background regions [default = %default]"),
  optparse::make_option(c("-t", "--testCovariate"), type = "character", default = NULL,
                        help = "Choose a test covariate [required]"),
  optparse::make_option(c("-a", "--adjustCovariate"), type = "character", default = NULL,
                        help = "Choose covariates to directly adjust [default = NULL]"),
  optparse::make_option(c("-m", "--matchCovariate"), type = "character", default = NULL,
                        help = "Choose covariate to balance permutations [default = NULL]"),
  optparse::make_option(c("-c", "--cores"), type = "integer", default = 20,
                        help = "Choose number of cores [default = %default]"),
  optparse::make_option(c("-e", "--cellComposition"), type = "logical", default = FALSE,
                        help = "Logical to estimate blood cell composition [default = %default]"),
  optparse::make_option(c("-k", "--sexCheck"), type = "logical", default = FALSE,
                        help = "Logical to confirm sex of each sample [default = %default]"),
  optparse::make_option(c("-d", "--EnsDb"), type = "logical", default = FALSE,
                        help = "Logical to select Ensembl transcript annotation database [default = %default]"),
  optparse::make_option(c("-f", "--GOfuncR"), type = "logical", default = TRUE,
                        help = "Logical to run GOfuncR GO analysis [default = %default]")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))


# DM.R --------------------------------------------------------------------

DMRichR::DM.R(genome = opt$genome,
              coverage = opt$coverage,
              perGroup = opt$perGroup,
              minCpGs =  opt$minCpGs,
              maxPerms =  opt$maxPerms,
              maxBlockPerms = opt$maxBlockPerms,
              cutoff = opt$cutoff,
              testCovariate = opt$testCovariate,
              adjustCovariate = if(!is.null(opt$adjustCovariate)){
                adjustCovariate <- opt$adjustCovariate %>%
                  strsplit(";") %>%
                  unlist() %>%
                  as.character()
                }else if(is.null(opt$adjustCovariate)){
                  adjustCovariate <- opt$adjustCovariate
                  },
              matchCovariate = opt$matchCovariate,
              cores = opt$cores,
              GOfuncR = opt$GOfuncR,
              sexCheck = opt$sexCheck,
              EnsDb = opt$EnsDb,
              cellComposition = opt$cellComposition)





