uniseq <- function(inname = "", rawname = "quaseq", rawdir = "", outdir = "", 
                      outnameN = "uniqueN_", outnameP = "uniqueP_", checkNNK = 1, indir = "", os = "win", barcodes = NULL) {
  
  # if (length(list(...)) > 0) {
  #   args <- as.list(match.call(expand.dots = TRUE))[-c(1:2)]
  #   names(args) <- NULL
  # } else {
  #   args <- NULL
  # }
  
  if (inname == "") {
    inname <- file.choose()
  }
  
  if (indir == "") {
    indir <- getwd()
    cat("using working R directory as current directory; Make sure your files are in the working R directory\n")
  }
  
  if (rawdir == "") {
    rawdir <- file.path(indir, paste0(substr(inname, 1, nchar(inname) - 4), "_QUA"))
  }
  
  if (outdir == "") {
    outdir <- file.path(indir, paste0(substr(inname, 1, nchar(inname) - 4), "_UNI"))
  }
  
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  
  ### find all files input files and cycle through all
  suff <- c("QF", "QR", "EF", "ER")
  
  
  plan(multisession)
  
  foreach (i = 1:length(barcodes)) %dofuture% {
    files <- list.files(path = rawdir, pattern = paste0(rawname,"[[:xdigit:]]{4}", barcodes[i], ".txt", collapse = ""), full.names = TRUE)
    All_A12 <- c()
    All_R36 <- c()
    iii = 0
    numLoaded <- 0
    for (currentname in files) {
      
      if(!isEmpty(readLines(currentname,1))) {
      iii = iii + 1
      
      start_time <- Sys.time()
      
      AllVar <- read.table(currentname, header = FALSE, colClasses = "character", fill = TRUE, comment.char = "", quote = "")
      R36 <- AllVar[, 1]
      R36Qual <- AllVar[, 2]
      A12 <- AllVar[, 3]
      ISNNK <- AllVar[, 4]
      rm(AllVar)
      
       if (checkNNK) {
         # eliminate the lines that do not satisfy NNK
         nonX <- which(grepl("[AC]", ISNNK)==FALSE)
         R36 <- R36[nonX]
         R36Qual <- R36Qual[nonX]
         A12 <- A12[nonX]
         ISNNK <- ISNNK[nonX]
       }
      
      # combine the pep and nucl seq into two  arrays
      All_A12 <- c(All_A12, A12)
      All_R36 <- c(All_R36, R36)
      
      numLoaded <- numLoaded + length(A12)
      
      time2 <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      print(paste0(currentname, "  loaded. Found ", numLoaded, " reads in ", time2, " sec"))
      }
    } # end of iii for-loop
    
    # save the data for current set of files
    # if (checkNNK){
    #   addNNK <- ""
    # } else {
    #   addNNK <- "noNNKchk"
    # }
    
    
    p_table = table(All_A12)
    # n_table = table(All_R36)
    
    
    All_A12 = names(p_table)
    FreqP = as.vector(p_table)
    
    # totalS <- length(All_A12)
    # totalN <- length(All_R36)
    
    time_start <- Sys.time()
    fhP <- file(paste0(outdir,"/", paste0(outnameP, barcodes[i],  ".txt")), "w") #addNNK,
    for (ii in 1:length(All_A12)){
      SP = strrep(" ", times = 14-nchar(All_A12[ii]))
    cat(sprintf("%s\n", paste0(All_A12[ii],SP, FreqP[ii])), file = fhP, append = TRUE)
      
      }
    
    time_end <- Sys.time()
    time_diff <- time_end - time_start
    print(paste("Wrote", length(FreqP), "unique seq in", time_diff, "sec"))
    
    
  }
  closeAllConnections()
      rm(list = ls())
}
        
                    
                    
                    
                    
    
      
      
      
      
      
        
        