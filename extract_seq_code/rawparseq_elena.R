rawparseq <- function(inname = "",chunk = 250000, outdir = "", outname = "parseq", indir = "", barcodes = c(), searchF = '',  searchR = '', Lins = NULL ) {
  
  # Usage:
  # rawseq()  will prompt for input file names/directories
  # or
  # rawseq('InName', File Name, 'InDir', Input Directory, 'Chunk', N,...
  #        'OutDir', Directory Name,'OutName', FileName )
  # rawseq takes raw illumina file in FASTAQ txt format, breaks it into
  # smaller chunks of N lines.
  # From each chunk, it extracts sequence string, coordinate string, and
  # quality string. The result is saved to out_dir/out_name000N.mat,
  # which contains 'AllSeq', 'AllCoord', 'AllQual' cell arrays of strings.
  # If the InName argument is not defined, the program will open a dialog
  # box to define the input file.
  # The other arguments have default values, which will be used if the user
  # did not pass the arguments.
  # The program checks for existing directory and file name. If files already
  # exist, the program asks whether files can be over-written.
  
  # check whether indir name was defined
  if (indir == "") {
    cat("using working R directory as current directory; Make sure your files are in the working R directory\n")
  }
  
  # check whether outdir name was defined
  if (outdir == "") {
    outdir <- file.path(indir, paste0(substring(inname, 1, nchar(inname) - 4), "_PAR")) # default save directory
  }
  
  # check whether output files already exist in the directory
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  seqfile = paste(indir, inname, sep="/")
  FID <- file(seqfile,"r")
  
  i <- 0
  #ii <- 0
  beg <- substr(readLines(FID,1), 1, 7)
  #f_len = length(readLines(FID))
  f_len = 951819232 #Exp12.fastq
  #f_len = 16000000
  
  ############################################### from parseq
  # make searchF and searchR
  LF <- unlist(nchar(searchF))
  LR <- unlist(nchar(searchR))
  BARL <- 5  # length of the barcode
  
  
  
  # make 1N-containing nucleotides %%%%%%%%%%%%%%
  searchRN <- c()
  searchFN <- c()
  
  for (i in 1:LF){
    searchFN[i] <- paste0(substr(searchF, 1, i-1), "N", substr(searchF, i+1, LF))
  }
  
  for (i in 1:LR){
    searchRN[i] <- paste0(substr(searchR, 1, i-1), "N", substr(searchR, i+1, LR))
  }
  
  
  # Make nucleotides with one point mutation
  searchRM <- list()
  searchFM <- list()
  
  for (i in 1:LF){
    
    temp2F <- c()
    
    
    tempF <- c(paste0(substr(searchF, 1, i-1), "A", substr(searchF, i+1, LF)),
               paste0(substr(searchF, 1, i-1), "T", substr(searchF, i+1, LF)),
               paste0(substr(searchF, 1, i-1), "C", substr(searchF, i+1, LF)),
               paste0(substr(searchF, 1, i-1), "G", substr(searchF, i+1, LF)))
    
    indexF <- sapply(tempF, function(x) grepl(searchF, x))
    
    for (j in 1:4){
      
      if (indexF[j]==FALSE){
        temp2F[length(temp2F)+1] <- tempF[j]
      }
    }
    
    searchFM[[i]] <- temp2F
  }
  
  for (i in 1:LR){
    temp2R <- c()
    
    
    tempR <- c(paste0(substr(searchR, 1, i-1), "A", substr(searchR, i+1, LR)),
               paste0(substr(searchR, 1, i-1), "T", substr(searchR, i+1, LR)),
               paste0(substr(searchR, 1, i-1), "C", substr(searchR, i+1, LR)),
               paste0(substr(searchR, 1, i-1), "G", substr(searchR, i+1, LR)))
    
    indexR <- sapply(tempR, function(x) grepl(searchR, x))
    
    
    for (j in 1:4){
      if (indexR[j]==FALSE){
        temp2R[length(temp2R)+1] <- tempR[j]
      }
      
    }
    searchRM[[i]] <- temp2R
  }
  searchFM = unlist(searchFM)
  searchRM = unlist(searchRM)
  
  
  # Make a set of nucleotides with one deletion inside
  
  LFD <- LF - 1
  LRD <- LR - 1
  sFD = unlist(strsplit(searchF,""))
  sRD = unlist(strsplit(searchR,""))
  
  searchFD = sapply(1:LF, function(i) paste0(sFD[-i], collapse = ""))
  searchRD = sapply(1:LR, function(i) paste0(sRD[-i], collapse = ""))
  
  searchFall = c(searchF, searchFM, searchFN)
  searchRall = c(searchR, searchRM, searchRN)
  
  
  ########################################################################################333333
  
  require(foreach)
  require(future)
  require(future.apply)
  require(doFuture, quietly = TRUE)
  plan(multisession)
  options(future.globals.maxSize = 1000 * 1024^2)  # 1000MB (1GB)
  
  foreach (i = 1:ceiling(f_len/(4*chunk))) %dofuture% {
    start_time <- Sys.time()
    FID <- file(seqfile,"r")
    currL <- scan(file = FID,what = "character", skip = (i-1)*4*chunk, nlines = 4*chunk, sep ="\n")
    
    
    #if (i == 1) {
    #  beg <- substr(currL[[1]], 1, 7)
    #print(beg)
    #}
    
    
    time1 <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    start_time <- Sys.time()
    
    A1 <- grepl(beg, currL, fixed = TRUE)
    endAl <- A1[(length(A1)-3):length(A1)]
    #print(A1)
    
    
    if (sum(endAl == c(1, 0, 0, 0)) != 4) {
      print("incomplete line set")
    }
    
    AllCord <- currL[A1]
    Slines <- which(A1) + 1
    AllSeq <- currL[Slines]
    Qlines <- which(A1) + 3
    AllQual <- currL[Qlines]
    
    AllBc = future_sapply(AllSeq, function(x) substr(x,1,5))
    isBC = which(future_sapply(AllBc, function(x) sum(x == barcodes))==1)
    
    AllCord = AllCord[isBC]
    AllSeq = AllSeq[isBC]
    AllQual = AllQual[isBC]
    
    remove(currL)
    
    time2 <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    start_time <- Sys.time()
    
    ########
    RemSeq <- AllSeq
    RemQual <- AllQual
    RemCord <- AllCord
    
    start_time <- Sys.time()
    currentname <- makename(i, outname)
    savename = paste(outdir,paste0(currentname, ".txt") , sep="/")
    fh <- file(savename, open="w")
    
    for (i in 1:length(searchFall))  {
      
      
      for (ii in 1:length(searchRall))  {
        # Search for adapters with perfect alignment/mutation/N
        
        
        
        res <- find_row_col(RemSeq, searchFall[i], searchRall[ii], Lins)
        
        rowF <- res[[1]]
        colF <- res[[2]]
        RemI <- res[[3]]
        linsert <- res[[4]]
        
        
        LF = nchar(searchFall[i])
        LR = nchar(searchRall[ii])
        
        if( i == 1 && ii == 1){
          tag = "<PERF>"
        } else{
          tag = "<1M/N>"
        }
        savePar(tag, fh, rowF, colF, RemSeq, RemQual, RemCord, LF, LR, BARL, linsert)
        
        # The information that remained after removing the found strings
        RemSeq <- RemSeq[RemI]
        RemQual <- RemQual[RemI]
        RemCord <- RemCord[RemI]
        
        if(length(RemSeq) == 0){
          break
        }
        
        cat("\n", file = fh)
        
      }
      if(length(RemSeq) == 0){
        break
      }
    }
    timeF <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    close(fh)
    if(length(readLines(savename)) == 0) {
      file.remove(savename)
    }
    remove(RemSeq,RemQual,RemCord)
    
    
    # cat(currentname, "F:", timeF, "s\n")
    ###########
    
    
    
    # for (jj in seq_along(AllSeq)) {
    #   writeLines(paste0(AllSeq[jj]," ", AllQual[jj]," ", substr(AllCoord[jj], 10, nchar(AllCoord[jj])-6)), fh)
    # }
    # 
    # close(fh)
    # if(length(readLines(savename)) == 0) {
    #   file.remove(savename)
    # }
    
    time3 <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    print(paste("wrote", currentname, "Load=", time1, "Convert=", time2, "Save=", time3))
    
    
    
    
    
  }
  
  closeAllConnections()
  
}
