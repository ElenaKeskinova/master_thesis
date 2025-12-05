rawseq <- function(inname = "",chunk = 250000, outdir = "", outname = "rawseq", indir = "", barcodes = c() ) {
  
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
    outdir <- file.path(indir, paste0(substring(inname, 1, nchar(inname) - 4), "_RAW")) # default save directory
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
      
    AllCoord <- currL[A1]
    Slines <- which(A1) + 1
    AllSeq <- currL[Slines]
    Qlines <- which(A1) + 3
    AllQual <- currL[Qlines]
    
    AllBc = future_sapply(AllSeq, function(x) substr(x,1,5))
    isBC = which(future_sapply(AllBc, function(x) sum(x == barcodes))==1)
    
    AllCoord = AllCoord[isBC]
    AllSeq = AllSeq[isBC]
    AllQual = AllQual[isBC]
    
    remove(currL)
      
    time2 <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    start_time <- Sys.time()
      
    currentname <- makename(i, outname)
    savename = paste(outdir,paste0(currentname, ".txt") , sep="/")
    fh <- file(savename, open="w")
      
    for (jj in seq_along(AllSeq)) {
      writeLines(paste0(AllSeq[jj]," ", AllQual[jj]," ", substr(AllCoord[jj], 10, nchar(AllCoord[jj])-6)), fh)
    }
    
    close(fh)
    if(length(readLines(savename)) == 0) {
      file.remove(savename)
    }
      
    time3 <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    print(paste("wrote", currentname, "Load=", time1, "Convert=", time2, "Save=", time3))
    
    
    
    
    
  }
  
  closeAllConnections()
  
}
