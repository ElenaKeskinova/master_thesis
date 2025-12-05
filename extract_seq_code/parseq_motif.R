parseq <- function(inname = '',  rawname = 'rawseq', rawdir = '' ,  outdir = '',  outname = 'parseq' ,  searchF = '',  searchR = '', indir = '/Users/ratmirderda/Documents/My PAPERS/Illumina Sequencing/', Lins = NULL) {
  # New version of Parseq script
  # it reads text file and saves the output to the text file.
  # so far works with one file, but its easy to for-loop it
  # to run it parseq('inname', 's_1_sequence.txt'); parseq('inname', 's_2_sequence.txt');
  # AP added Lins variable for the length of the seq
  # check whether indir name was defined
  
  if (inname == '') {
    message('using working R directory as current directory; Make sure your files are in the working R directory')
  }
  
  if (rawdir == '') {
    rawdir <- file.path(indir, paste0(substring(inname, 1, nchar(inname)-4), '_RAW')) # default directory that contains raw files
  }
  
  if (outdir == '') {
    outdir <- file.path(indir, paste0(substring(inname, 1, nchar(inname)-4), '_PAR')) # default save directory
  }
  ############################
  
  #outnameF <- paste0(outname, 'F')
  #outnameR <- paste0(outname, 'R')
  outnameREM <- paste0(outname, 'REM')
  
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  
  ############ End of the Input Section ##########
  ###############################################
  
  ###### Create nucleotide information in this section
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
  
  
  # Find all files input files and cycle through all
  files <- list.files(path = rawdir, pattern = 'rawseq[[:xdigit:]]{4}.txt', full.names = TRUE)
  
  require(foreach)
  require(future)
  require(doFuture, quietly = TRUE)
  
  plan(multisession)
  
  foreach (currentname = files) %dofuture% {
    N <- as.numeric(substr(currentname, nchar(currentname) - 7, nchar(currentname) - 4))
    savename <- file.path(outdir, paste0(makename(N, outname), ".txt"))
    
    
    
    AllVar <- read.table(currentname, header = FALSE, colClasses = "character", fill = TRUE, comment.char = "", quote = "")#*
    AllSeq <- AllVar[ ,1]
    AllQual <- AllVar[ ,2]
    AllCord <- AllVar[ ,3]
    rm(AllVar)
    
    searchFall = c(searchF, searchFM, searchFN)
    searchRall = c(searchR, searchRM, searchRN)
    
    RemSeq <- AllSeq
    RemQual <- AllQual
    RemCord <- AllCord
    fh <- file(savename, "w")
    
    start_time <- Sys.time()
    
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
    
    
    cat(currentname, "F:", timeF, "s\n")
  }
}
    