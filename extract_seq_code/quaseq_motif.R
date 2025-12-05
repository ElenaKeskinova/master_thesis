quaseq <- function(inname = "", rawname = "parseq", rawdir = "",
                      outdir = "", outname = "quaseq", indir = "/Users/ratmirderda/Documents/My PAPERS/Illumina Sequencing/",
                      RESCUE = 1, CUTOFF = "A", barcodes = "") {
  
  require(Biostrings)
  require(foreach)
  require(future)
  require(doFuture, quietly = TRUE)
  registerDoFuture()
  plan(multisession)
  
  if (inname == "") {
    inname <- file.choose() # If inname stays blank, call file.choose
  }
  
  if (rawdir == "") {
    rawdir <- file.path(indir, paste0(substr(inname, 1, nchar(inname) - 4), "_PAR")) # Default directory that contains raw files
  }
  
  if (outdir == "") {
    outdir <- file.path(indir, paste0(substr(inname, 1, nchar(inname) - 4), "_QUA")) # Default save directory
  }
  
  # Create outdir if it doesn't exist
  if (!file.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # Find all input files and cycle through all
  files <- list.files(path = rawdir, pattern = "parseq[[:xdigit:]]{4}.txt", full.names = TRUE)
  
  
  
   
  foreach (currentname = files) %dofuture% {
    start_time <- Sys.time()
    
    savename <- future_sapply(1:length(barcodes), function(x) file.path(outdir, paste0(outname,substr(basename(currentname),7,10), barcodes[x],".txt")))
    fh <- future_lapply(1:length(savename), function(x) file(savename[x], "w"))
    
    AllVar <- read.table(currentname, header = FALSE, colClasses = "character",comment.char = "", quote = "")
    R36Qual <- AllVar[ ,8]
    R36Bcode <- AllVar[ ,2]
    R36Seq <- AllVar[ ,4]
    R36Cord <- AllVar[ ,10]
    rm(AllVar)
    
    
    j <- 0
    j1 <- 0 #колко секвенции са преведени от всеки файл
    
    
    for (i in 1:length(R36Qual)) {
      # create temp variables to minimize access to the array
      n <- min(nchar(R36Seq[[i]]), nchar(R36Qual[[i]]))
      
      if (n >= 11) {
        BC <- R36Bcode[[i]]
        TS <- substr(R36Seq[[i]], 1, n)
        TQ <- substr(R36Qual[[i]], 1, n)
        CO <- R36Cord[[i]]
      } else {
        next  # if one of the strings is <11 bp, skip to the next
      }
      I <-  which((unlist(strsplit(TQ,"")))< CUTOFF) # позиции на некачествените
      DOT <- regexpr("\\.", TS, fixed = TRUE) # позиция на съвпадение или -1 ако няма
      
      if (length(I) == 0 && DOT == -1) {
        j <- j + 1
        SP <- paste(rep(" ", 36 - n + 3), collapse = "")  # space of variable Length
        SP2 <- strrep(" ", times = (14 - ceiling(n/3)))
        SP3 <- strrep(" ", times = (14 - floor(n/3)))
          if (n %% 3 == 0) {
            A12 <- as.character(Biostrings::translate(DNAString(TS)))
          } else if (n %% 3 == 2) {
            A12 <- Biostrings::translate(DNAString(paste0(TS, "K", collapse = "")),if.fuzzy.codon = "solve")
            A12 <- as.character(A12)
          
            # check if addition of K helped identification
            if (any(regexpr("X",A12) != -1)) {
              j = j-1
              next  # if codon is ambiguous, discard the seqeunce
            }
          } else {
            j <- j - 1  # ignore seqeunces shorter than 32 bp
            next
          }
        
        A12 <- unlist(strsplit(A12,""))
        # fix the TAG from stop to Glutamine (Q)
        # this stop codon is supressed in ER2738 (glnV or supE) strain
        STOP <- which(A12 == '*')
        if (length(STOP) > 0) {
          for (ii in 1:length(STOP)) {
            if (substr(TS,(STOP[ii] - 1) * 3 + 1,(STOP[ii] - 1) * 3 + 3) == 'TAG') {
              A12[STOP[ii]] <- 'Q'
              # ако има TAG кодон, се чете като глутамин
            }
          }
        }
       bcnum = match(BC, barcodes)
       if(is.na(bcnum) == TRUE){
         j = j-1
         next  
       } else {
         cat(sprintf("%s\n", paste(TS, SP, TQ, SP, paste0(A12, collapse = ""), SP2, isNNK(TS), SP3, CO)), file = fh[[bcnum]], append = TRUE)
      }
        }
        
        else if (sum(I %% 3) == 0 && DOT == -1) {
          # all bad reads are at the K positions or after 33rd base pair
          j1 <- j1 + 1
          TS = unlist(strsplit(TS,""))
          TS[I[I %% 3 == 0]] <- 'K'
          TS = paste(TS, collapse = "")
          SP <- strrep(' ', times = (36 - n + 3)) # space of variable Length
          SP2 <- strrep(" ", times = (14 - ceiling(n/3)))
          SP3 <- strrep(" ", times = (14 - floor(n/3)))
          if (n %% 3 == 0) {
            A12 <- as.character(translate(DNAString(TS), if.fuzzy.codon = "solve"))
          } else if (n %% 3 == 2) {
            A12 <- as.character(translate(DNAString(paste0(TS, 'K', collapse = "")), if.fuzzy.codon = "solve"))
          } else {
            j1 <- j1 - 1
            next
          }
          if (any(regexpr("X",A12) != -1)) {
              j1 = j1 - 1
              next # if any codon is ambiguous, discard the sequence
            }
        # fix TAG from stop to Glutamine (Q)
        # this stop codon is supressed in ER2738 (glnV or supE) strain
        A12 <- unlist(strsplit(A12,""))
        STOP <- which(A12 == '*')
        
        if (length(STOP) != 0) {
          for (ii in 1:length(STOP)) {
            if (substr(TS,(STOP[ii] - 1) * 3 + 1,(STOP[ii] - 1) * 3 + 3) == 'TAG') {
              A12[STOP[ii]] <- 'Q'
            } else {
              j1 = j1 - 1
              next
            }
          }
        }
        bcnum = match(BC, barcodes)
        if(is.na(bcnum) == TRUE){
          j1 = j1-1
          next  
        } else {
        cat(sprintf("%s\n", paste(TS, SP, TQ, SP, paste0(A12, collapse = ""), SP2, isNNK(TS), SP3, CO)), file = fh[[bcnum]], append = TRUE)
        }
        }
       
      }
    
        timeR <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        cat(sprintf("%s %d/%d seq (%.3f%%). Rescued %d seq (%.3f%%). Time: %s sec\n", currentname, j, i, (j/i)*100, j1, (j1/i)*100, timeR))
      
    } 
  closeAllConnections()
  }
        