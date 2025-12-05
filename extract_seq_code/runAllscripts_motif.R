# default names will be overwritten by file.choose() function
name <- c('s_1_sequence.txt', 's_2_sequence.txt')
indir <- "D:/Elena/ban/fastq_extract_seq"

# select files using file.choose() function
files <- file.choose()
name = basename(files)
indir = dirname(files)

setwd("D:/Elena/ban/Motifier_DataSet/extract_seq_code")
source("rawseq_motif.R")
source("parseq_motif.R")
source("quaseq_motif.R")
source("savePar_motif.R")
source("uniseq_motif.R")
source("find_row_col_motif.R")
source("isNNK.R")
source("makename.R")
source("rawparseq_elena.R")

os = 'win'
checkNNK = TRUE

# check if files are ok; open first 10 lines of each file 
# if only one file was selected, convert to character vector
if (!is.character(name)) {
  temp <- name
  name <- character(1)
  name[1] <- temp
}

for (i in seq_along(files)) {
  filepath <- files[i]
  cat('Checking file', filepath, '\n')
  cat('**********************************************\n')
  FID <- file(filepath, 'r')
  currL <- character(10)
  for (ii in 1:10) {
    currL[ii] <- readLines(FID, n = 1)
    cat(currL[ii], '\n')
    # fpout <- file.path(indir, paste0('test', ii, '.txt'))
    # FIDo <- file(fpout, 'w')
    # writeLines(currL[ii], FIDo)
  }
  cat('**********************************************\n')
  cat('\n')
}

# forward and reverse adapter sequences flanking the variable region
sF <- 'CAACGTGGC'
sR <- 'GCCT'
linsert <- 24
qC <- "D" # cutoff for quality of the nucleotides (1:40)
# barcodes = c("AACTG", "AAGCT", "AATGC","ACCGT", "ACGTT", "ACTGG",
#              "AGCCT", "AGGTC", "AGTCC", "ATCGG", "ATGCC", "ATTCG",
#              "CAAGT", "CAGTT", "CATGG", "CCATG", "CCGAT", "CCTGA",
#              "CGAAT", "CTAGC", "GAATC", "GATCA", "GCGTA", "GCTAG",
#              "GGACT", "GTCAA", "TAACG", "TGAAC", "TGCTA", "TTGCA") #all

barcodes = c("TGCTA","AAGCT")  #17b:"TGCTA","AAGCT" 

for (i in seq_along(files)) {
  
  # Start <- Sys.time()
  # rawparseq(inname= name[i], indir=indir, chunk=1500000, barcodes = barcodes, 'searchF'=sF, 'searchR' =sR,'Lins' = linsert) #1500000 leads to error
  # rawsTIME1 <- Sys.time() - Start
  # cat('************ rawparseq TIME =', as.numeric(rawsTIME1, units = 'secs'), '********************\n')

  Start <- Sys.time()
  rawseq(inname= name[i], indir=indir, chunk=1500000, barcodes = barcodes) #1500000 leads to error
  rawsTIME1 <- Sys.time() - Start
  cat('************ rawseq TIME =', as.numeric(rawsTIME1, units = 'secs'), '********************\n')

  Start <- Sys.time()
  parseq('inname'=name[i], 'indir'=indir, 'searchF'=sF, 'searchR' =sR,'Lins' = linsert )
  parTIME1 <- Sys.time() - Start
  cat('************ parseq TIME =', as.numeric(parTIME1, units = 'secs'), '********************\n')

  Start <- Sys.time()
  quaseq('inname' = name[i], 'indir' = indir, 'CUTOFF' = qC, barcodes = barcodes)
  parTIME1 <- Sys.time() - Start
  cat('************ quaseq TIME =', as.numeric(parTIME1, units = 'secs'), '********************\n')
  
  Start <- Sys.time()
  uniseq('inname' = name[i], 'indir' = indir, 'os' = os, 'checkNNK' = checkNNK, barcodes = barcodes)
  uniTIME1 <- Sys.time() - Start
  cat('************ uniseq TIME =', as.numeric(uniTIME1, units = 'secs'), '********************\n')
}
 
 closeAllConnections()

# optional block which generates truncated sequences
# for (i in seq_along(files)) {
#   truncseq('inname' = name[i], 'indir' = indir, 'os' = os)
# }
