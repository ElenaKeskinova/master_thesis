savePar <- function(TAG, fhF, rowF, colF, AllSeq, AllQual, AllCord, LF, LR, BARL, Linsert) {
  
  # for truncated nucleotides LF can be cell array containing original length
  # and truncated length: LF = {LFtruncated; LForiginal}
  if (is.list(LF)) {
    LFori <- LF[2]
    LF <- LF[1]
  } else {
    LFori <- LF
  }
  
  
    
    for (j in seq_along(rowF)) {
    # creating temporary static variables to minimize acess to array and
    # acelerate the program
    temp1 <- AllSeq[rowF[[j]]]
    temp2 <- AllQual[rowF[[j]]]
    temp3 <- AllCord[rowF[[j]]]
    addr <- colF[j]
    Lins <- Linsert[j]
    
    # this portion checks whether there are enough symbols before the search
    # strings for BAR and NKK. If not, the missing symbols are filled with a
    # line of dots '......'  (dot has an ASCII number of 46)
    if ((addr-BARL) < 1) {
      temp1 <- paste(strrep(".", (BARL+1-addr)), temp1, sep="")
      temp2 <- paste(strrep(".", (BARL+1-addr)), temp2, sep="")
      addr <- BARL+1
    }
    # if ((addr-NKKL-BARL+LF) < 1) {
    #   temp1 <- paste(rep(".", NKKL+BARL+1-addr), temp1, sep="")
    #   temp2 <- paste(rep(".", NKKL+BARL+1-addr), temp2, sep="")
    #   addr <- NKKL+BARL+1
    # }
    
    L1 <- nchar(temp1)
    L2 <- nchar(temp2) 
    BAR <- substr(temp1, addr-BARL, addr-1)
    St36 <- addr+LF #начало на секвенция
    En36 <- addr+LF+Lins-1 #край на секвенция
    FOR <- substr(temp1, addr, St36-1) # ляв адаптер
    FOR <- paste(strrep(".", LFori-LF), FOR, sep="")# допълване с точки отпред
    R36 <- substr(temp1, St36, pmin(En36, L1)) # секвенцията, или до края на стринга
    R36 <- paste(R36, strrep(" ", (37-nchar(R36))), sep="") # допълване с интервали до 36 букви
    REV <- substr(temp1, En36+1, pmin(En36+LR, L1)) # десен адаптер
    REV <- paste(REV, strrep(".", LR-nchar(REV)), sep="") # допълнен с точки до края
    
    QBAR <- substr(temp2, addr-BARL, addr-1)
    QFOR <- substr(temp2, addr, St36-1)
    QFOR <- paste(strrep(".", LFori-LF), QFOR, sep="")
    QR36 <- substr(temp2, St36, pmin(En36, L2))
    QR36 <- paste(QR36, strrep(" ", (37-nchar(QR36))), sep="")
    QREV <- substr(temp2, En36+1, pmin(En36+LR, L2))
    QREV <- paste(QREV, strrep(".", LR-nchar(QREV)), sep="")
    
    writeLines(paste0(TAG, ' ', BAR, ' ', FOR, ' ', R36, ' ', REV, strrep(' ', 10), 
                  QBAR, ' ', QFOR, ' ', QR36, ' ', QREV, '   ', temp3), 
               fhF)
    }
    
  
}