find_row_col <- function(AllSeq, searchF, searchR, Lins) {
  require(stringr)
  
 if (length(AllSeq) > 0)  
 {
   posF <- str_locate_all(AllSeq, unlist(searchF))
  for (ii in seq_along(posF)) { #редовете без съвпадение получават стойност 0, тези със получават започващата позиция
    if (length(posF[[ii]]) == 0) {
      posF[[ii]] <- 0
    } else if (length(posF[[ii]]) > 1) {
      posF[[ii]] <- posF[[ii]][1]
    }
  }
  LF <- unlist(nchar(searchF))
  Seq = sapply(1:length(AllSeq), function(x) substr(AllSeq[x],(5+LF+12),nchar(AllSeq[x]))) # поне 11 бази между двата адаптера
  linsert = c()
  posR <- str_locate_all(Seq, unlist(searchR))
  for (ii in seq_along(posR)) { #редовете без съвпадение получават стойност 0, тези със получават започващата позиция
    if (length(posR[[ii]]) == 0) {
      posR[[ii]] <- 0
      linsert[ii] <- 0
      posF[[ii]] <- 0
    } else if (length(posR[[ii]]) > 1) {
      posR[[ii]] <- posR[[ii]][1]
      linsert[ii] <- posR[[ii]] + 10
      if(!is.null(Lins)){
        if(linsert[ii]!=Lins && linsert[ii]!=(Lins-1)){
          posF[[ii]] <- 100
        }
        
      } 
    }
  }
  
  
  
  posF <- unlist(posF)
  col1 <- posF[posF == 6]
  row1 <- which(posF == 6)
  linsert <- linsert[row1]
  Remain <- which(posF == 0) #номерата на редове без съвпадение
  return(list(row1, col1, Remain, linsert))
 }
  else {
    return(list(0,0,c(),0))
  }
}
