isNNK <- function(instring) {
  I <- (1:12) * 3
  n <- floor(nchar(instring) / 3)
  outstring = sapply(I[1:n], function(x) if(substr(instring,x,x) == "A" | substr(instring,x,x) == "C") {substr(instring,x,x)} else { "."})
  outstring = paste0(outstring, collapse = "")
  return(outstring)
}
