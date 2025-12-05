makename <- function(num, inname) {
  Snum <- formatC(num, width = 4, flag = "0")
  outname <- paste0(inname, Snum)
  return(outname)
}
