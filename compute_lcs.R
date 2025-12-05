compute_lcs = function(a,b){ # return lcs string of 2 input strings
  require(qualV)
  lcs = qualV::LCS(strsplit(a, split = NULL)[[1]],strsplit(b,split = NULL)[[1]])$LCS
  lcs = paste0(lcs,collapse = "")
  return(lcs)
}