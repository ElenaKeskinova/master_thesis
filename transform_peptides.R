

# transform data from length 8 and 10 to length 7 ----
abs = c("Herceptin","b12","17b","21c")
for(ab in abs[1:4]){
  path = paste0("Exp12.f_UNI-10/", ab,"/")
  files <- list.files(path = path, pattern = "uniqueP_[[:upper:]]{5}.txt", full.names = TRUE)
  path8 = paste0("Exp12.f_UNI-8/", ab,"/")
  files8 <- list.files(path = path8, pattern = "uniqueP_[[:upper:]]{5}.txt", full.names = TRUE)
  
  #get all peps together for length 10
  AllPep = data.frame()
  for (f in files){
    Pep = read.table(file = f,
                     colClasses = c("character","integer"),col.names = c("Pep","Freq"))
    AllPep = rbind(AllPep,Pep)
    
  }
  AllPep_agg <- aggregate(AllPep[, 2], by = list(Pep = AllPep$Pep), FUN = sum)
  Peps_noC = sapply(AllPep_agg$Pep, function(p) {p1 = strsplit(p,"")[[1]]; p1[1]!="C" })
  Peps_noC = AllPep_agg[Peps_noC,]
  
  
  ## for 10 letters ----
  
  Peps7 = sapply(Peps_noC$Pep, \(p){
    ps = strsplit(p,"")[[1]]
    p1 = paste0(ps[1:7],collapse = "")
    p2 = paste0(ps[4:10],collapse = "")
    return(c(p1,p2))
  })
  
  freq10to7 = rep(Peps_noC$x,each = 2)
  
  ## for 8 letters ----
  AllPep = data.frame()
  for (f in files8){
    Pep = read.table(file = f,
                     colClasses = c("character","integer"),col.names = c("Pep","Freq"))
    AllPep = rbind(AllPep,Pep)
    
  }
  AllPep_agg <- aggregate(AllPep[, 2], by = list(Pep = AllPep$Pep), FUN = sum)
  
  
  
  
  Peps_noC8 = sapply(AllPep_agg$Pep, function(p) {p1 = strsplit(p,"")[[1]]; p1[1]!="C" })
  
  Peps_noC8 = AllPep_agg[Peps_noC8,]
  
  
  
  
  Peps87 = sapply(Peps_noC8$Pep, \(p){
    p1 = strsplit(p,"")[[1]]
    p1 = paste0(p1[1:7],collapse = "")
    p1
  })
  
  freq8to7 = Peps_noC8$x
  
  ## merge both dataframes ----
  pep = c(Peps7,Peps87)
  freq = c(freq10to7,freq8to7)
  all8107 = aggregate(freq,by = list(Pep = pep),FUN = sum)
  
  #write file with all peps and frequencies
  pathin = paste0("mixed-7graphs/", ab,"/")
  Fin = file(paste0(pathin,ab,"_allp810_noC7-fr.txt"),"w")
  for (ii in 1:length(all8107[[1]])){
    cat(sprintf("%s\n", paste0(all8107$Pep[ii], "    ", all8107$x[ii])), file = Fin, append = TRUE)
    
  }
  close(Fin)
  
  
  #only peps without frequencies
  Fin = file(paste0(pathin,ab,"_allp810_noC7-fr.txt"),"r")
  
  p = read.table(file = Fin,
                 colClasses = c("character","integer"),col.names = c("Pep","Freq"))$Pep
  
  Fout = file(paste0(pathin,ab,"_allp810_clean.txt"),"w")
  cat(p,file = Fout, sep = "\n")
  closeAllConnections()
  
  # only peps with frequencies >1 , for graphs
  
}  
 

