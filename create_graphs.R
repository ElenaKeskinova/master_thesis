library(future)
library(foreach)
library(igraph)
library(stringr)

source("newgraph.R")
source("graph_to_line.R")

abs = c("Herceptin","21c","17b","b12")
f = file("TUPs.txt") # created with TUPSCAN from all peptides
tups = read.table(file = f,colClasses = "character",col.names = c("Tups","motif","inf"), sep ="_")
tups$Tups = sapply(tups$Tups, \(t) str_remove(t,"  "))

ncores = parallel::detectCores()
plan(multisession(workers = ncores-2))

# make graphs ----
foreach(ab = abs) %dopar% {
  path = paste0("mixed-7graphs/",ab,"/")
  f = file(paste0(path,ab,"_allp810_noC7-fr.txt"))
  Peps = read.table(file = f,colClasses = c("character","integer"),col.names = c("Pep","Freq"))
  Peps = Peps[which(Peps$Freq>1),] # only with frequency>1
  ii = which(!(Peps$Pep %in% tups$Tups))
  Peps = Peps[ii,]
  peps = Peps$Pep
  
  G = newgraph(peps,2)
  print(paste(ab, vcount(G), ecount(G)))
  G = set_vertex_attr(G,name = "Freq",value = Peps$Freq)
  #print(components(G)$csize)
  save(G,file = paste0(path,ab,"big7or.RData"))
  
}
closeAllConnections()


# make line graphs ----
plan(multisession(workers = ncores-2))
foreach(ab = abs) %dopar% {
  source("compute_lcs.R")
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  
  lg = graph_to_line(G)
  
  print(paste(ab, vcount(lg), ecount(lg)))
  save(lg,file = paste0(path,ab,"big7-l.RData"))
  
}

plan(sequential)

# make data frame with all lcs for every edge ----

plan(multisession(workers = ncores-2))
foreach(ab = abs) %dopar% {
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  e = ends(G,E(G))
  enames = future_sapply(E(G), function(x) compute_lcs(ends(G,x)[1],ends(G,x)[2]))
  edg_lcs = data.frame(e,enames)
  
  save(edg_lcs,file = paste0(path,ab,"_all_lcs.RData"))
  
}

plan(sequential)


