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

# make df with all lcs for each sequence ----

plan(multisession(workers = 4))
foreach(ab = abs) %dopar% {
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  load(file = paste0(path,ab,"_all_lcs.RData"))
  
  p = V(G)$name
  plcs = lapply(p,\(pi){
    c(edg_lcs$enames[which(edg_lcs$X1 == pi)],edg_lcs$enames[which(edg_lcs$X2 == pi)])
  })
  pep_lcs = data.frame(p)
  pep_lcs$lcs = plcs
  
  plcs5 = lapply(plcs,\(set){
    c(unlist(sapply(set,\(pep){
      if(nchar(pep)==6){
        sapply(1:6,\(i) paste0(substr(pep,1,(i-1)),substr(pep,(i+1),6)))
      }
      else{
        pep
      }
    })))
  })
  
  pep_lcs$lcs5 = plcs5
  save(pep_lcs,file = paste0(path,ab,"_pep_lcs.RData"))
  
}

plan(sequential)

# add weights ----
source("add_weights.R")
for(ab in abs){
  lg = add_weights(ab,logw = T)
  
  path = paste0("mixed-7graphs/",ab,"/")
  save(lg, file = paste0(path,ab,"big7-logw.RData"))
}


# extract parameters and plot degree vs copy number ----
par(mfrow = c(2,2), mar = c(4,4,2,2))
for(ab in abs){
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  load(file = paste0(path,ab,"big7-logw.RData"))
  print(paste(ab, vcount(G), ecount(G), vcount(lg), ecount(lg), edge_density(G),edge_density(lg), V(G)$name[which.max(V(G)$Freq)]))
  plot(log2(degree(G)),log10(V(G)$Freq), xlab = "", ylab = "", main = ab, ylim = c(0,7))
  lines(c(0,max(log2(degree(G)))), rep(log10(2000),2),col = "red")
  #lines( rep(log2(median(degree(G))),2),range(log10(V(G)$Freq)),col = "red")
  if(ab %in% abs[3:4]){mtext("log2(vertex degree)", side = 1, line = 3) }
  if(ab %in% abs[c(1,3)]){mtext("log10(peptide copy number)", side = 2, line = 3) }
  
  i = which(V(G)$Freq>=2000)
  #print(length(which(degree(G)[i]>=8))/length(i))
  print(max(degree(G)[i]))
  mi = median(degree(G,i))
  mii = median(degree(G,which(V(G)$Freq<2000)))
  r = cor(log(degree(G))[i],log(V(G)$Freq)[i], method = "spearman")
  points(log2(degree(G))[i],log10(V(G)$Freq)[i],col = "red")
  text(0.1,7, paste("r =",round(r,2)),col = "red",adj = 0)
  text(0.1,6, paste("med(degree) =",round(mi,2)),col = "red",adj = 0)
  text(0.1,5, paste("med(degree) =",round(mii,2)),col = "black",adj = 0)
  x1 = degree(G)>=8
  x2 = V(G)$Freq>=2000
  #print(chisq.test(table(x1,x2)))
}

