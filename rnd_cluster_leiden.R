library(Biostrings)
library(igraph)
library(future.apply)
require(universalmotif)
source("newgraph.R")
source("graph_to_line.R")
source("rnd_graph.R")
source("cl_cores_functions.R")
require(parallel)
load(file = "lib_freqs.RData")

freq = rowSums(freqM7)/sum(freqM7)
n = 10 #all graphs 
vc = 10000


repfreqs = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  round(log2(V(G)$Freq))
})

rep_freqs = table(unlist(repfreqs))/sum(table(unlist(repfreqs))) # frequenciesof logarithmic copy numbers

r = 0.03
plan(multisession(workers = 20))
rnd_pepsets = future_lapply(1:n,\(i){ ### iterations for all graphs
  
  #freqi = sample(freq, 20)
  #names(freqi) = names(freq)
  
  lg = rnd_lgraph_repfreqs(vc,freqs = freq, repfreqs = rep_freqs)
  E(lg)$weight = E(lg)$weight/max(E(lg)$weight)
  gc()
  #print(c(i,"edges" ,vcount(lg)))
  # spec clust and dbscan
  print(paste(i,"graph"))
  
  ncl = 5
  cls = sapply(1:ncl,\(i) cluster_leiden(lg,resolution = r, objective_function = "CPM",n_iterations = 10)$membership)
  m = 3
  cls2 = sapply(1:ncl,\(i) cut_small_clusters(cls[,i],m))
  ari = (future_sapply(1:(ncl),\(i){sapply((1):ncl,\(j) {adjustedRandIndex(cls2[,i],cls2[,j]) })})) # mean ARI, better overlap if high
  cl = cls[,which.max(colSums(ari))]
  m = 3
  lcssets = cl_cores(lg,cl,m,0)
  lcssets = lapply(lcssets,\(set){
    c(unlist(sapply(set,\(pep){
      if(nchar(pep)==6){
        sapply(1:6,\(i) paste0(substr(pep,1,(i-1)),substr(pep,(i+1),6)))
      }
      else{
        pep
      }
    })))
  })
  
  pepsets = lcs_to_peps(lcssets,lg$pep_lcs)
  pepsets
},future.seed = T) #

plan(sequential)
rnd_pepsets = unlist(rnd_pepsets,recursive = F)
rndlen = sapply(rnd_pepsets,length)



save(rnd_pepsets,file = "rnd_pepsets_libfreqs_leiden.RData" )
#save(rnd_lcssets,file = "rnd_lcssets_rndfreqs.RData" )
rnd_pepsets_3 = rnd_pepsets[rndlen>2]
rnd_pssm = future_lapply(rnd_pepsets_3,\(set){
  mot = create_motif(set,type = "PWM",alphabet = "AA",bkg = rowSums(freqM7)/7, pseudocount = 1)@motif
  mot = rbind(mot,rep(0,ncol(mot)))
  rownames(mot)[nrow(mot)] = "-"
  mot
})
save(rnd_pepsets_3,rnd_pssm,file = "rnd_pepsets_pssms_libfreqs_leiden.RData")

