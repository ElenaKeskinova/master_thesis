library(Biostrings)
library(igraph)
library(future.apply)
require(universalmotif)
source("newgraph.R")
source("graph_to_line.R")
source("dbclust.R")
source("rnd_graph.R")
require(parallel)
load(file = "lib_freqs.RData")
#freq = consensusMatrix(allpep)

abs = abs = c("Herceptin","21c","17b","b12")
gsizes = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  vcount(G)
})

repfreqs = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  round(log2(V(G)$Freq))
})

rep_freqs = table(unlist(repfreqs))/sum(table(unlist(repfreqs))) # frequencies of logarithmic copy numbers

#bkg_m = create_motif(peps_7nC)    
bgmot = freqM7
freq = rowSums(freqM7)/sum(freqM7)

n = 100 #all graphs 
ncores = detectCores()-2
sizes = sample(min(gsizes):10000,n,replace = T)

require(dbscan)
require(universalmotif)
require(igraph)
plan(multisession(workers = 20))
rnd_pepsets = future_lapply(1:n,\(i){ ### iterations for all graphs
  vc= sizes[i]
  lg = rnd_lgraph_repfreqs(vc,freqs=freq,repfreqs=rep_freqs)
   
  gc()
  #print(c(i,"edges" ,vcount(lg)))
  # spec clust and dbscan
  print(paste(i,"graph"))
  
  pepsets = dbclust_peps(lg)
  
  print(paste(i,"dbcl ready"))
  
  pepsets
  
},future.seed = T, future.chunk.size = 5) #
plan(sequential)
rnd_pepsets = unlist(rnd_pepsets,recursive = F)
save(rnd_pepsets,file = "rnd_pepsets_bkgfreqs_2.RData" )
save(rnd_pepsets,file = "rnd_pepsets_bkgfreqs_3.RData" ) # 26/02/2026 for mthesis

load(file = "rnd_pepsets_bkgfreqs_2.RData")
# make pssm without alignment
library(universalmotif)
library(future.apply)
plan(multisession(workers = 20))
rnd_pssm_bkg = future_lapply(rnd_pepsets,\(set){
    mot = create_motif(set,type = "PWM",bkg = rowSums(freqM7)/7, pseudocount = 1)@motif
    mot = rbind(mot,rep(0,ncol(mot)))
    rownames(mot)[nrow(mot)] = "-"
    mot
}) # with pseudocount
save(rnd_pssm_bkg, file =  "rnd_pssm_bkgfr_2.RData")
save(rnd_pssm_bkg, file =  "rnd_pssm_bkgfr_3.RData")


kld = function(p,q){
  sum(p*log(p/q))
}

# background:

load(file = "lib_freqs.RData")
bgmot = freqM7

# calculate kld
plan(multisession)
rnd_len_kld = t(future_sapply(rnd_pepsets,\(set){
  m = create_motif(set,pseudocount = 1,alphabet = 'AA')
  c(length(set),kld(m@motif,bgmot),length(unique(set)))
  
}))
colnames(rnd_len_kld) = c("l","kld", "un l")

save(rnd_len_kld,file="rnd_len_kld_3.RData")
plot(log10(rnd_len_kld))

# rnd graphs for comparison to epitopes

ag_aaprob = read.csv(file = "AgAbIFprobs.csv")
aa_prob = ag_aaprob$p
names(aa_prob) = names(AMINO_ACID_CODE[1:20])
aa_prob = aa_prob[order(names(aa_prob))]

n = 100 #all graphs 
ncores = detectCores()-2
sizes = sample(min(gsizes):max(gsizes),100,replace = T)

require(dbscan)
require(universalmotif)
require(igraph)
plan(multisession(workers = 20))
rnd_pepsets = future_lapply(1:n,\(i){ ### iterations for all graphs
  vc= sizes[i]
  rg = rnd_graph_reps(vc,freqs=aa_prob,repfreqs=rep_freqs,s=i+200)
  
  g = rg[[1]]
  lg = rg[[2]]
  rm(rg)
  gc()
  #print(c(i,"edges" ,vcount(lg)))
  # spec clust and dbscan
  print(paste(i,"graph"))
  pepsets = dbclust_peps(g,lg)
  
  print(paste(i,"dbcl ready"))
  
  v_fr = V(g)$Freq
  names(v_fr) = V(g)$name
  freqs = sapply(pepsets,\(set) v_fr[set] )
  newsets = lapply(1:length(freqs),\(ii){
    rep(pepsets[[ii]],freqs[[ii]])
  })
  newsets
  
},future.seed = T) #
plan(sequential)
rnd_pepsets = unlist(rnd_pepsets,recursive = F)
save(rnd_pepsets,file = "rnd_pepsets_abagfreqs.RData")


rnd_align = sapply(rnd_pepsets,\(pepset){
    require(msa)
    l=msaClustalW(AAStringSet(pepset), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
    l= apply(as.matrix(AAStringSet(l)),1, paste,collapse="")
    l
    
  })

rnd_ppm = future_lapply(rnd_align,\(l){

    freq_matrix(l,AA_STANDARD,ps_c = 1)
  })


save(rnd_align,rnd_ppm,file = "rnd_ppm_abagfreqs.RData")

