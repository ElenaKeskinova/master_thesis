abs = c("Herceptin","21c","17b","b12")

library(future.apply)
library(igraph)
library(RandPro)
require(rgl)
source("add_weights.R")
# rnd library frequencies:
load(file = "4st_generation_library/lib_freqs.RData")
bgmot = freqM7

# add weights ----
for(ab in abs){
  lg = add_weights(ab,logw = T)
  
  path = paste0("mixed-7graphs/",ab,"/")
  save(lg, file = paste0(path,ab,"big7-logw.RData"))
}

# print all peptides to text files
for(ab in abs){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7or.RData"))
  
  writeLines(V(G)$name, file(paste0("allpeps_", ab,".txt")))
}


# spectral clustering ----
source("specclust.R")

for(ab in abs[1:3]){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-logw.RData"))
  specclust(ab,lg,bc = "logw")
}


## plot 3d graphics of clusters ----
ab = abs[2]
path = paste0("mixed-7graphs/",ab,"/")

load(paste0(path,ab,"big7-logw.RData"))

load(paste0(path,ab,"big7or.RData"))

cpl=colorRampPalette(c("#0050AA9F","#10AA109F","#50AF3055","#FFFF009F","#FFA0009F","#B50000"), alpha=T)
bc = "logw"
load(file = paste0(path,ab,"_",bc,"_","dbscan.RData"))
load(file = paste0(path,ab,"_",bc,"_","speccoord.RData"))
cldbsc = result_db[[1]]
ncl = result_db[[2]]
rndmx=form_matrix(ncol(coord),3, JLT=F)
mxsm=coord%*%rndmx
mxsm=mxsm/(sqrt(rowSums(mxsm^2)))
colnames(mxsm) = c("d1", "d2","d3")
plot3d(mxsm[cldbsc>0,],col=cpl(ncl-1)[cldbsc[cldbsc>0]])  
plot3d(mxsm,col=c("black",cpl(ncl-1))[cldbsc+1])  # non-clustered are black, clustered in colors
table(cldbsc)
 

# make pepsets ---- 
source("dbclust.R")

allpepsets_w = future_lapply(abs, \(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  bc = "logw"
  load(paste0(path,ab,"big7or.RData"))
  load(paste0(path,ab,"big7-",bc,".RData"))
  
  load(file = paste0(path,ab,"_",bc,"_","dbscan.RData"))
  load(file = paste0(path,ab,"_all_lcs.RData"))
  
  pepsets = gen_pepsets(lg,G,edg_lcs,result_db)
  save(pepsets, file = paste0(path,ab,"_peps_w.RData"))
  pepsets
  
})
names(allpepsets_w) = abs

allpepsets_w = lapply(abs, \(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(file = paste0(path,ab,"_peps_w.RData"))
  pepsets
})

allfreqs = future_lapply(1:4,\(i){
  ab = abs[i]
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7or.RData"))
  freqs= V(G)$Freq
  names(freqs) = V(G)$name
  
  sapply(allpepsets_w[[i]],\(set) freqs[set] )})

logfreqs = lapply(allfreqs,\(fr) sapply(fr,\(f) round(log2(f))))
logpepsets = lapply(1:4,\(i){
  lapply(1:length(logfreqs[[i]]),\(j){
    
    rep(allpepsets_w[[i]][[j]],logfreqs[[i]][[j]])
  })
})



clnames = sapply(1:4,\(i) sapply(1:length(logpepsets[[i]]),\(j) paste(abs[[i]],j,sep = "_")))
codes = unlist(sapply(1:4,\(i) rep(abs[i],length(logpepsets[[i]]))))
numcodes = unlist(sapply(1:4,\(i) rep(i,length(logpepsets[[i]]))))


for(i in 1:4){
  names(logpepsets[[i]]) = clnames[[i]]
}
for(i in 1:4){
  names(allpepsets_w[[i]]) = clnames[[i]]
}


# compare to clusters without weights ----
## tables with cluster distribution ----
for(ab in abs){
  bc = "logw"
  path = paste0("mixed-7graphs/",ab,"/")
  load(file = paste0(path,ab,"_",bc,"_","dbscan.RData"))
  cldbsc_w = result_db[[1]]
  load(file = paste0(path,ab,"dbscan.RData"))
  cldbsc = result_db[[1]]
  print(ab)
  print(table(cldbsc,cldbsc_w))
}



# lcs sets ----
all_lcssets = future_lapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  bc = "logw"
  
  load(paste0(path,ab,"big7-",bc,".RData"))
  
  load(file = paste0(path,ab,"_",bc,"_","dbscan.RData"))
  load(file = paste0(path,ab,"_all_lcs.RData"))
  
  require(vctrs)
  require(igraph)
  require(purrr)
  cldbsc = result_db[[1]]
  ncl = result_db[[2]]
  
  cls = sort(unique(cldbsc))
  i = which.max(components(lg)$csize)
  lgg = subgraph(lg,V(lg)[components(lg)$membership == i])
  
  lcssets = lapply(cls[which(cls!=0)], function(i) V(lgg)[which(cldbsc==i)]$name)
  lcsgsets = lapply(lcssets, \(set) {
    sg = subgraph(lg, set)
    c = which(igraph::components(sg)$csize >=10)
    lcscon = lapply(c, function(ci) V(sg)[which(igraph::components(sg)$membership == ci)]$name)
    lcscon
  })
  lcsgsets = list_flatten(lcsgsets)
  names(lcsgsets) = names(allpepsets_w[[ab]][1:length(lcsgsets)])
  lcsgsets
})


# merge clusters with high similarity, calculated froma scaled distance matrix ----

## b12: cluster 12,13,15----

peps121315 = c(allpepsets_w[[4]]$b12_12,allpepsets_w[[4]]$b12_13,allpepsets_w[[4]]$b12_15)
printlogo(peps121315)
allpepsets_w[[4]]$b12_1235m = peps121315

## 21c cluster 9 and 10

peps910 = c(allpepsets_w[[2]]$`21c_9`,allpepsets_w[[2]]$`21c_10`)
printlogo(peps910)
allpepsets_w[[2]]$'21c_910' = peps910


## list without original clusters

allpepsets_wm = allpepsets_w
allpepsets_wm[[2]][9:10] = NULL
allpepsets_wm[[4]][c(12,13,15)] = NULL


# logpepsets ----
allfreqs = future_lapply(1:4,\(i){
  ab = abs[i]
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7or.RData"))
  freqs= V(G)$Freq
  names(freqs) = V(G)$name
  
  sapply(allpepsets_wm[[i]],\(set) freqs[set] )})

logfreqs = lapply(allfreqs,\(fr) sapply(fr,\(f) round(log2(f))))
logpepsets = lapply(1:4,\(i){
  lapply(1:length(logfreqs[[i]]),\(j){
    
    rep(allpepsets_wm[[i]][[j]],logfreqs[[i]][[j]])
  })
})



clnames = sapply(1:4,\(i) sapply(1:length(logpepsets[[i]]),\(j) paste(abs[[i]],j,sep = "_")))
codes = unlist(sapply(1:4,\(i) rep(abs[i],length(logpepsets[[i]]))))
numcodes = unlist(sapply(1:4,\(i) rep(i,length(logpepsets[[i]]))))

for(i in 1:4){
  names(logpepsets[[i]]) = clnames[[i]]
}
for(i in 1:4){
  names(allpepsets_wm[[i]]) = clnames[[i]]
}

save(allpepsets_wm,logpepsets, file = "cluster_pepsets.RData")

# make motifs of new sets ----

require(universalmotif)
plan(multisession)
allm_w2 = future_lapply(logpepsets,\(pepsets){
  sapply(pepsets,\(set)create_motif(set))})
allm_wp2 =future_lapply(logpepsets,\(pepsets){
  sapply(pepsets,\(set)create_motif(set, pseudocount = 1))}) # with pseudocount

# pssm:
pssm_7_clean2 = future_lapply(logpepsets,\(pepsets){
  lapply(pepsets,\(set){
    mot = create_motif(set,type = "PWM",bkg = rowSums(freqM7)/7, pseudocount = 1)@motif
    mot = rbind(mot,rep(0,ncol(mot)))
    rownames(mot)[nrow(mot)] = "-"
    mot
  })
}) # with pseudocount

# calculate kld----
kld = function(p,q){
  sum(p*log(p/q))
}
all_kld = sapply(unlist(allm_wp2,recursive = F),\(mot) kld(mot@motif,bgmot))
lengths = sapply(unlist(logpepsets,recursive = F),length) 
