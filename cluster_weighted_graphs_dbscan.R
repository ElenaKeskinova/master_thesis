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

## results from spectral clustering ----
bc = "logw"
dbsc_results = lapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(file = paste0(path,ab,"_",bc,"_","dbscan.RData"))
  result_db
})

spect_results = lapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(file = paste0(path,ab,"_",bc,"_","speccoord.RData"))
  coord
})

### plots for thesis ----
par(mfrow = c(2,2), mar = c(3,3,2,2))
for(i in 1:4){
  plot(dbsc_results[[i]][[5]],dbsc_results[[i]][[4]], xlab =  "",ylab =  "",main =abs[i])
  j = which.max(dbsc_results[[i]][[4]])
  points(dbsc_results[[i]][[5]][j],dbsc_results[[i]][[4]][j], col =  "red",pch = 16)
  if(i>2){mtext("radius", side = 1, line = 3)}
  if(i%%2==1){
    mtext("clusters",side = 2,line = 3)
  }
}

for(a in 1:4){
  ab = abs[a]
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-logw.RData"))
  i = which.max(igraph::components(lg)$csize)
  g = subgraph(lg,V(lg)[igraph::components(lg)$membership == i])
  bc = "logw"
  load(file = paste0(path,ab,"_",bc,"_","dbscan.RData"))
  cldbsc = result_db[[1]]
  ncl = result_db[[2]]
  subgraphs = sapply(1:(ncl-1),\(i){
    subgraph(g, which(cldbsc == i))
  })
  compno = sapply(subgraphs, \(g) igraph::components(g)$no)
  hist(compno,main = ab,xlab =  "", breaks = seq(0,max(compno)+1))
  if(a>2){mtext("connected components per cluster", side = 1, line = 3)}
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

db_lcssets = lapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  
  load(paste0(path,ab,"big7-logw.RData"))
  bc = "logw"
  load(file = paste0(path,ab,"_",bc,"_","dbscan.RData"))
  load(file = paste0(path,ab,"_",bc,"_","speccoord.RData"))
  
  opdim = ncol(coord)
  gen_lcssets(lg,result_db,opdim)
})
db_lcssets5 = lapply(db_lcssets,\(absets){
  lapply(absets,\(set){
    c(unlist(sapply(set,\(pep){
      if(nchar(pep)==6){
        sapply(1:6,\(i) paste0(substr(pep,1,(i-1)),substr(pep,(i+1),6)))
      }
      else{
        pep
      }
    })))
  })
})

names(db_lcssets5) = abs
plan(multisession(workers = 4))
db_pepsets = future_lapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(file = paste0(path,ab,"_pep_lcs.RData"))
  lcs_to_peps(db_lcssets5[[ab]],pep_lcs)
})

## plot logos of biggest clusters ----
par(mfrow = c(1,4))
for(i in 1:4){
  l = sapply(db_pepsets[[i]],length)
  logoal(db_pepsets[[i]][[which.max(l)]])
  print(max(l))
}



allfreqs = future_lapply(1:4,\(i){
  require(igraph)
  ab = abs[i]
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7or.RData"))
  freqs= V(G)$Freq
  names(freqs) = V(G)$name
  
  sapply(db_pepsets[[i]],\(set) freqs[set] )})

logfreqs = lapply(allfreqs,\(fr) sapply(fr,\(f) round(log2(f))))

logpepsets = lapply(1:4,\(i){
  lapply(1:length(logfreqs[[i]]),\(j){
    
    rep(db_pepsets[[i]][[j]],logfreqs[[i]][[j]])
  })
})




clnames = sapply(1:4,\(i) sapply(1:length(logpepsets[[i]]),\(j) paste(abs[[i]],j,sep = "_")))
codes = unlist(sapply(1:4,\(i) rep(abs[i],length(logpepsets[[i]]))))
numcodes = unlist(sapply(1:4,\(i) rep(i,length(logpepsets[[i]]))))


for(i in 1:4){
  names(logpepsets[[i]]) = clnames[[i]]
}
for(i in 1:4){
  names(db_pepsets[[i]]) = clnames[[i]]
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


save(db_pepsets,db_lcssets,logpepsets, file = "cluster_pepsets.RData")

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
