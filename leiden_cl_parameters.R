library(future.apply)
library(igraph)
library(mclust)
source("cl_cores_functions.R")
source("printlogo.R")
library(universalmotif)
abs = c("Herceptin","21c","17b","b12")
# test parameters for leiden clustering and choose best parameters

# cluster statistics function ----
leidcl_stats = function(lg, reso, n = 5, m = 3 ){ # graph, resolution, n = iterations per resolution, m = minimal cluster size
  t = future_sapply(reso,\(res){ # change resolution from 0.1 to 1
    
    
    cls = sapply(1:n,\(i) cluster_leiden(lg,resolution = res, objective_function = "CPM",n_iterations = 10)$membership)
    
    cls2 = sapply(1:n,\(i) cut_small_clusters(cls[,i],m))
    ari = unlist(future_sapply(1:(n-1),\(i){sapply((i+1):n,\(j) {adjustedRandIndex(cls2[,i],cls2[,j]) })})) # mean ARI, better overlap if high
    clcount = mean(sapply(1:n,\(i)length(unique(cls2[,i]))))
    nontrash = mean(sapply(1:n,\(i) sum(cls2[,i]!=0)/nrow(cls2)))
    c(mean(ari),sd(ari),clcount,nontrash)
  },simplify =  "array", future.chunk.size = ceiling(length(reso)/20))	

  t = rbind(reso,t)
  rownames(t) = c(  "reso","mean ARI", "sd ARI",  "clusters",  "part > minsize")
  return(t)
}

# ARI with resolutions from 1e-5 to 1e-1 ----
resos = 10**(seq(-5,-1,1))

res1_nw = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-l.RData"))
  leidcl_stats(lg,resos)
},simplify =  "array")

res1_w = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-logw.RData"))
  E(lg)$weight = E(lg)$weight/max(E(lg)$weight)
  leidcl_stats(lg,resos)
},simplify =  "array")

res1_rw = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-logw.RData"))
  E(lg)$weight = E(lg)$weight/max(E(lg)$weight)
  E(lg)$weight = sample(E(lg)$weight,ecount(lg))
  leidcl_stats(lg,resos)
},simplify =  "array")

par(mfrow = c(1,1))
for(res1 in list(res1_nw,res1_w,res1_rw)){
  plot(log10(res1[1,,1]),res1[2,,1],col = 1, ylim = c(0,1),pch = 16)
  for(i in 2:4){
    points(log10(res1[1,,i]),res1[2,,i],col = i,pch = 16)
  }
  # for(i in 1:4){
  #   points(log10(res1[1,,i]),res1[5,,i],col = i,pch = 3)
  # }
  legend( "bottomleft",legend = abs,col = seq(1,4),pch = 1)
}

cls = lapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-l.RData"))
  cl = cluster_leiden(lg,resolution = 10**(-5), objective_function = "CPM",n_iterations = 10)$membership
  cut_small_clusters(cl,3)
})
par(mfrow = c(2,2), mar = c(4,4,2,2))
par(oma = c(0, 0, 4, 0))    # outer margins (top space for title)

for(i in 1:4){
  c = as.integer(factor(cls[[i]]))-1
  hist(c, breaks = seq(0,length(table(c)))-0.5, xlab =  "clusters",ylab =  "cluster size", main = abs[i])
}
mtext("Cluster sizes for stable clustering without weights \nwith resolution = 1e-5",
      outer = TRUE,
      cex = 1.3,
      line = 1)
for(i in 1:4){print(max(table(cls[[i]]))/length(cls[[i]]))}
par(mfrow = c(2,1), mar = c(4,4,2,2))
par(oma = c(0, 0, 0, 0))
l = list(res1_nw,res1_w)
tt = c("Clustering without weights", "Clustering with weights")
for(j in 1:2){
  res1 = l[[j]]
  plot(log10(res1[1,,1]),res1[2,,1],col = adjustcolor(1,0.7), ylim = c(0.3,1),pch = 16, ylab =  "ARI",xlab=  "",cex = 1.2, main = tt[j])
  if(j == 2){
    mtext( "log10(resolution)",side = 1, line = 3)
    legend( "bottomright",legend = abs,col = seq(1,4),pch = 16)
    }
  for(i in 1:4){
    points(log10(res1[1,,i]),res1[2,,i],col = adjustcolor(i,0.7),pch = 16,cex = 1.2)
  }
  # for(i in 1:4){
  #   points(log10(res1[1,,i]),res1[5,,i],col = i,pch = 3)
  # }
 
}
# Distances with resolutions from 1e-5 to 1e-1 ----

# asses ARI with CPM as metric ----
## leiden clustering without weights ----
plan(multisession(workers = 20))
minres_nw = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-logw.RData"))
  edge_density(lg)
})
reso_nw = sapply(minres_nw,\(r){
  seq(r/2,0.2,by = (0.2-r/2)/100)
})

res_cpm_nw = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-l.RData"))
   
  source("cl_cores_functions.R")
  library(igraph)
  leidcl_stats(lg,reso_nw[,ab])
  
},simplify =  "array")
res = res_cpm_nw
x = reso
plot(x[,1],y=res[1,,1],ylim = c(0,1),xlim = (range(reso)),col = 1,main =  "CPM no weights",ylab ="ARI",xlab = "resolution")
for(i in 2:4){
  points(x[,i],res[1,,i],col = i)
}

legend( "bottomright",legend = abs,col = seq(1,4),pch = 1)

plot(x[,1],res[3,,1],col = 1,main = "nmb of clusters",ylim = c(1,max(res[3,,])))
for(i in 1:4){
  points(x[,i],res[3,,i],col = i,)
}
legend( "right",legend = abs,col = seq(1,4),pch = 1)

## with weights ----

minres_w = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-logw.RData"))
  ew = E(lg)$weight/max(E(lg)$weight)
  sum(ew)/(vcount(lg)*(vcount(lg)-1)/2)
})
reso = sapply(minres_w,\(r){
  seq(r/2,0.2,by = (0.2-r/2)/100)
})
plan(multisession(workers = 20))
#reso = seq(0.01,0.1,by = 0.01)
res_cpm_w = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-logw.RData"))
  E(lg)$weight = E(lg)$weight/max(E(lg)$weight)
  
  leidcl_stats(lg,reso[,ab])
  
},simplify =  "array")

res = res_cpm_w
x = reso
plot(x[,1],y=res[1,,1],ylim = c(0,1),xlim = (range(reso)),col = 1,main =  "CPM weights",ylab ="ARI",xlab = "resolution")
for(i in 2:4){
  points(x[,i],res[1,,i],col = i)
}
for(i in 1:4){ # fraction of saved vertices
  points(x[,i],res[4,,i],col = i,pch = 3)
}
legend( "bottomright",legend = abs,col = seq(1,4),pch = 1)

par(mfrow = c(1,1))
plot(x[,1],res[3,,1],col = 1,main = "Number of clusters",ylim = c(1,max(res[3,,])),)
for(i in 1:4){
  points(x[,i],res[3,,i],col = i,)
}
legend( "right",legend = abs,col = seq(1,4),pch = 1)

plot(x[,1],(res[4,,1]/max(res[4,,1]))*res[1,,1],col = 1,main = "ari vs clusters",xlim = (range(reso)),ylim = c(0,1))
for(i in 1:4){
  points(x[,i],(res[4,,i]/max(res[4,,i]))*res[1,,i],col = i)
}
## rnd_w ----
res_cpm_rw = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-logw.RData"))
  source("cl_cores_functions.R")
  library(igraph)
  E(lg)$weight = E(lg)$weight/max(E(lg)$weight)
  res_5r = sapply(1:5,\(rn){
    E(lg)$weight = sample(E(lg)$weight,ecount(lg))
    leidcl_stats(lg,reso[,ab])
  
  },simplify =  "array")
  apply(res_5r,c(1,2),mean)
},simplify =  "array")

res = res_cpm_rw
x = reso
plot(x[,1],y=res[1,,1],ylim = c(0,1),xlim = (range(reso)),col = 1,main =  "CPM rnd weights",ylab ="ARI",xlab = "resolution")
for(i in 2:4){
  points(x[,i],res[1,,i],col = i)
}
for(i in 1:4){ # fraction of saved vertices
  points(x[,i],res[4,,i],col = i,pch = 3)
}
legend( "bottomright",legend = abs,col = seq(1,4),pch = 1)

plot(x[,1],res[3,,1],col = 1,main = "nmb of clusters",ylim = c(1,max(res[3,,])))
for(i in 1:4){
  points(x[,i],res[3,,i],col = i,)
}
legend( "right",legend = abs,col = seq(1,4),pch = 1)

# check cluster distances ----
## with weights ----
reso2 = c(0.02,0.04,0.06)
plan(multisession(workers = 20))
dist_w = lapply(abs,\(ab){
  library(universalmotif)
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-logw.RData"))
  E(lg)$weight = E(lg)$weight/max(E(lg)$weight)
  
  lapply(reso2,\(r){ # change resolution 
    
    cls = cluster_leiden(lg,resolution = r, objective_function = "CPM",n_iterations = 10)
    m = 3
    lcssets = clust_big(V(lg)$name,memb = cls$membership,mins = m)
    lcssets5 = lapply(lcssets,\(set){
      c(unlist(sapply(set,\(pep){
        if(nchar(pep)==6){
          sapply(1:6,\(i) paste0(substr(pep,1,(i-1)),substr(pep,(i+1),6)))
        }
        else{
          pep
        }
      })))
    })
    leiden_mots = lapply(lcssets5,\(set)  create_motif(set,alphabet = "AA",pseudocount = 1)@motif)
    nn = length(leiden_mots)
    dm = unlist(future_sapply(1:(nn-1),\(ii){ # between random motifs
      sapply((ii+1):nn, \(jj){
        m1 = leiden_mots[[ii]]
        m2 = leiden_mots[[jj]]
        sqrt(sum((m1-m2)**2))
      })
    }))
    dm
  })
  
})
for(i in 1:4){
  boxplot(dist_w[[i]],main = abs[i])
}
## with random weights ----
reso3 = sapply(minres_w,\(r){
  seq(r/2,0.1,by = (0.1-r/2)/5)
})
plan(multisession(workers = 20))
dist_rw = lapply(abs,\(ab){
  library(universalmotif)
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-logw.RData"))
  E(lg)$weight = E(lg)$weight/max(E(lg)$weight)
  E(lg)$weight = sample(E(lg)$weight,ecount(lg))
  print(ab)
  lapply(reso2,\(r){ # change resolution 
    
    cls = cluster_leiden(lg,resolution = r, objective_function = "CPM",n_iterations = 10)
    m = 3
    lcssets = clust_big(V(lg)$name,memb = cls$membership,mins = m)
    lcssets5 = lapply(lcssets,\(set){
      c(unlist(sapply(set,\(pep){
        if(nchar(pep)==6){
          sapply(1:6,\(i) paste0(substr(pep,1,(i-1)),substr(pep,(i+1),6)))
        }
        else{
          pep
        }
      })))
    })
    leiden_mots = lapply(lcssets5,\(set)  create_motif(set,alphabet = "AA",pseudocount = 1)@motif)
    nn = length(leiden_mots)
    dm = unlist(future_sapply(1:(nn-1),\(ii){ # between random motifs
      sapply((ii+1):nn, \(jj){
        m1 = leiden_mots[[ii]]
        m2 = leiden_mots[[jj]]
        sqrt(sum((m1-m2)**2))
      })
    }))
    dm
  })
  
})
for(i in 1:4){
  boxplot(dist_rw[[i]],main = abs[i])
}

##  no weights ----
reso4 = sapply(minres_nw,\(r){
  seq(r/2,0.1,by = (0.1-r/2)/5)
})
reso4 = reso3
plan(multisession(workers = 20))
dist_nw = lapply(abs,\(ab){
  library(universalmotif)
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-l.RData"))
  print(ab)
  lapply(reso2,\(r){ # change resolution 
    
    cls = cluster_leiden(lg,resolution = r, objective_function = "CPM",n_iterations = 10)
    m = 3
    lcssets = clust_big(V(lg)$name,memb = cls$membership,mins = m)
    lcssets5 = lapply(lcssets,\(set){
      c(unlist(sapply(set,\(pep){
        if(nchar(pep)==6){
          sapply(1:6,\(i) paste0(substr(pep,1,(i-1)),substr(pep,(i+1),6)))
        }
        else{
          pep
        }
      })))
    })
    leiden_mots = lapply(lcssets5,\(set)  create_motif(set,alphabet = "AA",pseudocount = 1)@motif)
    nn = length(leiden_mots)
    dm = unlist(future_sapply(1:(nn-1),\(ii){ # between random motifs
      # print(ii)
      sapply((ii+1):nn, \(jj){
        m1 = leiden_mots[[ii]]
        m2 = leiden_mots[[jj]]
        
        sqrt(sum((m1-m2)**2))
      })
    }))
    dm
  })
  
})

## compare distance distributions----

for(i in 1:4){
  boxplot(dist_w[[i]],xlim = c(0,0.1),at = reso2,boxwex = 0.005,main = abs[i],names = reso2)
  boxplot(dist_rw[[i]],add = T,xlim = c(0,0.1),at = reso2,boxwex = 0.005,col = adjustcolor("pink",0.6),names =  F)
  boxplot(dist_nw[[i]],add = T,xlim = c(0,0.1),at = reso2,boxwex = 0.005,col = adjustcolor("yellow",0.6),names =  F)
  
}

### p values ----
allp = lapply(1:4,\(i){
  w_vs_nw = sapply(1:3,\(j){
    wilcox.test(dist_w[[i]][[j]],dist_nw[[i]][[j]],alternative =  "greater")$p.value
  })
  w_vs_rw = sapply(1:3,\(j){
    wilcox.test(dist_w[[i]][[j]],dist_rw[[i]][[j]],alternative =  "greater")$p.value
  })
  rw_vs_nw = sapply(1:3,\(j){
    wilcox.test(dist_rw[[i]][[j]],dist_nw[[i]][[j]],alternative =  "greater")$p.value
  })
  comp = cbind(w_vs_nw,w_vs_rw,rw_vs_nw)
  comp
})
sapply(allp,p.adjust,method =  "BH")
### plots ----

par(mfrow = c(2, 2), mar = c(2,2,2,2))
par(oma = c(2, 2, 4, 0))
for(i in 1:4){
  par(xpd = NA, bty = "n")
  plot(1, type = "n",bty = "n", xlab = "", ylab = "",ylim = range(unlist(dist_w)),xlim = c(0,10),axes = F, frame.plot = FALSE,main = abs[i])
  if(i>2){mtext("resolution", side = 1, line = 3)}
  if(i%%2==1){
    mtext("Inter-motif distances",side = 2,line = 3)
    
    }
  offset = c(0,0.5,1)
  axis(2, labels = T) 
  axis(1,at = seq(0,2)*3 +1+offset,labels = round(reso2,2),tick = FALSE, lwd = 0)
  boxplot(dist_w[[i]],add = T,bty =  "n",notch = T,at = seq(0,2)*3 +offset,boxwex = 0.5,names =  F,axes = F, outline = F)
  boxplot(dist_rw[[i]],add = T,bty =  "n", notch = T,xlim = c(0,0.1),at = seq(0,2)*3 +1 + offset,boxwex = 0.5,col = adjustcolor("pink",0.6),names =  F,axes = F, outline = F)
  boxplot(dist_nw[[i]],add = T,bty =  "n", notch = T,xlim = c(0,0.1),at = seq(0,2)*3 +2 + offset,boxwex = 0.5,col = adjustcolor("yellow",0.6),names =  F,axes = F, outline = F)
  
  if(i==4){
    legend( "bottomright", legend = c( "weights", "random weights",  "no weights"),col = c( "grey",  "pink",  "yellow"), pch = 15, cex = 1.2,bty =  "n")
  }
}


## plot weights vs random weights----
par(mfrow = c(4, 1), mar = c(4,4,2,2))
res = res_cpm_w
x = reso
plot(x[,1],y=res[2,,1],pch = 16,bty = "n",ylim = c(0.5,1.05),xlim = (range(reso)),col = 1,main =  "Clustering stability, \n graphs with weights",ylab ="ARI",xlab = "")
for(i in 2:4){
  points(x[,i],res[2,,i],col = i, pch = 16)
}
#legend( "bottomright",legend = abs,col = seq(1,4),pch = 16,bty = "n")
legend( "right",legend = abs,col = seq(1,4),pch = 16,bty =  "n")
par(mar = c(4,4,2,2))
res = res_cpm_rw
x = reso
plot(x[,1],y=res[2,,1],pch = 16,bty = "n",ylim = c(0.5,1.05),xlim = (range(reso)),col = 1,main =  "Clustering stability, \n graphs with random weights",ylab ="ARI",xlab = "")
for(i in 2:4){
  points(x[,i],res[2,,i],pch = 16,col = i)
}
#legend( "bottomright",legend = abs,col = seq(1,4),pch = 16,bty = "n")
res = res_cpm_w
plot(x[,1],res[4,,1],pch = 16,col = 1,main = "Number of clusters, \n graphs with weights",ylim = c(1,max(res[4,,])),xlab =  "",ylab =  "clusters")
for(i in 1:4){
  points(x[,i],res[4,,i],col = i,pch = 16)
}
res = res_cpm_w
plot(x[,1],res[5,,1],pch = 16,col = 1,main = "Vertices in cluster with size>3, \n graphs with weights",ylim = c(0,1),xlab =  "resolution",ylab =  "proportion in clusters")
for(i in 1:4){
  points(x[,i],res[5,,i],col = i,pch = 16)
}

## test difference ----
for(i in 1:4){
  w = wilcox.test(res_cpm_w[2,25:100,i], res_cpm_rw[2,25:100,i], paired = T, alternative =  "greater")
  print(w)
}


# choose best parameters ----

# similarity to next 10 points
sim10 = sapply(abs,\(ab){
  sapply(1:90,\(i){
    a = res_cpm_w[2,i,ab]
    b = res_cpm_w[2,(i+1):(i+10),ab]
    abs((a-mean(b))/sd(b))
  })
})
bestreso_w = sapply(abs,\(ab){
  res = res_cpm_w[,1:90,ab]
  r = reso[1:90,ab]
  rb = r[which.max(res[4,])]
  k =which(res[2,]>0.99*max(res[2,]) & sim10[,ab]<3 & r>rb & res[4,]>0.6*max(res[4,]) & res[5,]>0.5)[1]
  print(res[1,k])
  r[k]
})
bestreso_w = round(bestreso_w,3)

## cluster with best parameters ----
leiden_best_lcs = lapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-logw.RData"))
  E(lg)$weight = E(lg)$weight/max(E(lg)$weight) # norm weights
  r = bestreso_w[ab]
  n = 10
  cls = sapply(1:n,\(i) cluster_leiden(lg,resolution = r, objective_function = "CPM",n_iterations = 10)$membership)
  m = 3
  cls2 = sapply(1:n,\(i) cut_small_clusters(cls[,i],m))
  ari = (future_sapply(1:(n),\(i){sapply((1):n,\(j) {adjustedRandIndex(cls2[,i],cls2[,j]) })})) # mean ARI, better overlap if high
  cl = cls[,which.max(colSums(ari))]
  m = 3
  lcssets = cl_cores(lg,cl,m,0)
  })
names(leiden_best_lcs) = abs
leid_lengths = lapply(leiden_best_lcs,\(sets) sapply(sets,length))
clnames = sapply(1:4,\(i) sapply(1:length(leiden_best_lcs[[i]]),\(j) paste(abs[[i]],j,sep = "_"))) # unique cluster names
for(i in 1:4){ names(leiden_best_lcs[[i]]) = clnames[[i]]  }

leiden_lcssets5 = lapply(leiden_best_lcs,\(absets){
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
names(leiden_lcssets5) = abs

plan(multisession(workers = 4))
leiden_pepsets = future_lapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(file = paste0(path,ab,"_pep_lcs.RData"))
  lcs_to_peps(leiden_lcssets5[[ab]],pep_lcs)
})
names(leiden_pepsets) = abs
for(i in 1:4){ names(leiden_pepsets[[i]]) = clnames[[i]]  }

save(leiden_best_lcs,bestreso_w,leid_lengths,leiden_lcssets,leiden_lcssets5,leiden_pepsets,file = ("leiden_clustering_result_090126_bestreso_weights.RData"))

### prune small clusters ----

pepsets_len = lapply(leiden_pepsets,\(sets) sapply(sets,length)); names(pepsets_len) = abs
leiden_pepsets_2 = lapply(abs,\(ab){
  leiden_pepsets[[ab]][pepsets_len[[ab]]>1]
})
leiden_pepsets_3 = lapply(abs,\(ab){
  leiden_pepsets[[ab]][pepsets_len[[ab]]>2]
})
names(leiden_pepsets_2) = abs
names(leiden_pepsets_3) = abs

### remove fully overlapping clusters ----
leiden_pepsets_un = future_lapply(leiden_pepsets_3,\(sets){
  n = length(sets)
  to_remove = c()
  for(s1 in 1:n){
    for(s2 in 1:n){
      is12 = length(intersect(sets[[s1]],sets[[s2]]))
      if(s1 != s2 & is12 == length(sets[[s1]])){
        to_remove = c(to_remove,s1)
        break
      }
    }
  }
  sets[-to_remove]
})

save(leiden_pepsets_un,file =  "unique_leiden_pepsets.RData")
### plot kld ----
len = lapply(leiden_pepsets_un,\(sets) sapply(sets, length))
kld = function(p,q){
  sum(p*log(p/q))
}
leid_ppm = lapply(leiden_pepsets_un,\(sets){ lapply(sets,\(cl) create_motif(cl, alphabet = "AA", pseudocount = 1)@motif)})
leid_kld = lapply(leid_ppm,\(sets){ sapply(sets, kld, bgmot)})

load(file = "rnd_pepsets_libfreqs_leiden.RData")
rndlen = sapply(rnd_pepsets,length)
rnd_pepsets_3 = rnd_pepsets[rndlen>2]
cntr_len = rndlen[rndlen>2]
control_ppm = lapply(rnd_pepsets_3,\(peps) create_motif(peps, alphabet = "AA", pseudocount = 1)@motif)
control_kld = sapply(control_ppm, kld, bgmot)
load(file="rnd_len_kld_2.RData")
rnd_s = rnd_len_kld[rnd_len_kld[1,]<=1.2*max(unlist(len)),]
control_kld = c(control_kld,rnd_s[,2])
cntr_len = c(cntr_len, rnd_s[,1])

plot(log2(cntr_len[cntr_len<130]),log2(control_kld[cntr_len<130]),ylim = c(0,log2(max(unlist(leid_kld)))),col = adjustcolor("grey",0.7),pch = 16)
for(i in 1:4){
  points(log2(len[[i]]),log2(leid_kld[[i]]) , col = colors[i],xlim = c(0,1.2*max(unlist(len))))
}

wilcox.test(control_kld[cntr_len<10], unlist(leid_kld), alternative = "less")
boxplot(list(control_kld[cntr_len<10], unlist(leid_kld)), axes = F,ylim = c(min(control_kld[cntr_len<10])-2,max(unlist(leid_kld))+2),ylab = "KLD",main = "KLD compared to random clusters")
axis(2)
axis(1, labels = c("random","real"),at = c(1,2),tick = FALSE, lwd = 0)
points(1.5,20,pch = 8)
lines(c(1,2),c(19,19))

## plot largest motifs 
for(i in 1:4){
  ii = which.max(sapply(leiden_pepsets_un[[i]], length))
  logoal(leiden_pepsets_un[[i]][[ii]])
  print(length(leiden_pepsets_un[[i]][[ii]]))
}
### also for lcssets
leiden_lcssets_un = future_lapply(leiden_lcssets5,\(sets){
  n = length(sets)
  to_remove = c()
  for(s1 in 1:n){
    for(s2 in 1:n){
      is12 = length(intersect(sets[[s1]],sets[[s2]]))
      if(s1 != s2 & is12 == length(sets[[s1]])){
        to_remove = c(to_remove,s1)
        break
      }
    }
  }
  sets[-to_remove]
})
leiden_lcsmots = lapply(leiden_lcssets_un,\(lcssets){
  lapply(lcssets,\(set)  create_motif(set,alphabet = "AA",pseudocount = 1)@motif)
})
save(leiden_lcssets_un,leiden_lcsmots,file = "unique_leiden_lcssets.RData")  # with 4 replicates of b12
save(leiden_lcssets_un,leiden_lcsmots,file = "unique_leiden_lcssets_m.RData")  # with 3 replicates of b12, in thesis
