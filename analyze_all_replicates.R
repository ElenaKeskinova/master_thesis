abs = c("Herceptin","21c","17b","b12")

for(ab in abs){
  path = paste0("Exp12.f_UNI-10/", ab,"/")
  files <- list.files(path = path, pattern = "uniqueP_[[:upper:]]{5}.txt", full.names = TRUE)
  path8 = paste0("Exp12.f_UNI-8/", ab,"/")
  files8 <- list.files(path = path8, pattern = "uniqueP_[[:upper:]]{5}.txt", full.names = TRUE)
  
  fnames = sapply(files, \(f){
    n = unlist(strsplit(f,"uniqueP_"))[2]
    substr(n,1,5)
  })
  names(files) = fnames
  names(files8) = fnames
  
  #get all peps together for length 10
  AllPep = lapply(files,\(f) read.table(file = f,
                                        colClasses = c("character","integer"),col.names = c("Pep","Freq")))
  
  
  Peps_noC = lapply(AllPep,\(df){
    noC = sapply(df$Pep, function(p) {p1 = strsplit(p,"")[[1]]; p1[1]!="C" })
    df[noC,]
  })
  
  
  
  ### for 10 letters ----
  
  Peps7 = lapply(Peps_noC, \(dfp){
    pp = sapply(dfp$Pep, \(p){
      ps = strsplit(p,"")[[1]]
      p1 = paste0(ps[1:7],collapse = "")
      p2 = paste0(ps[4:10],collapse = "")
      return(c(p1,p2))
    })
    c(pp)
  })
  
  
  freq10to7 = lapply(Peps_noC,\(df){rep(df$Freq,each = 2)})
  
  ### for 8 letters ----
  
  AllPep = lapply(files8,\(f) read.table(file = f, colClasses = c("character","integer"),col.names = c("Pep","Freq")))
  Peps_noC8 = lapply(AllPep,\(df){
    noC = sapply(df$Pep, function(p) {p1 = strsplit(p,"")[[1]]; p1[1]!="C" })
    df[noC,]
  })
  
  Peps87 = lapply(Peps_noC8, \(dfp){
    sapply(dfp$Pep,\(p){
      p1 = strsplit(p,"")[[1]]
      p1 = paste0(p1[1:7],collapse = "")
      p1
    })
  })
  
  freq8to7 = lapply(Peps_noC8,\(df)df$Freq)
  
  ### merge both dataframes ----
  peps = lapply(fnames,\(fn) c(Peps7[[fn]],Peps87[[fn]]))
  names(peps) = fnames
  freqs = lapply(fnames,\(fn) c(freq10to7[[fn]],freq8to7[[fn]]))
  names(freqs) = fnames
  
  all7 = lapply(fnames,\(fn) aggregate(freqs[[fn]],by = list(Pep = peps[[fn]]),FUN = sum)) # 4 data frames with freqs
  names(all7) = fnames
  
  #write file with all peps and frequencies
  pathin = paste0("replicate_graphs/")
  save(all7, file =  paste0(pathin,ab,"_repl_pepfreqs.RData"))
  Fin = paste0(pathin,ab,"_",fnames,"_allp7-fr.txt")
  for(i in 1:4){
    for (ii in 1:length(all7[[i]][[1]])){
      cat(sprintf("%s\n", paste0(all7[[i]]$Pep[ii], "    ", all7[[i]]$x[ii])), file = Fin[i], append = TRUE)
      
    }
  }
  
  #only peps without frequencies
  Fin = paste0(pathin,ab,"_",fnames,"_allp7.txt")
  for(i in 1:4){
    
    p =all7[[i]]$Pep
    cat(p,file = Fin[i], sep = "\n")
  }
  
  closeAllConnections()
  
}

# put all peps in one big data frame with clusters of big graphs occurences in every replicate ----
load(file =  "unique_leiden_pepsets.RData")

proportions = lapply(abs,\(ab){ # proportions of vertices from each replicate in the big graph
  load(file =  paste0( "replicate_graphs/",ab,"_repl_pepfreqs.RData"))
  reppeps   = sapply(all7,\(l) (l$Pep))
  path = paste0("mixed-7graphs/",ab,"/")
  load(file = paste0(path,ab,"big7or.RData"))
  gpeps = V(G)$name
  sapply(reppeps,\(rp) length(intersect(rp,gpeps))/length(gpeps))
})

# which cluster size has over 95% probability to contain all 4 replicates ----

minsize = sapply(proportions,\(p){
  p = min(p)
  ps = sapply(1:50,\(k){
    pbinom(0,k,p)
  })
  which(ps<0.05)[1]
})

# distribution of replicate peps per leiden cluster ----
cluster_distr = lapply(1:4,\(i){
  ab = abs[i]
  load(file =  paste0( "replicate_graphs/",ab,"_repl_pepfreqs.RData"))
  reppeps   = lapply(all7,\(l) l$Pep)
  cls = leiden_pepsets_un[[i]]
  m = sapply(reppeps,\(rp){
    sapply(cls,\(cp){
      length(intersect(cp,rp))
    }, simplify =  "array")
  })
})
cluster_distr10 = lapply(1:4,\(i){
  cluster_distr[[i]][which(sapply(leiden_pepsets_un[[i]],length)>=10),]
})
chisq.test(cluster_distr10[[1]])

cluster_distr_TF = lapply(1:4,\(i) {
  cluster_distr[[i]]>=sapply(leiden_pepsets_un[[i]],length)/4
})
llen = lapply(leiden_pepsets_un,\(ps) sapply(ps,length))
sapply(llen,median)

cluster_distr_ns = lapply(1:4,\(i){
  cld = cluster_distr[[i]]
  l = llen[[i]]
  sapply(1:4,\(j){
    sapply(1:length(l),\(ci){
      cld[ci,j] > qbinom(0.05,l[ci],proportions[[i]][j])
    })
  })
})
par(mfrow = c(2,2))
for(i in 1:4){
  plot(llen[[i]], rowSums(cluster_distr_ns[[i]]),ylab =  "represented number \nof replicates", xlab =  'cluster size', cex = 1.5 ,main = abs[i],pch = 16, col = adjustcolor( "grey",0.5))
}

t2 = lapply(leiden_pepsets_un,\(lp) table(unlist(lp)))

# cluster big graphs with leiden and modularity with res = 1 ----

## with all four replicates: ----
pepleid = lapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  if(ab== "b12"){
    load(file = paste0(path,ab,"big7or_old4.RData"))
  }
  else{load(file = paste0(path,ab,"big7or.RData"))}
  
  cluster_leiden(G, objective_function = "modularity")
})

pepleidm = sapply(pepleid, \(pl) pl$membership)
pepleid_pepsets = lapply(1:4,\(i){
  mm = pepleidm[[i]]
  ms = minsize[i]
  mmbig = which(table(mm)>=ms)
  ab = abs[i]
  path = paste0("mixed-7graphs/",ab,"/")
  if(ab== "b12"){
    load(file = paste0(path,ab,"big7or_old4.RData"))
  }
  else{load(file = paste0(path,ab,"big7or.RData"))}
  
  pepsets = lapply(mmbig,\(mi){
    V(G)$name[which(mm==mi)]
  })
})
cluster_distrp = lapply(1:4,\(i){
  ab = abs[i]
  load(file =  paste0( "replicate_graphs/",ab,"_repl_pepfreqs.RData"))
  reppeps   = lapply(all7,\(l) l$Pep)
  cls = pepleid_pepsets[[i]]
  m = sapply(reppeps,\(rp){
    sapply(cls,\(cp){
      length(intersect(cp,rp))
    }, simplify =  "array")
  })
})

library(chisq.posthoc.test)
chip = lapply(1:4,\(i){
  chisq.posthoc.test(cluster_distrp[[i]], simulate.p.value = T, method =  "BH")
})
chi_p = lapply(chip,\(df){
  df[seq(2,nrow(df),2),]
})
chi_r = lapply(chip,\(df){
  df[seq(1,nrow(df)-1,2),]
})
for(i in 1:4){heatmap(as.matrix(chi_r[[i]][,3:6]), main = abs[i])}

library(pheatmap)
for(i in 1:4){pheatmap(as.matrix((chi_r[[i]][,3:6])),color = colorRampPalette(c("blue","white","firebrick"))(11), breaks = c(-15,-10,-7,-5,-4,-3,3,4,5,7,10,15),legend = F,display_numbers = ifelse(chi_p[[i]][,3:6]<0.05, "*", ""),fontsize_number = 20,main = abs[i], scale =  "none", cluster_rows = F,treeheight_row = 0, cluster_cols = F, labels_col = c(1,2,3,4),labels_row = seq(1,nrow(chi_p[[i]])), xlab =  "replicates", ylab =  "clusters")}

sapply(chi_p, \(p) which(p<0.05,arr.ind = T))

# remove peps from first replicate of b12 where they are overrepresented in clusters -----
badcl4 = which(chi_p[[4]]$AGGTC<0.05 & chi_r[[4]]$AGGTC>3)
badpeps = unlist(pepleid_pepsets[[4]][badcl4])
ab = abs[4]
load(file =  paste0( "replicate_graphs/",ab,"_repl_pepfreqs.RData"))
reppeps4   = lapply(all7,\(l) l$Pep)[[1]] # replicate with most sequences
badpeps = intersect(badpeps,reppeps4)

path = paste0("mixed-7graphs/",ab,"/")
load(file = paste0(path,ab,"big7or.RData")) # graph with 3 replicates



i = 4

reppeps   = lapply(all7[2:4],\(l) l$Pep)
gpeps = V(G)$name
props_new = sapply(reppeps,\(rp) length(intersect(rp,gpeps))/length(gpeps))
p = min(props_new)
ps = sapply(1:50,\(k){
  pbinom(0,k,p)
})
mins = which(ps<0.05)[1]

#ms = minsize[i]
cl = cluster_leiden(G, objective_function = "modularity")
mm = cl$membership
mmbig = which(table(mm)>=mins)
pepsets4 = lapply(mmbig,\(mi){
  V(G)$name[which(mm==mi)]
})



distr4 =  sapply(reppeps,\(rp){
  sapply(pepsets4,\(cp){
    length(intersect(cp,rp))
  }, simplify =  "array")
})

chi4 = chisq.posthoc.test(distr4, simulate.p.value = T, method =  "BH")
chi_p4 = chi4[seq(2,nrow(chi4),2),3:5]

chi_r4 = chi4[seq(1,nrow(chi4)-1,2),3:5]


library(pheatmap)
pheatmap(as.matrix((chi_r4)),color = colorRampPalette(c("blue","white","firebrick"))(11), breaks = c(-15,-10,-7,-5,-4,-3,3,4,5,7,10,15),display_numbers = ifelse(chi_p4<0.05, "*", ""),fontsize_number = 20,main =  "b12 clusters \nwithout replicate 1", scale =  "none", cluster_rows = F,treeheight_row = 0, cluster_cols = F, labels_col = c(2,3,4),labels_row = seq(1,nrow(chi_p4)), xlab =  "replicates", ylab =  "clusters")

for(i in 1:4){
  logoal(all7[[i]]$Pep)
}
