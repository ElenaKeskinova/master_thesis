load("mixed-7graphs/lib_allpeps_7nC.RData")
library(igraph)
source("newgraph.R")
source("graph_to_line.R")
library(universalmotif)
source("cmd_add.R")
library(uwot)

# plot library motifs over map of rnd motifs ----
load("20k_lib_lcsmots.RData")
load(file = "20K_rndmotifs.RData")
load(file =  "cmd_rndmotifs20k_l.RData")

smots = sample(1:length(lib_motifs),size = 3000)
smots = sample((lib_motifs),size = 3000)
lib_dst = future_sapply(rnd_motifs,\(m1){ # between library and random motifs
  sapply(lib_motifs, \(m2){
    sqrt(sum((m1-m2)**2))
  })
})

libmot_y = cmd_add(D_train = distm_r,d_new = lib_dst,V = coord_r$points[,1:95], coord_r$eig[1:95])
cdu = umap(coord_r$points[,1:95],n_components = 2, ret_model = T)
lib_umap = umap_transform(X = libmot_y,model = cdu)

par(mfrow = c(1,1), mar = c(3,3,1,1))
plot(cdu$embedding[,1:2], pch = 16, col =  adjustcolor( "orange",0.5),xlab = "dimension 1",ylab = "dimension 2",cex = 0.5, ylim = range(lib_umap[,2]) + c(-0.5,0.5),xlim = range(lib_umap[,1]))
points(lib_umap, pch = 16,col =  adjustcolor( "black",0.6),cex = 0.5)
legend("topleft",legend = c("random","from library"),pch = 16, col = c("orange","black"),bty = "n", cex = 1.3)


cdd = umap(rbind(coord$points[,1:95],libmot_y),n_components = 2)
plot(cdd[1:20000,], pch = 16, col =  adjustcolor( "orange",0.5),cex = 0.8,xlab = "dimension 1",ylab = "dimension 2", xlim = range(cdd[,1]),ylim = range(cdd[,2]))
points(cdd[20001:nrow(cdd),], pch = 16,col =  adjustcolor( "black",0.6),cex = 0.8)
legend("topleft",legend = c("random","from library"),pch = 16, col = c("orange","black"),bty = "n")

# calculate cmdscale of library motifs ----

library(future.apply)
plan(multisession(workers = 20))
distm = future_sapply(lib_motifs,\(m1){ # between random library motifs
  sapply(lib_motifs,\(m2){
    sqrt(sum((m1-m2)**2))
  })
})

save(distm,lib_motifs,file = "20Klib_mots.RData")

# only library
coord = cmdscale(distm, k = 100, eig = T)
save(coord,file =  "cmd_libmots20k.RData")
save(distm,lib_motifs,coord,file =  "cmd_libmotifs_20k.RData") # 13/01/26 leiden motifs with resolution 0.05

load(file = "cmd_libmotifs_20k.RData")
plot(coord$eig[1:100])
cd2 = umap(coord$points[,1:95],n_components = 2, ret_model = T)
plot(cd2$embedding)

## add real leiden motifs to the embedding of library motifs ----
plan(multisession(workers = 20))
load(file = "unique_leiden_lcssets_m.RData")
motifs = unlist(leiden_lcsmots,recursive = F)
names(motifs) = unlist(sapply(leiden_lcsmots,names)) # motifs from lcs clusters from the 4 graphs
d_new = future_sapply(lib_motifs,\(m1){ # between random motifs
  sapply(motifs, \(m2){
    sqrt(sum((m1-m2)**2))
  })
})

mot_y = cmd_add(D_train = distm,d_new = d_new,V = coord$points[,1:95], coord$eig[1:95])
cdmot = umap_transform(X = mot_y[,1:95],model = cd2)


save(d_new,mot_y,file = "cmd_libmots20k_lcsmotcoords.RData") # 13/01/26 with four replicates for b12
save(d_new,mot_y,file = "cmd_libmots20k_lcsmotcoords_m.RData") # 04/03/26



codes_nm = sapply(names(motifs),\(nm) substr(nm,1,3))
numcodes = unlist(sapply(1:4,\(i) rep(i,length(leiden_lcsmots[[i]]))))

plot(cd2$embedding, pch = 16, col =  adjustcolor( "grey",0.5))
points(cdmot, pch = 16,col =  adjustcolor( "black",0.7))

plot(cd2$embedding, pch = 16, col =  adjustcolor( "grey",0.5))

# plots with best motifs in red ----
load(file = "sim_score_pvals.RData")
load(file = "similarity_scores_leiden.RData")
par(mfrow = c(2,2))
bestmots = sapply(1:4,\(i) which(max_scores_all[which(codes==i),i] == max(max_scores_all[which(codes==i),i]), arr.ind = T))
largest = sapply(leiden_pepsets_un,\(set) names(which.max(sapply(set,length))))
for(i in 1:4){
  plot(cd2$embedding[,1:2], pch = 16, col =  adjustcolor( "grey",0.4), xlab = "",ylab = "",xlim = range(cd2$embedding[,1:2])+c(0,0.5))
  title(main = abs[i],cex.main = 1.7)
  if(i%%2==1){
  title(ylab = "dimension 2",cex.lab = 1.2)
  }
  if(i>2){
    title(xlab = "dimension 1",cex.lab = 1.2)
  }
  pnts = cdmot[numcodes == i,1:2]
  points(pnts, pch = 16,col =  adjustcolor("black" ,0.6),cex = 1)
  points(pnts[names(which(p[[i]]<0.05)),], pch = 16,col =  adjustcolor("orangered" ,1),cex = 1.2)
  points(pnts[names(bestmots[i]),1],pnts[names(bestmots[i]),2], pch = 16,col =  adjustcolor("orangered" ,1),cex = 2.2)
  #points(pnts[names(which(phyper_adj[[i]]<=0.05)),], pch = 0,col =  adjustcolor("orange" ,1),cex = 2.2,lwd = 2.5)
  #points(pnts[largest[i],1],pnts[largest[i],2], pch = 0,col =  adjustcolor("green" ,1),cex = 2.2,lwd = 2.5)
  
  #legend("bottomright", col = c("grey","black","orangered","orangered","orange"),pch = c(16,16,16,16,0),bty = "n",legend = c("library","antibody clusters","z-score < 0.05","closest to epitope","significantly \nidiotypic"),pt.cex = c(1.2,1.2,1.2,2.2,2.2))
  # points(lib_umap, pch = 16,col =  adjustcolor( "red",0.3))
  
}

plot.new()
legend("center", col = c("grey","black","orangered","orangered"),pch = c(16,16,16,16),bty = "n",legend = c("library","antibody clusters","epitope similarity \np-value< 0.05","closest to epitope"),cex = 1.2,pt.cex = c(1.2,1.2,1.2,2.2))

# see neighbour counts vs similarity to epitope ----

dist1 = apply(d_new,1,min)

dist100 = sapply(1:nrow(d_new),\(di) sort(d_new[di,],decreasing = F)[100])
hist(dist100)
names(dist100) = names(motifs)
codes = unlist(sapply(1:4,\(i) rep(i,length(leiden_pepsets_un[[i]]))))
par(mfrow = c(2,2),mar = c(3,3,2,2))
for(i in 1:4){
  plot(dist100[rownames(max_scores_all)[codes==i]],max_scores_all[codes==i,i])
}

hist(d_new)
hist(sample(distm,5000))

rd = quantile(dist1,0.99)
#rd = max(dist1)
pntsd = rowSums(d_new<=rd)
for(i in 1:4){
  plot(log2(pntsd[rownames(max_scores_all)[codes==i]]),max_scores_all[codes==i,i], xlab = "neighbours", ylab = "score", main = abs[i])
}

par(mfrow = c(1,4))
for(i in 1:4){
  scores = max_scores_all[codes==i,i]
  pn = log10(pntsd[rownames(max_scores_all)[codes==i]]+1)
  igood = which(scores>12)
  ibad = which(scores<=12)
  if(i==3){igood = which(scores>12)}
  boxplot(list(pn[igood],pn[ibad]), boxwex = 0.6,ylim = log10(range(pntsd+1)), main = abs[i], outline = F, col = 0,axes = F)
  points(1+rnorm(length(igood),0,0.02),pn[igood])
  axis(1, labels = c("score>12","score<=12"),at = c(1,2),tick = FALSE, lwd = 1, cex.axis = 1.1)
  
  pw = wilcox.test(pn[igood],pn[ibad],alternative = "less")
  print(pw)
  lines(x=c(1,2),y=c(3.3,3.3))
  if(pw$p.value<0.05){
    points(1.5,3.5,pch = 8)
  }
  else{
    text(1.5,3.5, labels = "n.s.")
  }
  if(i==1){
    mtext("log10(number of close library motifs+1)",side = 2, line = 3)
    axis(2)
  }
}

