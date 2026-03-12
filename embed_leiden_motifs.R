# generate rnd motifs and embed them with multidimensional scaling
library(uwot)
library(future.apply)
require(ggseqlogo)
require(ggplot2)
require(Biostrings)
source("printlogo.R")


# rnd motifs with distributions like real motifs and like library ----
load(file = "unique_leiden_lcssets.RData")
load(file = "20k_lib_lcsmots.RData")

## 10 000 motifs with identical distributions in columns
motifs = unlist(leiden_lcsmots,recursive = F)
lmotifs = sample(lib_motifs,length(motifs))
motcols = do.call(cbind,c(motifs,lmotifs))
aanm = rownames(motcols)
nc = ncol(motcols)

rnd_motifs4 = unlist(lapply(1:1000,\(i){
  motc1 = motcols[sample(1:20,20),]
  rownames(motc1) = aanm
  
  lapply(1:10,\(j){
    mot = motc1[,sample(1:nc,5)]
    #mot[mot<0.2] = 0
    #t(t(mot)/colSums(mot))
  })
}),recursive = F)

# 5000 motifs with reduced entropy, only letters>0.2
lmotifs = sample(lib_motifs,length(motifs))
motcols = do.call(cbind,c(motifs,lmotifs))
rnd_motifs5 = unlist(lapply(1:1000,\(i){
  motc1 = motcols[sample(1:20,20),]
  rownames(motc1) = aanm
  
  lapply(1:5,\(j){
    mot = motc1[,sample(1:nc,5)]
    mot[mot<0.2] = min(mot)
    t(t(mot)/colSums(mot))
  })
}),recursive = F)

# 5000 motifs with reduced entropy, only letters>0.1
lmotifs = sample(lib_motifs,length(motifs))
motcols = do.call(cbind,c(motifs,lmotifs))
rnd_motifs6 = unlist(lapply(1:1000,\(i){
  motc1 = motcols[sample(1:20,20),]
  rownames(motc1) = aanm
  
  lapply(1:5,\(j){
    mot = motc1[,sample(1:nc,5)]
    mot[mot<0.1] = min(mot)
    t(t(mot)/colSums(mot))
  })
}),recursive = F)

rnd_motifs = c(rnd_motifs4,rnd_motifs5,rnd_motifs6)


## distance matrix between all random motifs ----
library(future.apply)
plan(multisession(workers = 20))
distm_r = future_sapply(rnd_motifs,\(m1){ # between random motifs
  sapply(rnd_motifs,\(m2){
    sqrt(sum((m1-m2)**2))
  })
},future.chunk.size = 1000)

save(distm_r,rnd_motifs,file = "20K_rndmotifs.RData") # with motifs from library and from antibodies, 12/01/2026


# only random
coord_r = cmdscale(distm_r, k = 100, eig = T)
save(coord_r,distm_r,rnd_motifs,file =  "cmd_rndmotifs20k_l.RData") # with motifs from library, 12/01/2026


load(file = "20K_rndmotifs.RData")
load(file =  "cmd_rndmotifs20k_l.RData")

# d_rnd = distm[(length(motifs)+1):length(allmot),(length(motifs)+1):length(allmot)]
# mx1 = distm[1:length(motifs),(length(motifs)+1):length(allmot)] # ori motifs vs random
# mx2 = coord$points # rnd motifs embedding
# motc= mx1%*%mx2
distmot = future_sapply(motifs,\(m1){ # between  motifs
  sapply(motifs, \(m2){
    sqrt(sum((m1-m2)**2))
  })
})

plan(multisession(workers = 20))
d_new = future_sapply(rnd_mots,\(m1){ # between random motifs
  sapply(motifs, \(m2){
    sqrt(sum((m1-m2)**2))
  })
})

mot_y = cmd_add(D_train = distm,d_new = d_new,V = coord$points[,1:95], coord$eig[1:95])

save(distmot,d_new,mot_y,file = "embed_abmots_20krnd.RData")
barplot(apply(coord$points[,1:95],2,var))

dist2 = apply(mx2[1:50,],1,\(v1) apply(mot_y[,1:95],1,\(v2){ sqrt(sum((v1-v2)**2))}))
d_dist = dist2 - mx1[,1:50]

cd2 = umap(t(t(coord$points[,1:95])/coord$eig[1:95]), ret_model = T)
cd2 = umap(coord$points[,1:95],n_components = 2, ret_model = T)
cdmot = umap_transform(X = mot_y[,1:95],model = cd2)
cdmot2 = umap_transform(X = motc,model = cd2)

library(rgl)
plot3d(cd2$embedding)
plot(cd2$embedding, pch = 16, col =  "grey")

points(cdmot, pch = 16,col =  "black")
points(cd2[(length(motifs)+1):length(allmot),], pch = 16,col =  "pink")
points(cd2[1:length(motifs),], pch = 16,col =  "black")

# kmeans

kmn_cls = kmeans(cd2$embedding, centers = 20)
for(i in 1:20){
  
  plot(cd2$embedding, pch = 16, col =   "grey")
  points(cd2$embedding[which(kmn_cls$cluster==i),],pch = 16,col =  "pink")
}

for(i in 1:20){
  cl1 = which(kmn_cls$cluster==i)
  cl1m = rnd_motifs[cl1]
  ori_order = rownames(motifs[[1]])
  cl1m = lapply(cl1m,\(m) {rownames(m) = ori_order; m})
  
  m = Reduce('+', cl1m)/length(cl1m)
  print(ggseqlogo(m))
}


save(K_new,V,lambda,y_new,file = "matrix_coord.RData")

