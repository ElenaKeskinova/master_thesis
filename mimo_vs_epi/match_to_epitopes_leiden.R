load(file = "epitope_graphs_6.RData")
source("mimo_vs_epi/paths_scores_functions.R")
source("cl_cores.R")
load(file = "lib_freqs.RData")
bgmot = freqM7
load( "leiden_lcssets.RData")
abs = c("Herceptin","21c","17b","b12")

library(universalmotif)
library(igraph)
library(future.apply)

# load pepsets from leiden clusters ----
load(file =  "unique_leiden_pepsets.RData")

plan(multisession(workers = 20))
leid_pssms = future_lapply(leiden_pepsets_un,\(pepsets){
  lapply(pepsets,\(set){
    mot = create_motif(set,type = "PWM",alphabet = "AA",bkg = rowSums(freqM7)/7, pseudocount = 1)@motif
    mot = rbind(mot,rep(0,ncol(mot)))
    rownames(mot)[nrow(mot)] = "-"
    mot
  })
}) # with pseudocount

codes = unlist(sapply(1:4,\(i) rep(i,length(leiden_pepsets_un[[i]]))))

require(ggseqlogo)
require(ggplot2)
require(Biostrings)



# all scores
library(vctrs)
allpaths = lapply(g_full,\(g){
  paths = longpaths(g)
  paths = unlist(paths,recursive = F)
  paths
})

pssms = unlist(leid_pssms, recursive = F)
scores_all = future_lapply(1:4,\(ab_i){
  paths = allpaths[[ab_i]]
  sapply(paths,\(path) {
    path = unlist(strsplit(path,""))
    
    sapply(pssms,\(mot){
      score(mot,path)
      
    })
    
  })
})

for(i in 1:4){
  colnames(scores_all[[i]]) = allpaths[[i]]
}
alln = unlist(sapply(leid_pssms,names))
for(i in 1:4){
  rownames(scores_all[[i]]) = alln
}

max_scores_all = future_sapply(scores_all,\(mat){
  apply(mat,1,max)
})

hist(max_scores_all)
save(max_scores_all, file = "similarity_scores_leiden.RData")

# plot score distributions
library(vioplot)
colors = c("red","green","blue","orange")

codes = unlist(sapply(1:4,\(i) rep(i,length(leid_pssms[[i]]))))
for(i in 1:4){
  
  
  sc = max_scores_all[which(codes==i),i]
  control = max_scores_all[,i]
  
  hist(control,  col=adjustcolor( "grey", 0.5), main = abs[i], breaks = 20)
  
  
  hist(sc,add = T, col = adjustcolor(colors[i],0.6))
  
  
}
# rnd scores on big computer ----
# scores with pssm from background frequencies ----
load("rnd_pepsets_pssms_libfreqs_leiden.RData")
plan(multisession(workers = 20))
rnd_scores_bkg = future_sapply(allpaths,\(paths){
  sapply(rnd_pssm,\(mot){
    scores = unlist(sapply(paths,\(path) {
      path = unlist(strsplit(path,""))
      score(mot,path)
      
    }))
    
    max(scores)
    
  })
})
# save(lib_pepsets,rnd_pssm,rnd_scores_bkg, file = "rnd_leiden_motif_scores_bkg.RData")  # no weights, with peptides from library

save(rnd_pssm,rnd_scores_bkg, file = "rnd_leiden_motif_scores_bkg_w.RData") # weights, with probabilities

hist(sample(rnd_scores_bkg,12000))
hist(max_scores_all, col = adjustcolor("green", alpha = 0.5),add = T)



#plot with random scores as controls----

#plot with random scores as controls----
library(vioplot)
colors = c("red","green","blue","orange")
load(file = "rnd_leiden_motif_scores_bkg_w.RData")
sizes = lapply(leiden_lcssets,\(sets) sapply(sets,length))
# plot on one plot
par(mfrow = c(1,1))
q1 = 0
q2 = 0.5
par(xpd = NA, bty = "n")
plot(1, type = "n",bty = "n", xlab = "", ylab = "similarity score",ylim = range(max_scores_all)+c(0,1),xlim = c(0,9.5),axes = F, frame.plot = FALSE)
#axis(1, labels = FALSE)  # x-axis without labels
axis(2, labels = T) 
axis(1,at = seq(1,4)*2 - 0.5,labels = abs,tick = FALSE, lwd = 0)
for(i in 1:4){
  
  sc = max_scores_all[which(codes==i),i]
  #sc = sc[which(sizes[[i]]>quantile(sizes[[i]],q1) & sizes[[i]]<quantile(sizes[[i]],q2))]
  control =rnd_scores_bkg[,i]
  
  pos = i*2
  
  colin = adjustcolor(colors[i],0.4)
  colbd = colors[i]
  vioplot(sc,add = T,at = pos,wex = 1, col=colin, border=colbd, rectCol=colbd, lineCol=colbd)
  
  
  colin = "lightgrey"
  colbd = "grey"
  at = i*2-1
  vioplot(control,add = T,at = at,wex = 1, col=colin, border=colbd, rectCol=colbd, lineCol=colbd,bty = "n")
  
  # p = p.adjust(pnorm(zsc[[i]],lower.tail = F),method = "BH")
  c = min(sc[which(p[[i]]<=0.05)])
  lines(c(pos-0.5,pos+0.5),c(c,c))
  
  c = min(sc[which(p[[i]]<=0.1)])
  lines(c(pos-0.5,pos+0.5),c(c,c),lty = 3)
  
  
  pp = (wilcox.test(sc, control,"greater")$p.value)
  if(pp<0.05){
    points(pos-0.5,18,pch = 8)
  }
}
box(bty = "n")
legend("bottomright",legend = c("0.05","0.1"),cex = 1,lty = c(1,3),title= "Corrected p \nfrom control \ndistribution ",bty = "n")
legend("topleft", pch = 8,legend = "p<0.05",title = "Mann-Whitney test",bty = "n")


# z-scores
codes = unlist(sapply(1:4,\(i) rep(i,length(leiden_pepsets_un[[i]]))))
p = lapply(1:4,\(i){
  sc = max_scores_all[which(codes==i),i]
  control =rnd_scores_bkg[,i]
  p = ecdf(control)(sc)
  p = p.adjust((1-p),method = "BH")
  names(p) = names(leiden_pepsets_un[[i]])
  p
})

sapply(p, min)
save(p,file = "sim_score_pvals.RData")


bestmots = sapply(p,\(ps) which(ps==min(ps)))
bestpl = c("WHPKDPD","RQIYSLE","HGCQMVS","NREPVID") # Herceptin_292 21c_279 17b_356 b12_808 , chosen from the clusters with best scores
# VCCTHSS from b12_356 is also good
bestpdb = c("PWGPDAA","RQIYSLE","MVPRSVC","DDERMRC") # from dbscan
which_cl = sapply(1:4,\(i) sapply(leiden_pepsets_un[[i]],\(set) length(intersect(bestpdb,set))))
which_cl = sapply(which_cl,\(s) names(which(s>0)))

worstmots = sapply(p,\(ps) which(ps>0.95*max(ps)))
worstpeps = c("SINHGLM","QSHPDLF","SSSNYFA","LFLRTYL") # Herceptin_354 21c_126 17b_811 b12_970
