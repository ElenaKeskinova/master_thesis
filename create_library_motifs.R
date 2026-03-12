load("mixed-7graphs/lib_allpeps_7nC.RData")
source("newgraph.R")
source("graph_to_line.R")
source("rnd_graph.R")
source( "add_weights.R")
source("cl_cores_functions.R")

# test parameters ----

## rnd graphs with rnd peps ----
load(file = "lib_freqs.RData")
freq = rowSums(freqM7)/sum(freqM7)
repfreqs = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  round(log2(V(G)$Freq))
})
rep_freqs = table(unlist(repfreqs))/sum(table(unlist(repfreqs))) # frequencies of logarithmic copy numbers



lib_lg = rnd_lgraph_repfreqs(10000,freqs = freq, repfreqs = rep_freqs)
E(lib_lg)$weight = E(lib_lg)$weight/max(E(lib_lg)$weight)
ew = E(lib_lg)$weight

rl = sum(ew)/(vcount(lib_lg)*(vcount(lib_lg)-1)/2)
l_res = seq(rl,0.1,by=(0.1-rl)/20)

plan(multisession(workers = 20))
result_lib = sapply(1:5,\(i){
  lib_lg = rnd_lgraph_repfreqs(10000,freqs = freq, repfreqs = rep_freqs)
  E(lib_lg)$weight = E(lib_lg)$weight/max(E(lib_lg)$weight)
  print(sum(ew)/(vcount(lib_lg)*(vcount(lib_lg)-1)/2))
  future_sapply(l_res,\(res){
    n = 5
    cls = sapply(1:n,\(i) cluster_leiden(lib_lg,resolution = res, objective_function = "CPM",n_iterations = 10)$membership)
    m = 3
    cls2 = sapply(1:n,\(i) cut_small_clusters(cls[,i],m))
    ari = unlist(future_sapply(1:(n-1),\(i){sapply((i+1):n,\(j) {adjustedRandIndex(cls2[,i],cls2[,j]) })})) # mean ARI, better overlap if high
    clcount = mean(sapply(1:n,\(i)length(unique(cls2[,i]))))
    nontrash = mean(sapply(1:n,\(i) sum(cls2[,i]!=0)/nrow(cls2)))
    c(mean(ari),sd(ari),clcount,nontrash)
  },simplify =  "array")
  
},simplify =  "array")
  
## plot test results ----

plot(l_res,result_lib[1,,1],main = "mean ARI")
for(i in 2:5){
  points(l_res,result_lib[1,,i],col = i)
}
plot(l_res,result_lib[3,,1],main =  "clusters")
for(i in 2:5){
  points(l_res,result_lib[3,,i],col = i)
}

r_ari1 = cbind(l_res,result_lib[1,]) # best resolution is 0.07
r_ari = sapply(1:5,\(i){
  l_res[which(result_lib[1,,i]>=0.9*max(result_lib[1,,i]))[i]]
})
r = mean(r_ari) # 0.07

# with real peps ----
all_weights = unlist(sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-logw.RData"))
  E(lg)$weight
}))

hist(all_weights)
r = mean(bestreso_w)

peps = sample(peps_7nC,15000)
lib_g = newgraph(peps,2)
lib_lg = graph_to_line_plcs(lib_g)
rl = edge_density(lib_lg)
l_res = seq(rl,0.15,by=(0.15-rl)/20)


plan(multisession(workers = 20))
result_lib = future_sapply(l_res,\(res){
  n = 5
  cls = sapply(1:n,\(i) cluster_leiden(lib_lg,resolution = res, objective_function = "CPM",n_iterations = 10)$membership)
  m = 3
  cls2 = sapply(1:n,\(i) cut_small_clusters(cls[,i],m))
  ari = unlist(future_sapply(1:(n-1),\(i){sapply((i+1):n,\(j) {adjustedRandIndex(cls2[,i],cls2[,j]) })})) # mean ARI, better overlap if high
  clcount = mean(sapply(1:n,\(i)length(unique(cls2[,i]))))
  nontrash = mean(sapply(1:n,\(i) sum(cls2[,i]!=0)/nrow(cls2)))
  c(mean(ari),sd(ari),clcount,nontrash)
},simplify =  "array")

## plot test results ----

plot(l_res,result_lib[1,],main = "mean ARI")
plot(l_res,result_lib[3,],main =  "clusters")

res_ml1 = l_res[which.max(result_lib[1,])]
res_ml2 = l_res[which.max(result_lib[1,])]
res_ml3 = l_res[which.max(result_lib[1,])]
res_ml4 = l_res[which.max(result_lib[1,])]
r1 = l_res[which(result_lib[1,2:21]>0.99*max(result_lib[1,2:21]))[1]+1]
r2 = l_res[which(result_lib[1,2:21]>0.99*max(result_lib[1,2:21]))[1]+1]
# result: plateau starts at about 0.05

# library graph without weights with chosen resolution ----
r = 0.05
plan(multisession(workers = 20))
lib_lcssets = future_sapply(1:10,\(i){
  peps = sample(peps_7nC,15000)
  lib_g = newgraph(peps,2)
  lib_lg = graph_to_line_plcs(lib_g)
  
  ncl = 5
  cls = sapply(1:ncl,\(i) cluster_leiden(lib_lg,resolution = r, objective_function = "CPM",n_iterations = 10)$membership)
  m = 3
  cls2 = sapply(1:ncl,\(i) cut_small_clusters(cls[,i],m))
  ari = (future_sapply(1:(ncl),\(i){sapply((1):ncl,\(j) {adjustedRandIndex(cls2[,i],cls2[,j]) })})) # mean ARI, better overlap if high
  cl = cls[,which.max(colSums(ari))]
  m = 3
  lcssets = cl_cores(lib_lg,cl,m,0)
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
  
  pepsets = lcs_to_peps(lcssets,lib_lg$pep_lcs)
  save(pepsets,file = paste("rnd_pepsets",i,"g.RData"))
  lcssets
  
},future.seed = T)

lib_lcssets = unlist(lib_lcssets, recursive = F)
lib_motifs = lapply(lib_lcssets,\(set)  create_motif(set,alphabet = "AA",pseudocount = 1)@motif) #ppms
save(lib_lcssets,lib_motifs, file = "20k_lib_lcsmots.RData")

lib_pepsets = lapply(1:10,\(i) {load(file = paste("rnd_pepsets",i,"g.RData")); pepsets}) # to compare to epitopes
lib_pepsets = unlist(lib_pepsets,recursive = "F")
lib_pssm = future_lapply(lib_pepsets,\(set){
  mot = create_motif(set,type = "PWM",alphabet = "AA",bkg = rowSums(freqM7)/7, pseudocount = 1)@motif
  mot = rbind(mot,rep(0,ncol(mot)))
  rownames(mot)[nrow(mot)] = "-"
  mot
})

# save(lib_pepsets,rnd_pssm,rnd_scores_bkg, file = "rnd_leiden_motif_scores_bkg.RData")

