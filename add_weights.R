# put weights on the edges of line graphs, according to the copy number of the respective sequences
add_weights = function(ab,logw = F){
  require(future.apply)
  path = paste0("mixed-7graphs/",ab,"/")
  
  load(paste0(path,ab,"big7or.RData"))
  load(paste0(path,ab,"big7-l.RData"))
  load(file = paste0(path,ab,"_all_lcs.RData"))
  
  v_freqs = V(G)$Freq
  names(v_freqs) = V(G)$name
  
  if(logw == T){
    v_freqs = log2(v_freqs)
  }
  
  lcs = V(lg)$name
  plan(multisession)
  weights = weights(lg,edg_lcs,v_freqs)
  set_edge_attr(lg,name = "weight",value = weights)
}

add_weights_corr = function(ab,bg = bgmot){ # correct for closeness to background
  require(future.apply)
  require(scales)
  path = paste0("mixed-7graphs/",ab,"/")
  
  load(paste0(path,ab,"big7or.RData"))
  load(paste0(path,ab,"big7-l.RData"))
  load(file = paste0(path,ab,"_all_lcs.RData"))
  
  v_freqs = V(G)$Freq
  names(v_freqs) = V(G)$name
  
  oddsP = sapply(V(G)$name,\(p){
    p = unlist(strsplit(p,""))
    prob = prod(sapply(1:7,\(i) bgmot[p[i],i]))
    prob/(0.05**7) # how probability is increased compared to totally rnd peptide
  })
  v_freqs = v_freqs/oddsP
  m = 1-min(v_freqs)
  v_freqs = log2(v_freqs+m)
  
  
  lcs = V(lg)$name
  plan(multisession)
  weights = weights(lg,edg_lcs,v_freqs)
  set_edge_attr(lg,name = "weight",value = weights)
}

weights = function(lg,edg_lcs,v_freqs){ # compute edge weights from copy number of sequenced phages
  
  future_sapply(E(lg),\(e){
    ee = ends(lg,e)
    ee1 = ee[1]
    ee2 = ee[2]
    
    seqs1 = unlist(edg_lcs[which(edg_lcs[,3]==ee1),1:2])
    seqs2 = unlist(edg_lcs[which(edg_lcs[,3]==ee2),1:2])
    common = intersect(seqs1,seqs2)
    sum(v_freqs[common])
  })
}