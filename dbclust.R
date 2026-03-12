
chooseep_db = function(coord,knnd,nn){
  require(dbscan)
  q1 = 0.2*length(knnd)
  q2 = 0.95*length(knnd)
  
  dj=knnd[seq(q1,q2,1)]
  rdj=unique(round(dj,4))
  proct=proc.time()
  # ncores=detectCores()-2
  # plan(multisession,workers=ncores)
  len = sapply(rdj,\(di){
    cldbsc=dbscan(coord,eps=di,minPts=nn)$cluster
    length(unique(cldbsc))
  })
  ij=which.max(len)
  
  #print(proc.time()-proct)
  epi=rdj[ij]
  #plot(rdj,len)
  #plot(knnd)
  #points(which.min(abs(knnd-epi)),epi,col=2, pch=16)
  
  cldbsc=dbscan(coord,eps=epi,minPts=nn)$cluster
  ncl=length(unique(cldbsc))
  # print(ncl)
  # print(unique(cldbsc))
  return(list(cldbsc,ncl,epi,len,rdj))
}

#cluster
dbclust_peps = function(lg)
{
  require(vctrs)
  require(dbscan)
  edg_lcs = lg$edg_lcs
  #spectral embedding
  arpopt=list(maxiter=100000, tol=1e-6)
  
  ic = which.max(igraph::components(lg)$csize)
  lgg = subgraph(lg,V(lg)[which(igraph::components(lg)$membership == ic)])
  L=embed_laplacian_matrix(lgg, no=35, which="sa", type="I-DAD", options=arpopt)#
  
  print(paste("laplacian ready"))
  
  opdim=dim_select(L$D[3:35])+2 # we want at least two dimensions (2nd and third)
  coord = L$X[,2:opdim]
  #proj on sphere
  
  projS=coord/(sqrt(rowSums(coord^2)))
  #dbscan
  opdim = opdim-1 # first eigenvector is not used
  nn = opdim*2 
  knnd=kNNdist(projS,nn)
  knnd=sort(knnd)
  result_db = chooseep_db(projS,knnd,nn)
  cldbsc = result_db[[1]]
  ncl = result_db[[2]]
  
  cls = sort(unique(cldbsc))
  lcssets = lapply(cls[which(cls!=0)], function(i) V(lgg)[which(cldbsc==i)]$name)
  lcsgsets = lapply(lcssets, \(set) {
    sg = subgraph(lgg, set)
    c = which(igraph::components(sg)$csize >=nn)
    lcscon = lapply(c, function(ci) V(sg)[which(igraph::components(sg)$membership == ci)]$name)
    lcscon
  })
  lcsgsets = list_flatten(lcsgsets)
  # add the other components
  ic = which(igraph::components(lg)$csize >=nn & igraph::components(lg)$csize<max(igraph::components(lg)$csize))
  lcs2 = lapply(ic, function(c) V(lg)[which(igraph::components(lg)$membership == c)]$name)
  if(length(lcs2)==1){
    lcsgsets = append(lcs2,lcsgsets)
  } else { lcsgsets = c(lcsgsets,lcs2)}
  
  
  lcssets5 = lapply(lcsgsets,\(set){
    c(unlist(sapply(set,\(pep){
      if(nchar(pep)==6){
        sapply(1:6,\(i) paste0(substr(pep,1,(i-1)),substr(pep,(i+1),6)))
      }
      else{
        pep
      }
    })))
  })
  
  pepsets = lcs_to_peps(lcssets5,lg$pep_lcs)
  
  
  pepsets = list_drop_empty(pepsets)
  
  freqs = sapply(pepsets,\(set) lg$pep_reps[set] )
  newsets = lapply(1:length(freqs),\(ii){
    rep(pepsets[[ii]],freqs[[ii]])
  })
  
  return(newsets)
  
}



gen_lcssets = function(lg,result_db,opdim){
  require(vctrs)
  require(igraph)
  require(purrr)
  cldbsc = result_db[[1]]
  ncl = result_db[[2]]
  
  cls = sort(unique(cldbsc))
  i = which.max(igraph::components(lg)$csize)
  lgg = subgraph(lg,V(lg)[igraph::components(lg)$membership == i])
  
  lcssets = lapply(cls[which(cls!=0)], function(i) V(lgg)[which(cldbsc==i)]$name)
  lcsgsets = lapply(lcssets, \(set) {
    sg = subgraph(lg, set)
    c = which(igraph::components(sg)$csize >=opdim*2)
    lcscon = lapply(c, function(ci) V(sg)[which(igraph::components(sg)$membership == ci)]$name)
    lcscon
  })
  lcsgsets = list_flatten(lcsgsets)
  
  
  ic = which(igraph::components(lg)$csize >=opdim*2 & igraph::components(lg)$csize<max(igraph::components(lg)$csize))
  lcs2 = lapply(ic, function(c) V(lg)[which(igraph::components(lg)$membership == c)]$name)
  if(length(lcs2)==1){
    lcsgsets = append(lcs2,lcsgsets)
  } else { lcsgsets = c(lcsgsets,lcs2)}
  lcsgsets
}

lcs_to_peps = function(lcssets, pep_lcs){
  
  lapply(lcssets,\(set) {
    numlcs = sapply(pep_lcs$lcs5,\(l) length(intersect(set,l)))
    peps = pep_lcs$p[numlcs>=2]
    peps
  })
}

