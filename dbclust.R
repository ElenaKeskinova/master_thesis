
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
dbclust_peps = function(g,lg)
{
  require(vctrs)
  edg_lcs = lg$edg_lcs
  #spectral embedding
  arpopt=list(maxiter=100000, tol=1e-6)
  
  i = which.max(components(lg)$csize)
  lg = subgraph(lg,V(lg)[which(components(lg)$membership == i)])
  L=embed_laplacian_matrix(lg, no=35, which="sa", type="I-DAD", options=arpopt)#
  
  print(paste(i,"laplacian ready"))
  
  opdim=dim_select(L$D)
  coord = L$X[,2:opdim]
  #proj on sphere
  projS=coord/(sqrt(rowSums(coord^2)))
  #dbscan
  nn = opdim*2
  knnd=kNNdist(projS,nn)
  knnd=sort(knnd)
  result_db = chooseep_db(projS,knnd,nn)
  cldbsc = result_db[[1]]
  ncl = result_db[[2]]
  # 
  # cls = sort(unique(cldbsc))
  # lcssets = lapply(cls[which(cls!=0)], function(i) V(lg)[which(cldbsc==i)]$name)
  # pepsets = lapply(lcssets,\(set) {
  #   peps = unlist(sapply(set,\(lcs) c(edg_lcs[which(edg_lcs[,3]==lcs),1:2])))
  #   unique(c(peps))
  # })
  # ic = which(igraph::components(g)$csize >=10 & igraph::components(g)$csize<max(igraph::components(g)$csize))#
  # peps2 = lapply(ic, function(c) V(g)[which(igraph::components(g)$membership == c)]$name)
  # if(length(peps2)==1){
  #   pepsets = append(peps2,pepsets)
  # } else { pepsets = c(pepsets,peps2)}
  # 
  # pepsets = list_drop_empty(pepsets)
  # return(pepsets)
  
  gen_pepsets(lg,g,edg_lcs,result_db)
  
}

dbclust_peps_cl = function(g,lg) # return pepsets and clustering results
{
  require(vctrs)
  require(dbscan)
  require(igraph)
  edg_lcs = lg$edg_lcs
  #spectral embedding
  arpopt=list(maxiter=100000, tol=1e-6)
  
  i = which.max(components(lg)$csize)
  lg = subgraph(lg,V(lg)[which(components(lg)$membership == i)])
  L=embed_laplacian_matrix(lg, no=35, which="sa", type="I-DAD", options=arpopt)#
  
  print(paste(i,"laplacian ready"))
  
  opdim=dim_select(L$D)
  coord = L$X[,2:opdim]
  #proj on sphere
  projS=coord/(sqrt(rowSums(coord^2)))
  #dbscan
  nn = opdim*2
  knnd=kNNdist(projS,nn)
  knnd=sort(knnd)
  result_db = chooseep_db(projS,knnd,nn)
  cldbsc = result_db[[1]]
  ncl = result_db[[2]]
  
  pp = gen_pepsets(lg,g,edg_lcs,result_db)
  return(list(pp,result_db))
  
}

gen_pepsets = function(lg,g,edg_lcs,result_db){
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
  pepsets = lapply(lcsgsets,\(set) {
    peps = unlist(sapply(set,\(lcs) c(edg_lcs[which(edg_lcs[,3]==lcs),1:2])))
    unique(c(peps))
  })
  ic = which(igraph::components(g)$csize >=10 & igraph::components(g)$csize<max(igraph::components(g)$csize))#
  peps2 = lapply(ic, function(c) V(g)[which(igraph::components(g)$membership == c)]$name)
  if(length(peps2)==1){
    pepsets = append(peps2,pepsets)
  } else { pepsets = c(pepsets,peps2)}
  
  pepsets = list_drop_empty(pepsets)
  return(pepsets)
  
}

gen_lcssets = function(lg,result_db){
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
    c = which(igraph::components(sg)$csize >=5)
    lcscon = lapply(c, function(ci) V(sg)[which(igraph::components(sg)$membership == ci)]$name)
    lcscon
  })
  lcsgsets = list_flatten(lcsgsets)
  
  
  ic = which(igraph::components(lg)$csize >=5 & igraph::components(lg)$csize<max(igraph::components(lg)$csize))#
  lcs2 = lapply(ic, function(c) V(lg)[which(igraph::components(lg)$membership == c)]$name)
  if(length(lcs2)==1){
    lcsgsets = append(lcs2,lcsgsets)
  } else { lcsgsets = c(lcsgsets,lcs2)}
  
  lcsgsets = list_drop_empty(lcsgsets)
  return(lcsgsets)
  
}

gen_pepsets_old = function(lg,g,edg_lcs,result_db){
  require(vctrs)
  require(igraph)
  require(purrr)
  cldbsc = result_db[[1]]
  ncl = result_db[[2]]
  
  cls = sort(unique(cldbsc))
  i = which.max(components(lg)$csize)
  lgg = subgraph(lg,V(lg)[components(lg)$membership == i])
  # 
  lcssets = lapply(cls[which(cls!=0)], function(i) V(lgg)[which(cldbsc==i)]$name)
  # lcsgsets = lapply(lcssets, \(set) {
  #   sg = subgraph(lg, set)
  #   c = which(igraph::components(sg)$csize >=10)
  #   lcscon = lapply(c, function(ci) V(sg)[which(igraph::components(sg)$membership == ci)]$name)
  #   lcscon
  # })
  # lcsgsets = list_flatten(lcsgsets)
  pepsets = lapply(lcssets,\(set) {
    peps = unlist(sapply(set,\(lcs) c(edg_lcs[which(edg_lcs[,3]==lcs),1:2])))
    unique(c(peps))
  })
  ic = which(igraph::components(g)$csize >=10 & igraph::components(g)$csize<max(igraph::components(g)$csize))#
  peps2 = lapply(ic, function(c) V(g)[which(igraph::components(g)$membership == c)]$name)
  if(length(peps2)==1){
    pepsets = append(peps2,pepsets)
  } else { pepsets = c(pepsets,peps2)}
  
  pepsets = list_drop_empty(pepsets)
  return(pepsets)
  
}