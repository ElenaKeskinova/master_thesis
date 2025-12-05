specclust <- function(ab, lg,bc) {
  # name of antibody, line graph to use, barcode for saving results
  require(igraph)
  require(parallel)
  require(pbapply)
  require(future.apply)
  require(qualV)
  require(cluster)
  require(clusterCrit)
  require(scales)
  require(pROC)
  require(rgl)
  require(dbscan)
  require(RandPro)
  
  path = paste0("mixed-7graphs/",ab,"/")
  
  i = which.max(components(lg)$csize)
  g = subgraph(lg,V(lg)[components(lg)$membership == i])
  
  
  
  ## Make coordinates
  
  arpopt=list(maxiter=100000, tol=1e-6)
  
  ncores=detectCores()-2
  cpl=colorRampPalette(c("#000000FF","#0050AA9F","#10AA109F","#50AF3055","#FFFF009F","#FFA0009F","#B50000"), alpha=T)
  
  L=embed_laplacian_matrix(g, no=35, which="sa", type="I-DAD", options=arpopt)#
  opdim=dim_select(L$D)
  coord = L$X[,2:opdim]
  #proj on sphere
  projS=coord/(sqrt(rowSums(coord^2)))
  #dbscan
  nn = opdim*2
  knnd=kNNdist(projS,nn)
  knnd=sort(knnd)
  
  result_db = chooseep_db(projS,knnd,nn)
  save(result_db,file = paste0(path,ab,"_",bc,"_","dbscan.RData"))
  save(coord,file = paste0(path,ab,"_",bc,"_","speccoord.RData") )
  cldbsc = result_db[[1]]
  ncl = result_db[[2]]
  hist(cldbsc, main = ab)
  
  # rndmx=form_matrix(ncol(coord),3, JLT=F)
  # mxsm=coord%*%rndmx
  # mxsm=mxsm/(sqrt(rowSums(mxsm^2)))  
  # plot3d(mxsm[cldbsc>0,],col=cpl(ncl-1)[cldbsc[cldbsc>0]])  
  # 
  # # kmeans
  # 
  # result_km = choose_kmeans(projS,opdim)
  # save(result_km,file = paste0(path,ab,"_kmeans.RData"))
  # clkmn = result_km[[1]]
  # hist(clkmn)
  # plot3d(mxsm,col=cpl(ncl-1)[clkmn])  
  # 
  # tcl=table(cldbsc,clkmn)
  # tcl
}

chooseep_db = function(coord,knnd,nn){
  q1 = 0.2*length(knnd)
  q2 = 0.95*length(knnd)
  
  dj=knnd[seq(q1,q2,1)]
  rdj=unique(round(dj,4))
  proct=proc.time()
  ncores = detectCores()-1
  plan(multisession,workers=ncores)
  len = future_sapply(rdj,\(di){
    cldbsc=dbscan(coord,eps=di,minPts=nn)$cluster
    length(unique(cldbsc))
  }, future.seed=T)
  ij=which.max(len)
  plan(sequential)
  print(proc.time()-proct)
  epi=rdj[ij]
  plot(rdj,len)
  plot(knnd)
  points(which.min(abs(knnd-epi)),epi,col=2, pch=16)
  
  cldbsc=dbscan(coord,eps=epi,minPts=nn)$cluster
  ncl=length(unique(cldbsc))
  print(ncl)
  return(list(cldbsc,ncl,epi,len,rdj))
}


choose_kmeans = function(coord,opdim){
  rng=round(opdim/2):round(2*opdim)
  projMx=round(coord,4)
  proct=proc.time()
  plan(multisession, workers=ncores)
  kmN=rng[which.max(rowMeans(future_sapply(1:ncores,\(i){
    x=t(sapply(rng, \(n){
      cl=kmeans(projMx,n,iter.max = 1000, nstart = 10)
      y=as.numeric(unlist(intCriteria(projMx,cl$cluster,c("Dunn","Calinski_Harabasz","PBM","Silhouette","Point_biserial"))))
      if (any(is.na(y))) y[is.na(y)]=mean(y[!is.na(y)])
      return(y)
    }))
    #print(sum(is.na(x)))
    sqrt(rowSums(apply(x,2,\(co) (co-min(co))/diff(range(co)))^2))
  }, future.seed=T)))]
  plan(sequential)
  print(proc.time()-proct)
  
  clkmn=kmeans(projMx,kmN,iter.max = 1000, nstart = 50)$cluster
  return(list(clkmn, kmN))
}
