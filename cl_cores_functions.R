# functions for cluster filtering that work with membership vectors


cl_cores = function(g, memb, size_cutoff,deg_cutoff){ # find the cores of clusters and discard small cores
  # graph, membership vector, minimum cluster size
  require(igraph)
  require(vctrs)
  cores = lapply(unique(memb),\(c){
    subg = subgraph(g,which(memb==c))
    deg_ratio =  igraph::degree(subg)/igraph::degree(g, which(memb==c))
    
    core = which(deg_ratio>deg_cutoff)
    if(length(core)>=size_cutoff){
      return(unlist(V(subg)[core]$name))
    }
  })
  
   
    
  
  cores = list_drop_empty(cores)
  return(cores) # pepsets of the cluster cores
 
}

cl_cores_v = function(g, memb, size_cutoff,deg_cutoff){ # find the cores of clusters and discard small cores
  # graph, communit object, minimum cluster size
  # return a membership vector
  require(igraph)
  membc = memb
  for(c in unique(memb)){
    cli = which(memb==c)
    #if(length(cli)>2){
      subg = subgraph(g,cli)
      deg_ratio =  igraph::degree(subg)/igraph::degree(g, cli)
      
      core = which(deg_ratio>deg_cutoff)
      membc[cli] = 0
      if(length(core)>size_cutoff){
        membc[cli[core]] = c
      }
    #}
    
  }
  

  
  return(membc) # membership vector with zeros where in periphery
  
}

cl_indegree = function(g, memb){ # return degree ratio for every vertex for its assugned cluster
  # graph, communit object
  # return a membership vector
  require(igraph)
  degr = array(0, dim = length(memb))
  for(c in unique(memb)){
    cli = which(memb==c)
    if(length(cli)>2){
      subg = subgraph(g,cli)
      deg_ratio =  igraph::degree(subg)/igraph::degree(g, cli)
      degr[cli] = deg_ratio
    
    }
    
  }
  return(degr)
}

cut_small_clusters= function(memb, mins){
  # membership vector and minimal size
  require(igraph)
  memb2 = array(0, dim = length(memb))
  for(c in unique(memb)){
    cli = which(memb==c)
    if(length(cli)>=mins){
      memb2[cli] = c
    }
  }
  return(memb2)
}

clust_big = function(vgn, memb, mins){ # names of graph vertices, membership vector, minimal size
  mmbig = which(table(memb)>=mins)
 
  lapply(mmbig,\(mi){
    vgn[which(memb==mi)]
  })
}

lcs_to_peps = function(lcssets, pep_lcs){
  
  lapply(lcssets,\(set) {
    numlcs = sapply(pep_lcs$lcs5,\(l) length(intersect(set,l)))
    peps = pep_lcs$p[numlcs>=2]
    peps
  })
}