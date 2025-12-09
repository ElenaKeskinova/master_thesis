
# functions to compute paths and scores ----

longpaths = function(g){ # graph
  vs = V(g)[which(degree(g)==1)] # take paths only between terminal nodes
  
  peps = sapply(vs,\(v){
    paths = all_simple_paths(g,from =v,to = vs)
    peps = (sapply(paths,\(path){
      p = V(g)[names(path)]$AA
      np = length(p)
      nx = length(which(p=="-"))
      if(np-nx>4 & nx<np/3 ){ 
        paste0(p,collapse = "")
      }
      # return an array with string paths for each antibody
      
    }))
  })
  unlist(peps)
}

paths7 = function(g){ # all paths of length 7
  
  peps = sapply(V(g),\(v){
    paths = all_simple_paths(g,from = v, cutoff = 7)
    peps = (sapply(paths,\(path){
      p = V(g)[names(path)]$AA
      np = length(p)
      nx = length(which(p=="-"))
      if(np==7 & nx<np/3 ){ 
        paste0(p,collapse = "")
      }
      # return an array with string paths for each antibody
      
    }))
  })
  unlist(peps)
  
}

paths7_sc = function(g){ # sum of contacts of all paths of length 7, graph with node attribute "contacts" is needed
  
  sc = sapply(V(g),\(v){
    paths = all_simple_paths(g,from = v, cutoff = 7)
    scv = (sapply(paths,\(path){
      p = V(g)[names(path)]$AA
      np = length(p)
      nx = length(which(p=="-"))
      if(np==7 & nx<np/3 ){ 
        scp = sum(V(g)[names(path)]$contacts) # sum of contacts of path
        names(scp) = paste0(p,collapse = "")
        scp
      }
      # return an array with path contacts for each antibody
      
    }))
    scv
  })
  unlist(sc)
  
}

score = function(mot,epi){ # provide a pssm and an epitope as a character array, sum of positive scores
  n = length(epi)
  if(n<=ncol(mot)){
    max(sapply(1:(ncol(mot)-n+1),\(i){
      d = (diag(mot[epi,i:(i+n-1)]))
      sum(d[which(d>0)])
    }))
  } 
  else{
    max(sapply(1:(n - ncol(mot)+1),\(i){
      d = (diag(mot[epi[i:(i+ncol(mot)-1)],]))
      sum(d[which(d>0)])
    }))
  }
  
}

score_prob = function(mot,epi){ # provide a pssm and an epitope as a character array, sum of all positions
  n = length(epi)
  if(n<=ncol(mot)){
    max(sapply(1:(ncol(mot)-n+1),\(i){
      d = (diag(mot[epi,i:(i+n-1)]))
      sum(d)
    }))
  } 
  else{
    max(sapply(1:(n - ncol(mot)+1),\(i){
      d = (diag(mot[epi[i:(i+ncol(mot)-1)],]))
      sum(d)
    }))
  }
  
}

score_binary = function(mot,epi){ # provide a pssm and an epitope as a character array, number of scires>1
  n = length(epi)
  if(n<=ncol(mot)){
    max(sapply(1:(ncol(mot)-n+1),\(i){
      d = (diag(mot[epi,i:(i+n-1)]))
      sum(d>=1)
    }))
  } 
  else{
    max(sapply(1:(n - ncol(mot)+1),\(i){
      d = (diag(mot[epi[i:(i+ncol(mot)-1)],]))
      sum(d>=1)
    }))
  }
  
}

score_blosum = function(mot,path,blosum){ # provide a ppm and an epitope as a character array
  n = length(path)
  if(n<=ncol(mot)){
    max(sapply(1:(ncol(mot)-n+1),\(i){
      
      m = mot[,i:(i+n-1)]*cbind(blosum[,path]) # part of motif times substitution frequency of path aas
      mean(colSums(m)[which(path!="-")])
    }))
  } 
  else{
    max(sapply(1:(n - ncol(mot)+1),\(i){
      usedpath = path[i:(i+ncol(mot)-1)]
      m = mot*cbind(blosum[,usedpath])
      mean(colSums(m)[which(usedpath!="-")])
    }))
  }
}


score_pep = function(pep,path,similarity_m){ # provide a mimotope , epitope path as a string, and matrix of amino acid similarities
  path = unlist(strsplit(path,""))
  pep = unlist(strsplit(pep,""))
  n = length(path)
  np = length(pep)
  m = similarity_m[pep,path]
  m[,which(path == "-")] = NA
  if(np<n){
    d = split(m, row(m) - col(m))[as.character(c((np-n):0))]
  }
  else{
    d = split(m, row(m) - col(m))[as.character(c(0:(np-n)))]
  }
  max(sapply(d,mean,na.rm = T))
}

