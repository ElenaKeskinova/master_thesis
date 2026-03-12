rnd_graph = function(N,lett = 7,freqs=freq){
  
  a = sapply(1:N, function(x){
    paste(sample(AA_STANDARD, lett, prob = freqs[AA_STANDARD], replace = TRUE), collapse = "")
  })
  a = unique(a)
  g = newgraph(a,2)
  gl = graph_to_line(g)
  return(list(g,gl)) 
}

rnd_lgraph = function(N,lett = 7,freqs=freq){
  
  a = sapply(1:N, function(x){
    paste(sample(AA_STANDARD, lett, prob = freqs[AA_STANDARD], replace = TRUE), collapse = "")
  })
  a = unique(a)
  g = newgraph(a,2)
  gl = graph_to_line(g)
  return(gl) 
}

rnd_lgraph_repfreqs = function(N,lett = 7,freqs=freq, repfreqs){
  # provide a table with copy number frequencies
  a = sapply(1:N, function(x){
    paste(sample(AA_STANDARD, lett, prob = freqs[AA_STANDARD], replace = TRUE), collapse = "")
  })
  a = unique(a)
  reps = sample(as.numeric(names(repfreqs)),length(a),prob = repfreqs,replace = T)
  names(reps) = a
  g = newgraph(a,2)
  g = set_vertex_attr(g,name = "Freq",value = reps)
  gl = graph_to_line_plcs(g)
  w = weights(gl,gl$edg_lcs,reps)
  gl = set_edge_attr(gl,name = "weight",value = w)
  gl = set_graph_attr(gl, name =  "pep_reps", value = reps)
  return(gl) 
}

lgraph_repfreqs = function(a, repfreqs){
  # provide selected set of peptides a table with copy number frequencies
  
  a = unique(a)
  reps = sample(as.numeric(names(repfreqs)),length(a),prob = repfreqs,replace = T)
  names(reps) = a
  g = newgraph(a,2)
  g = set_vertex_attr(g,name = "Freq",value = reps)
  gl = graph_to_line_plcs(g)
  w = weights(gl,gl$edg_lcs,reps)
  gl = set_edge_attr(gl,name = "weight",value = w)
  return(gl) 
}





rnd_graph_w = function(N,lett = 7,freqs=freq,s){
  
  a = sapply(1:N, function(x){
    paste(sample(sort(AA_STANDARD), lett, prob = freqs, replace = TRUE), collapse = "")
  })
  fr = table(a)
  f = fr[which((fr)>1)]
  aa = names(fr[which((fr)>1)])
  g = newgraph(aa,2)
  g = set_vertex_attr(g,name="Freq",value = f)
  gl = graph_to_line(g)
  
  weights = weights(gl,gl$edg_lcs,f)
  gl = set_edge_attr(gl,name = "weight",value = weights)
  
  return(list(g,gl)) 
}

rnd_graph_wcopy = function(N,lett = 7,freqs=freq,pepfr,s){
  
  a = sapply(1:N, function(x){
    paste(sample(sort(AA_STANDARD), lett, prob = freqs, replace = TRUE), collapse = "")
  })
  
  g = newgraph(a,2)
  g = set_vertex_attr(g,name="Freq",value =pepfr)
  gl = graph_to_line(g)
  
  names(pepfr) = V(g)$name
  pepfr = log2(pepfr)
  weights = weights(gl,gl$edg_lcs,pepfr)
  gl = set_edge_attr(gl,name = "weight",value = weights)
  
  return(list(g,gl)) 
}

rnd_graph_reps = function(N,lett = 7,freqs=freq,repfreqs,s){ # needs frequencies of aa in alphabetical order(of 1 letter code), and frequencies of repetitions = log2(copy number of peptides)
  
  a = sapply(1:N, function(x){
    paste(sample(AA_STANDARD, lett, prob = freqs[AA_STANDARD], replace = TRUE), collapse = "")
  })
  a = unique(a)
  reps = sample(as.numeric(names(repfreqs)),length(a),prob = repfreqs,replace = T)
  names(reps) = a
  g = newgraph(a,2)
  g = set_vertex_attr(g,name = "Freq",value = reps)
  gl = graph_to_line(g)
  w = weights(gl,gl$edg_lcs,reps)
  gl = set_edge_attr(gl,name = "weight",value = w)
  
  return(list(g,gl)) 
}
