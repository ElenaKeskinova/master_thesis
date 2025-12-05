# create graph from an array of sequences and a threshhold of mismatches to form an edge


newgraph = function(seqdata, mistakes){
  require(future.apply)
  require(igraph)
  require(comparator)
  m=pairwise(comparator::LCS(),seqdata,return_matrix = TRUE)
  m[]=m<=(mistakes*2)&m>0
  newgraph=graph_from_adjacency_matrix(m,"undirected")
  newgraph = set_vertex_attr(newgraph, "name", value=seqdata)
  return(newgraph)
}

