graph_to_line = function(g){
  # convert a graph to line graph
  source("compute_lcs.R")
  require(igraph)
  require(future.apply)
  plan(multisession(workers = 4))
  enames = future_sapply(E(g), function(x) compute_lcs(ends(g,x)[1],ends(g,x)[2]))
  all_lcs = cbind(ends(g,E(g)),enames)

  lg = make_line_graph(g)
  lg = set_vertex_attr(lg,name = "name",value = enames)
  lg = set_graph_attr(lg,all_lcs,name = "edg_lcs")
  lg_s = simplify(contract(lg, as.numeric(as.factor(enames)),vertex.attr.comb = list(name = "first") ))
  return(lg_s)
}
