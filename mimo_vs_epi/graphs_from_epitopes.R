library(Biostrings)
library(igraph)

# take files with results from Pymol for contact amino acids on antigen surfaces, files are made with the function analyze_contacts2 from extract_epi_dist.py 
# obtained from structures 1yyl,3lqa,2ny7,1n8z from pdb structure, structures are cleaned from solvent, NAG, SO4 and alternative structures before extracting distances
# with contact cutoff 5A between epitope and paratope residues
# make graphs from epitope amino acids, introduce gap nodes where distance is over 6 A 

AA_code = sapply(AMINO_ACID_CODE,toupper)
AA_code = c(AA_code,"XXX")
names(AA_code)[27] = "-"
abs = c("Herceptin","21c","17b","b12")

distances = lapply(abs,\(ab){ # distances between amino acids closer than 5 A to antibody, from PDB, calculated with Pymol
  df = read.table(paste0(ab,"_epi_dist.txt"),header = T, sep = ",")
  df$aan1 = sapply(df$Residue1, substr, 1,3)
  df$aan2 = sapply(df$Residue2, substr, 1,3)
  df = df[which(df$aan1 %in% AA_code & df$aan2 %in% AA_code),] # only valid amino acids
  
})



# create distance matrix

dist_mx = lapply(distances,\(df){
  alaa = unique(unlist(df[,1:2]))
  m = matrix(0, nrow=length(alaa), ncol=length(alaa),
             dimnames=list(alaa, alaa))
  for(i in 1:nrow(df)){
    e1 = df[i,1]
    e2 = df[i,2]
    d = df[i,3]
    m[e1,e2] = d
    m[e2,e1] = d
  }
  m
})

# make graphs ----
### minimal graphs ----
# start with small graphs and the add edges to connect all amino acids of the epitopes

d = 6
ds = lapply(distances,\(df){
  df$Distance<=d
  
})

graphs_min = lapply(1:4,\(i){
  require(igraph)
  e = ds[[i]]
  edges = cbind(distances[[i]]$Residue1,distances[[i]]$Residue2)[e,]
  g = graph_from_edgelist(edges, directed = F)
  aas = sapply(V(g)$name,\(nm){
    c3 = substr(nm,1,3) 
    names(AA_code[which(AA_code==c3)])
    
  } )
  set_vertex_attr(g,name = "AA", value = unlist(aas))
  
})
names(graphs_min) = abs
for(i in 1:4){plot(graphs_min[[i]],main = abs[i])}




### select edges so that graphs are connected with minimal number of edges ----

disttemp = lapply(distances,\(df){
  df[df$Distance<=d,] 
  
})

distscan = lapply(distances,\(df){ # between 6 and 12
  df[df$Distance >d & df$Distance <=d*2 ,]
  
})
distscan = lapply(distscan,\(df){ # between 6 and 12
  o = order(df$Distance) # put in ascending order
  df[o,]
  
})


gtemp = graphs_min # graphs to add edges to



# iterate through nodes over 6 A to include the necessary ones for a connected graph
for(i in 1:4){
  df = distscan[[i]]
  for(j in 1:nrow(df)){
    vs = V(gtemp[[i]])$name
    row = df[j,]
    if((row$Residue1 %in% vs == F) | (row$Residue2 %in% vs == F) | components(gtemp[[i]])$membership[row$Residue1] != components(gtemp[[i]])$membership[row$Residue2] ){
      di = row$Distance
      # if the chosen edge connects 2 disconnected components, include it with a gap node in the middle
      print(di)
      if(di%/%d == 1){ # 1 gap node
        xnum = nrow(disttemp[[i]])
        n1 = paste0("XXX",xnum) # "gap" node 
        disttemp[[i]] = rbind(disttemp[[i]],c(row$Residue1,n1,di))
        disttemp[[i]] = rbind(disttemp[[i]],c(n1,row$Residue2,di))
        
      }else if (row$Distance%/%d == 2){ # 2 gap nodes
        xnum = nrow(disttemp[[i]])
        n1 = paste0("XXX",xnum)
        n2 = paste0("XXX",xnum+1)
        disttemp[[i]] = rbind(disttemp[[i]],c(row$Residue1,n1,di))
        disttemp[[i]] = rbind(disttemp[[i]],c(n1,n2,di))
        disttemp[[i]] = rbind(disttemp[[i]],c(n2,row$Residue2,di))
      }
      gtemp[[i]] = graph_from_data_frame(disttemp[[i]],directed = F)
    }
  }
  
}
for(i in 1:4){plot(gtemp[[i]],main = abs[i])}

g_full = gtemp

# make AA 1 letter code as attribute
g_full = lapply(g_full,\(g){
  aas = sapply(V(g)$name,\(nm){
    c3 = paste0(unlist(strsplit(nm,""))[1:3],collapse = ""); 
    names(AA_code[which(AA_code==c3)])
    
  } )
  set_vertex_attr(g,name = "AA", value = aas)
  
})
save(g_full,file = "epitope_graphs_6.RData") #with  distance 6
