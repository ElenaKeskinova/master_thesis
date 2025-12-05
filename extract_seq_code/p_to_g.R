library(igraph)


f = "D:/Elena/ban/Motifier_DataSet/Exp12.f_UNI-6/b12/uniqueP_AGGTC.txt"
name = basename(f)
name = substr(name,9,nchar(name)-4)
dir = dirname(f)

AllPep = read.table(file = f, colClasses = c("character","integer"),col.names = c("Pep","Freq"))
hist(AllPep$Freq)

F1 = file(paste0(substr(f,1,nchar(f)-4),"_p.txt"),"w")
cat(AllPep$Pep,file = F1, sep = "\n")
closeAllConnections()

# clean manually with TUPSCAN and create a "_clean" file 
####

Pepclean = sort(read.table(file = paste0(dir,"/",name,"_clean.txt"),
                           colClasses = "character",col.names = "Pep")$Pep)
Freqclean = sapply(Pepclean, function(x) AllPep$Freq[which(AllPep$Pep == x)])


g = newgraph(Pepclean)
g = set_vertex_attr(g,"Freq",value = Freqclean)
g_d = subgraph(g,which(degree(g)>0))
#distrcheck(g_d)
write.graph(g,paste0(dir,"/",name,".graphml"),"graphml")
g = read.graph(paste0(dir,"/",name,".graphml"),format = "graphml")

g_1c = subgraph(g,which(components(g)$membership==which(sizes(components(g)) >= vcount(g)/2)))

cluster = cluster_leiden(g_1c, objective_function = "modularity",resolution_parameter = 1)
g_1c = set_vertex_attr(g_1c,"Cl1",value = cluster$membership)
write.graph(g_1c,paste0(dir,"/",name,"_1c.graphml"),"graphml")

library(entropy)
entr = sapply(V(g_1c)$name, function(x) { entropy(table(strsplit(x,"")), method = "ML")})
boxplot(log10(degree(g_1c))~entr)


cl_centr = sapply(1:cluster$nb_clusters, function(x) which(eigen_centrality(subgraph(g_1c, which(cluster$membership == x)))$vector == 1))
cl_cen = unlist(cl_centr)
entr_cl_centr1 = sapply(names(cl_cen), function(x) { entropy(table(strsplit(x,"")), method = "ML")})
mean(entr_cl_centr1)
mean(entr)


