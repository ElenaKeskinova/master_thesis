
move_renamed = function(src, dst){
  # move structures from top 2 clusters in different folder and rename them
  if(!dir.exists(dst)) dir.create(dst)
  prefix = unlist(strsplit(basename(src),"-"))[2]
  
  energy = read.table(paste0(src,"/clusters_haddock-sorted.stat_best4"))
  topmod = sapply(1:2, \(i) unlist(strsplit(energy[i,1],"clust"))[2])
  
  
  for(j in 1:length(topmod)){ # go through each of the top clusters
    if(!dir.exists(dst)) dir.create(dst)
    
    
    # List files with full paths
    cluster_num = topmod[j]
    files <- list.files(src,pattern = paste0("^cluster", cluster_num, "_[0-9]+\\.pdb$"), full.names = TRUE)
    
    # Extract just the file names
    old_names <- basename(files)
    modnums = substr
    
    # Add prefix
    new_names <- paste0(prefix, "_top", j,"_",seq(1,4), ".pdb")
    
    # Move & rename
    file.copy(files, file.path(dst, new_names))
  }
  
}


haddock_folders = list.files(path = "dock_motifier_mimotopes/haddock_results",full.names = T)
dst = "dock_motifier_mimotopes/best_haddock_models"

for(src in haddock_folders){
  move_renamed(src,dst)
}

# leiden best clusters for her, 17b and b12
haddock_folders = list.files(pattern = "*best_leiden",path = "dock_motifier_mimotopes/haddock_results",full.names = T)
dst = "dock_motifier_mimotopes/best_haddock_models_leiden"

for(src in haddock_folders){
  move_renamed(src,dst)
}

# leiden worst clusters , no file with cluster haddock scores
best2cl = list(c(1,2),c(7,1),c(3,2),c(1,3))

haddock_folders = list.files(pattern = "*worst_leiden",path = "dock_motifier_mimotopes/haddock_results",full.names = T)
dst = "dock_motifier_mimotopes/best_haddock_models_leiden"

for(src in haddock_folders){
  if(!dir.exists(dst)) dir.create(dst)
  prefix = unlist(strsplit(basename(src),"-"))[2]
  
  
  topmod = best2cl[[i]]
  
  
  for(j in 1:length(topmod)){ # go through each of the top clusters
    if(!dir.exists(dst)) dir.create(dst)
    
    
    # List files with full paths
    cluster_num = topmod[j]
    files <- list.files(src,pattern = paste0("^cluster", cluster_num, "_[0-9]+\\.pdb$"), full.names = TRUE)
    
    # Extract just the file names
    old_names <- basename(files)
    modnums = substr
    
    # Add prefix
    new_names <- paste0(prefix, "_top", j,"_",seq(1,4), ".pdb")
    
    # Move & rename
    file.copy(files, file.path(dst, new_names))
  }
}

