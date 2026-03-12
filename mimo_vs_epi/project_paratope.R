AA_code = sapply(AMINO_ACID_CODE,toupper)
AA_code = c(AA_code,"XXX")
names(AA_code)[27] = "-"
abs = c("Herceptin","21c","17b","b12")

epi_contacts = lapply(abs,\(ab){
  read.delim(file = paste0("dock_motifier_mimotopes/",ab,"_contacts_epi.txt"))
})
epi_contacts = lapply(epi_contacts,\(ab){
  ab$AB_contacts = lapply(ab$AB_contacts,\(c_list){aas = unlist(strsplit(c_list, split = ","));sapply(aas,\(a) substr(a,1,nchar(a)-1))})
  ab
})
names(epi_contacts) = abs

# contacts of best mimos
mimo_contacts = lapply(abs[1:4],\(ab){
  read.delim(file = paste0("dock_motifier_mimotopes/best_haddock_models/",ab,"_best_top1_1.pdb_contacts.txt"))
})
mimo_contacts = lapply(mimo_contacts,\(ab){
  ab$AB_contacts = lapply(ab$AB_contacts,\(c_list){aas = unlist(strsplit(c_list, split = ","));sapply(aas,\(a) substr(a,1,nchar(a)-1))})
  ab
})
names(mimo_contacts) = abs[1:4]

# data frame of distances ----
pdistances = lapply(abs,\(ab){ # distances between amino acids closer than 5 A to antibody, from PDB, calculated with Pymol
  df = read.table(paste0(ab,"_pardist.txt"),header = T, sep = ",")
  df$aan1 = sapply(df$Residue1, substr, 1,3)
  df$aan2 = sapply(df$Residue2, substr, 1,3)
  df = df[which(df$aan1 %in% AA_code & df$aan2 %in% AA_code),] # only valid amino acids
  df
})
names(pdistances) = abs
allaas = lapply(pdistances,\(df){
  alaa = unique(unlist(df[,1:2]))
  sapply(alaa,\(a) substr(a,1,nchar(a)-1))})


dist_mx = lapply(pdistances,\(df){
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

twodim_coord = lapply(dist_mx,\(dmat){
  c = cmdscale(dmat,k = 2)
  c
})


# plot paratope map----
for(i in 2){
  cnt = which(allaas[[i]] %in% unlist(epi_contacts[[i]]$AB_contacts))
  cntm = which(allaas[[i]] %in% unlist(mimo_contacts[[i]]$AB_contacts))
  plot(twodim_coord[[i]], main = paste("Paratope amino acids of",abs[i]), cex = 1.6,bty = "n", axes = F,xlab = "",ylab = "", ylim = range(twodim_coord[[i]][,2])+c(0,7), xlim = range(twodim_coord[[i]][,1])+c(0,7))
  points(twodim_coord[[i]][cnt,],pch = 16, col = "red",cex = 2)
  points(twodim_coord[[i]][cntm,],pch = 1, col = "blue",cex = 2.6,lwd = 2)
  #text(twodim_coord[[i]],labels = allaas[[i]], cex = 0.8)
  text(x = (min(twodim_coord[[i]][,1])+3.5),y = max(twodim_coord[[i]][,2])+2, labels = paste("mimotope:\n",bestpl[i]),cex = 1.3)
  legend("topright",bty = "n",title = "Contact amino acids",legend = c("with epitope","with mimotope"),col = c("red","blue"),pt.cex = c(2,2.6),pch = c(16,1),cex = 1.2)
  
  #text(twodim_coord[[i]],labels = allaas[[i]], cex = 0.8)
}
i =2

# for best of leiden clusters
mimo_contactsl = lapply(abs,\(ab){
  read.delim(file = paste0("dock_motifier_mimotopes/best_haddock_models_leiden/",ab,"_best_leiden_top1_1.pdb_contacts.txt"))
})
mimo_contactsl = lapply(mimo_contactsl,\(ab){
  ab$AB_contacts = lapply(ab$AB_contacts,\(c_list){aas = unlist(strsplit(c_list, split = ","));sapply(aas,\(a) substr(a,1,nchar(a)-1))})
  ab
})
names(mimo_contactsl) = abs
bestpl = c("WHPKDPD","RQIYSLE","HGCQMVS","NREPVID")
names(bestpl) = abs
par(mfrow = c(2,2))
for(i in  abs){
  cnt = which(allaas[[i]] %in% unlist(epi_contacts[[i]]$AB_contacts))
  cntm = which(allaas[[i]] %in% unlist(mimo_contactsl[[i]]$AB_contacts))
  #,xlab = "dimension one",ylab = "dimension two"
  plot(twodim_coord[[i]], main = paste("Paratope amino acids of",i), cex = 1.6,bty = "n", axes = F,xlab = "",ylab = "", ylim = range(twodim_coord[[i]][,2])+c(0,7), xlim = range(twodim_coord[[i]][,1])+c(0,7))
  points(twodim_coord[[i]][cnt,],pch = 16, col = "red",cex = 2)
  points(twodim_coord[[i]][cntm,],pch = 1, col = "blue",cex = 2.6,lwd = 2)
  #text(twodim_coord[[i]],labels = allaas[[i]], cex = 0.8)
  text(x = (min(twodim_coord[[i]][,1])+3.5),y = max(twodim_coord[[i]][,2])+2, labels = paste("mimotope:\n",bestpl[i]),cex = 1.3)
  legend("topright",bty = "n",title = "Contact amino acids",legend = c("with epitope","with mimotope"),col = c("red","blue"),pt.cex = c(2,2.6),pch = c(16,1),cex = 1.2)
}
i = 4

# for worst leiden clusters

mimo_contactslw = lapply(abs,\(ab){
  read.delim(file = paste0("dock_motifier_mimotopes/best_haddock_models_leiden/",ab,"_worst_leiden_top1_1.pdb_contacts.txt"))
})
mimo_contactslw = lapply(mimo_contactslw,\(ab){
  ab$AB_contacts = lapply(ab$AB_contacts,\(c_list){aas = unlist(strsplit(c_list, split = ","));sapply(aas,\(a) substr(a,1,nchar(a)-1))})
  ab
})
names(mimo_contactslw) = abs

worstpl = c("SINHGLM","QSHPDLF","SSSNYFA","LFLRTYL")
names(worstpl) = abs
for(i in  abs){
  cnt = which(allaas[[i]] %in% unlist(epi_contacts[[i]]$AB_contacts))
  cntm = which(allaas[[i]] %in% unlist(mimo_contactslw[[i]]$AB_contacts))
  #,xlab = "dimension one",ylab = "dimension two"
  plot(twodim_coord[[i]], main = paste("Paratope amino acids of",i), cex = 1.6,bty = "n", axes = F,xlab = "",ylab = "", ylim = range(twodim_coord[[i]][,2])+c(0,7), xlim = range(twodim_coord[[i]][,1])+c(0,7))
  points(twodim_coord[[i]][cnt,],pch = 16, col = "red",cex = 2)
  points(twodim_coord[[i]][cntm,],pch = 1, col = "blue",cex = 2.6,lwd = 2)
  #text(twodim_coord[[i]],labels = allaas[[i]], cex = 0.8)
  text(x = (min(twodim_coord[[i]][,1])+3.5),y = max(twodim_coord[[i]][,2])+2, labels = paste("mimotope:\n",worstpl[i]),cex = 1.3)
  legend("topright",bty = "n",title = "Contact amino acids",legend = c("with epitope","with mimotope"),col = c("red","blue"),pt.cex = c(2,2.6),pch = c(16,1),cex = 1.2)
}

# big plot with best vs worst peptides
par(mfrow = c(4,2), mar = c(0.5,3,3,0))
mm = list(mimo_contactsl,mimo_contactslw)
mimos = list(bestpl,worstpl)
for(i in  abs){
  for(j in 1:2){
    cnt = which(allaas[[i]] %in% unlist(epi_contacts[[i]]$AB_contacts))
  cntm = which(allaas[[i]] %in% unlist(mm[[j]][[i]]$AB_contacts))
  #,xlab = "dimension one",ylab = "dimension two"
  plot(twodim_coord[[i]], cex = 1.5,bty = "n", axes = F,xlab = "",ylab = "", ylim = range(twodim_coord[[i]][,2])+c(0,7), xlim = range(twodim_coord[[i]][,1])+c(0,7))
  points(twodim_coord[[i]][cnt,],pch = 16, col = "red",cex = 2)
  points(twodim_coord[[i]][cntm,],pch = 1, col = "blue",cex = 2.6,lwd = 2)
  #text(twodim_coord[[i]],labels = allaas[[i]], cex = 0.8)
  text(x = (min(twodim_coord[[i]][,1])+3.5),y = max(twodim_coord[[i]][,2])+2, labels = paste("mimotope:\n",mimos[[j]][i]),cex = 1.3)
  if(j ==2& i =="b12"){
    legend("topright",bty = "n",title = "Contact amino acids",legend = c("with epitope","with mimotope"),col = c("red","blue"),pt.cex = c(2,2.6),pch = c(16,1),cex = 1.2)

  }
  if(j == 1 ){
    mtext(i,side = 2, line = 1)
  }
  if(i=="Herceptin" & j==1){
    title(main = "High sequence similarity \n to epitope",cex.main = 1.3)
  }
  if(i=="Herceptin" & j==2){
    title(main = "Low sequence similarity \n to epitope",cex.main = 1.3)
  }
  }
  }


# second best structure
mimo_contactsl2 = lapply(abs,\(ab){
  read.delim(file = paste0("dock_motifier_mimotopes/best_haddock_models_leiden/",ab,"_best_leiden_top2_1.pdb_contacts.txt"))
})
mimo_contactsl2 = lapply(mimo_contactsl2,\(ab){
  ab$AB_contacts = lapply(ab$AB_contacts,\(c_list){aas = unlist(strsplit(c_list, split = ","));sapply(aas,\(a) substr(a,1,nchar(a)-1))})
  ab
})
names(mimo_contactsl2) = abs

mimo_contactslw2 = lapply(abs,\(ab){
  read.delim(file = paste0("dock_motifier_mimotopes/best_haddock_models_leiden/",ab,"_worst_leiden_top2_1.pdb_contacts.txt"))
})
mimo_contactslw2 = lapply(mimo_contactslw2,\(ab){
  ab$AB_contacts = lapply(ab$AB_contacts,\(c_list){aas = unlist(strsplit(c_list, split = ","));sapply(aas,\(a) substr(a,1,nchar(a)-1))})
  ab
})
names(mimo_contactslw2) = abs

par(mfrow = c(2,2), mar = c(0.5,3,3,0))
mm = list(mimo_contactsl2,mimo_contactslw2)
mimos = list(bestpl,worstpl)
for(i in  abs[2:3]){
  for(j in 1:2){
    cnt = which(allaas[[i]] %in% unlist(epi_contacts[[i]]$AB_contacts))
    cntm = which(allaas[[i]] %in% unlist(mm[[j]][[i]]$AB_contacts))
    #,xlab = "dimension one",ylab = "dimension two"
    plot(twodim_coord[[i]], cex = 1.5,bty = "n", axes = F,xlab = "",ylab = "", ylim = range(twodim_coord[[i]][,2])+c(0,7), xlim = range(twodim_coord[[i]][,1])+c(0,7))
    points(twodim_coord[[i]][cnt,],pch = 16, col = "red",cex = 2)
    points(twodim_coord[[i]][cntm,],pch = 1, col = "blue",cex = 2.6,lwd = 2)
    #text(twodim_coord[[i]],labels = allaas[[i]], cex = 0.8)
    text(x = (min(twodim_coord[[i]][,1])+3.5),y = max(twodim_coord[[i]][,2])+2, labels = paste("mimotope:\n",mimos[[j]][i]),cex = 1.3)
    if(j ==2& i =="21c"){
      legend(x = (max(twodim_coord[[i]][,1])-12),y = min(twodim_coord[[i]][,2])-2,bty = "n",title = "Contact amino acids",legend = c("with epitope","with mimotope"),col = c("red","blue"),pt.cex = c(2,2.6),pch = c(16,1),cex = 1.1)
      
    }
    if(j == 1 ){
      mtext(i,side = 2, line = 1)
    }
    if(i=="21c" & j==1){
      title(main = "High sequence similarity \n to epitope",cex.main = 1.3)
    }
    if(i=="21c" & j==2){
      title(main = "Low sequence similarity \n to epitope",cex.main = 1.3)
    }
  }
}


# contact overlap between top1 and top2 conformation as Jaccard Index

for(i in 1:4){
  nx = length(unique(c(unlist(mimo_contactsl[[i]]$AB_contacts),unlist(mimo_contactsl2[[i]]$AB_contacts))))
  print(length(intersect(unlist(mimo_contactsl[[i]]$AB_contacts),unlist(mimo_contactsl2[[i]]$AB_contacts)))/nx)
}

for(i in 1:4){
  nx = length(unique(c(unlist(mimo_contactslw[[i]]$AB_contacts),unlist(mimo_contactslw2[[i]]$AB_contacts))))
  print(length(intersect(unlist(mimo_contactslw[[i]]$AB_contacts),unlist(mimo_contactslw2[[i]]$AB_contacts)))/nx)
}
for(i in 1:4){
  nx = length(unique(c(unlist(mimo_contactslw[[i]]$AB_contacts),unlist(mimo_contactsl[[i]]$AB_contacts))))
  print(length(intersect(unlist(mimo_contactslw[[i]]$AB_contacts),unlist(mimo_contactsl[[i]]$AB_contacts)))/nx)
}

# contact overlap btw best and epi as prcentage of epitope contacts

for(i in 1:4){
  nx = length(unique(unlist(epi_contacts[[i]]$AB_contacts)))
  print(length(intersect(unlist(mimo_contactsl[[i]]$AB_contacts),unlist(epi_contacts[[i]]$AB_contacts)))/nx)
}
for(i in 1:4){
  nx = length(unique(unlist(epi_contacts[[i]]$AB_contacts)))
  print(length(intersect(unlist(mimo_contactslw[[i]]$AB_contacts),unlist(epi_contacts[[i]]$AB_contacts)))/nx)
}

