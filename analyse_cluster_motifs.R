
# compare kld to kld from the random clusters, calculate z-scores, make graphics
library(stats)
library(gplm)


# get random results
gsizes = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  vcount(G)
})
load("rnd_len_kld_3.RData" )
par(mfrow = c(2,2))
for(i in 1:4){
  plot(rnd_len_kld, main = abs[i],log =  "xy",xlab =  "", col = adjustcolor( "grey",0.6), xlim = range(lengths),ylim = c(min(rnd_len_kld[,2]),max(all_kld)))
  points(lengths[numcodes==i], all_kld[numcodes==i],col =  "red",log =  "xy",pch = 16)
  if(i>2){mtext("cluster size", side = 1, line = 3)}
  
}

real_len = lapply(db_pepsets,\(sets) sapply(sets, length))
for(i in 1:4){names(real_len[[i]]) = names(logpepsets[[i]])}
  

data = data.frame(log10(rnd_len_kld[,c(3,2)]))
  
  model=kreg(data$un.l,data$kld,grid = data$un.l)
  ii = order(data$un.l)
  res = (model$y-data$kld[ii])
  plot(model$x,res,col = adjustcolor("black",alpha.f = 0.2),pch = 16,main = paste(abs[i], "log"))
  colnames(data) = c( "l", "kld")

# calculate z scores from distribution of the kld of random clusters----
kl = all_kld
loglen = log10(unlist(real_len))
len = unlist(real_len)
colcodes = numcodes
par(mfrow = c(2,2))
z_sc = lapply(1:4,\(i){
  
  real_l = loglen[which(colcodes==i)]
  names(real_l) = names(len[which(colcodes==i)])
  real_k = log10(kl[which(colcodes==i)])
  names(real_k)= names(kl[which(colcodes==i)])
 
  pred =kreg(data$l,data$kld,grid = real_l,bandwidth = model$bandwidth)
  ii = order(data$l)
  sq_res = res**2
  #hist(sq_res)
  mod_res = kreg(model$x,sq_res,grid = data$l)
  pred_res = kreg(model$x,sq_res,grid = real_l,bandwidth = mod_res$bandwidth) # var predictions for real data
  sd = sqrt(pred_res$y)
  
  ij = order(real_l)
  plot((data), col =  "grey",xlim = c(min(loglen), max(model$x)),ylim = c(min((data$kld)),max(real_k)))
  points(model$x, model$y + 2*sqrt(mod_res$y), col =  "pink")
  points(model$x, model$y - 2*sqrt(mod_res$y), col =  "pink")
  points(pred,col =  "red",pch = 16)
  points((real_l), log10(all_kld[numcodes==i]),col =  "blue",pch = 16)
  z = as.vector((real_k[ij]-pred$y)/sd)
  
  big = which(is.na(pred$y)) # out of the range of the randomly generated clusters
  l = length(model$x)
  
 
  m = which.max(pred$x[-big])
  z[big] = as.vector((real_k[ij][big]-pred$y[-big][m])/sd[-big][m])
  names(z) = names(sort(real_l))
  z
  
})


# p- values
allp = sapply(unlist(z_sc),\(z){
  pnorm(z,lower.tail = F)
})
allp = p.adjust(allp,method = "BH")

# mthesis plots kld vs length
par(mfrow = c(2,2))
for(i in 1:4){
  plot((data), col =  "grey",xlim = c(min(model$x), max(model$x)),ylim = c(min((data$kld)),max(real_k)), xlab =   "", ylab =  "", main = abs[i])
  points(log10(real_len[[i]]), log10(all_kld[numcodes==i]),col =  "red",pch = 16)
  pp = names(which(allp[numcodes==i]>0.01))
  points(log10(real_len[[i]])[pp], log10(all_kld[numcodes==i])[pp],col =  "blue",pch = 16)
  
  if(i>2){mtext("log10(cluster size)", side = 1, line = 3)}
  if(i%%2==1){mtext("log10(KLD)", side = 2, line = 3)}
  if(i==4){
    legend( "topright", pch = c(1,16,16), col = c( "grey", "blue", "red"), legend = c( 'control', 'p >= 0.01', 'p < 0.01' ),bty = "n")
  }
}

# kld vs length leiden vs dbscan clusters -----
load(file =  "unique_leiden_pepsets.RData")
### plot kld ----

kld = function(p,q){
  sum(p*log(p/q))
}
leid_logfreqs = future_lapply(1:4,\(i){
  require(igraph)
  ab = abs[i]
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7or.RData"))
  freqs= V(G)$Freq
  names(freqs) = V(G)$name
  
  sapply(leiden_pepsets_un[[i]],\(set) round(log2(freqs[set])) )})

leid_logpepsets = lapply(1:4,\(i){
  lapply(1:length(leid_logfreqs[[i]]),\(j){
    
    rep(leiden_pepsets_un[[i]][[j]],leid_logfreqs[[i]][[j]])
  })
})
for(i in 1:4){names(leid_logpepsets[[i]]) = names(leiden_pepsets_un[[i]])}
llen = lapply(leid_logpepsets,\(sets) sapply(sets, length))

leid_ppm = lapply(leid_logpepsets,\(sets){ lapply(sets,\(cl) create_motif(cl, alphabet = "AA", pseudocount = 1)@motif)})
leid_kld = lapply(leid_ppm,\(sets){ sapply(sets, kld, bgmot)})

## scatterplot leiden and dbscan together
for(i in 1:4){
  plot((data), col =  "grey",xlim = c(log10(5), max(model$x)),ylim = c(min((data$kld)),log10(max(unlist(leid_kld))+1)), xlab =   "", ylab =  "", main = abs[i])
  points(log10(lengths[numcodes==i]), log10(all_kld[numcodes==i]),col =  adjustcolor("blue", 0.7),pch = 16)
  # pp = names(which(allp[numcodes==i]>0.01))
  # points(log10(lengths[numcodes==i])[pp], log10(all_kld[numcodes==i])[pp],col =  "blue",pch = 16)
  points(log10(llen[[i]]), log10(leid_kld[[i]]),col =  adjustcolor("pink",0.7),pch = 16)
  if(i>2){mtext("log10(cluster size)", side = 1, line = 3)}
  if(i%%2==1){mtext("log10(KLD)", side = 2, line = 3)}
  if(i==4){
    legend( "topright", pch = c(1,16,16), col = c( "grey", "blue", "red"), legend = c( 'control', 'p >= 0.01', 'p < 0.01' ),bty = "n")
  }
}

## with number of unique peptides
real_len = lapply(db_pepsets,\(sets) sapply(sets, length))

for(i in 1:4){
  plot(log10(real_len[[i]]),log2(all_kld[numcodes==i]),main = abs[i], xlab =  "", ylab =  "",col = adjustcolor("blue",0.8),pch = 16, xlim = log10(range(unlist(c(real_len,len)))), ylim = log2(range(unlist(c(all_kld,leid_kld)))))
  points(log10(len[[i]]),log2(leid_kld[[i]]), col = adjustcolor("pink",0.7), pch = 16) 
  mtext(side = 1,"log10(size)", line = 3)
  if(i == 1){
    mtext(side = 2,"log2(KLD)", line = 3)
    legend( "bottomleft",legend = c( "Leiden",  "DBSCAN"),bty =  "n", cex = 1.2, pch = 16, col = c(adjustcolor("pink",0.7),adjustcolor("blue",0.8)))
    }
}

## bowplot for leiden vs dbscan kld 
par(mfrow = c(1,4), mar = c(3,3,2,0))
for(i in 1:4){
  
  boxplot(list(log2(leid_kld[[i]]),log2(all_kld[numcodes==i])), boxwex = 0.6, col = c( "pink", "lightblue") ,xlab =  "", ylab =  "" ,axes = F, main = abs[i], ylim = log2(c(min((all_kld)), max(unlist(leid_kld))+12)))
  if(i == 1){axis(2); mtext(side = 2,"log2(KLD)", line = 3)}
  axis(1, labels = c("Leiden","DBSCAN"),at = c(1,2),tick = FALSE, lwd = 0, cex.axis = 1.4)
  
  p = wilcox.test(log2(leid_kld[[i]]),log2(all_kld[numcodes==i]))$p.value
  print(p)
  print(sum(ib))
  if(p<0.05/4){
    points(1.5,log2(30),pch = 8)
    lines(c(1,2),log2(c(26,26)))
    
  }
}

