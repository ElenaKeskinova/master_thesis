# obtain aa freq distributions from the whole library

# only cyclic 10-mers are available
lib4 = readLines(file("D:/Elena/ban/Motifier_DataSet/4st_generation_library/4st_generation_random_library_Filtered_AA.fs"))
libnames = lib4[seq(1, length(lib4), by = 2)]
libpeps = lib4[seq(2, length(lib4), by = 2)]
tpeps = table(libpeps)
hist(log10(tpeps))
# remove cystein
peps_nC = sapply(libpeps, \(p){
  toupper(paste0(unlist(strsplit(p,""))[2:11],collapse = ""))
})

peps_7nC = sapply(libpeps, \(p){
  p = unlist(strsplit(p,""))[2:11]
  c(toupper(paste0(p[1:7],collapse = "")), toupper(paste0(p[4:10],collapse = "")))
})

peps_7nC = c(peps_7nC)
save(peps_7nC,file = "4st_generation_library/lib_allpeps_7nC.RData")
save(peps_nC,file = "4st_generation_library/lib_allpeps_nC.RData")

library(Biostrings)
freqM = consensusMatrix(peps_nC,as.prob = T)
freqM7 = consensusMatrix(peps_7nC,as.prob = T)
freqM7
kld(bgmot,freqM7)
freql = rowSums(freqM7)/7
rowSums(freqM)/10
bkg_m = create_motif(peps_7nC)
bgmot = bkg_m@motif
freq = rowSums(freqM)/ncol(freqM)

save(freqM,freqM7,file="4st_generation_library/lib_freqs.RData")
load(file="4st_generation_library/lib_freqs.RData")

