#### NB> Script to produce initial frequencies
# 1. Discards admixed specimens from dataset
# 2. Discards non-informative bands (deltafreq = 0.15): working only with bands having at least 15% of allele freq difference between wheat and Aegilops
# 3. Computes population-specific frequencies
# 4. One population was made essentially of admixed specimens,
#    -> we use the frequencies of its closest neighbor.


#### functions
merging = function(A, B, linkA, linkB){
  indsA = linkA
  indsB = linkB
  common = intersect(indsA, indsB)
  matA = A[ match(common, indsA),]
  matB = B[ match(common, indsB),]
  cbind(matA, matB)
  }


# function for swapping AFLP signal (1->0 and 0->1)
swapfun = function(x, focus){
  x[, focus] = 1 - x[, focus]
  x
  }


## Open Infos
data = read.delim("TGCTAATCBinary.txt",header=T,row.names=1)
matinfo = read.delim("FocalSpecimens.txt",header=T)

# recode species column and keep specimens that are well differentiated
matinfo.wh = matinfo[ matinfo$TypAeg != 'triun',]

# matinfo.ae = matinfo[ matinfo$Aeg>0.99,]
matinfo.ae = matinfo[ matinfo$TypAeg == 'triun',]
matinfo = rbind(matinfo.wh, matinfo.ae)
matinfo[ matinfo$TypAeg != 'triun', ]$TypAeg = "Taes"


### Join Both tables and extract 01 matrix of it
x = merging(matinfo, data, matinfo$Tag, rownames(data))


### Remove admixed specimens
x.clean = x[ x$TypAeg == "Taes" | x$Aeg > 0.99, ]
matm.clean = x.clean[, (ncol(matinfo) + 1):ncol(x.clean)]
rownames(matm.clean)=x.clean$Tag


### Recode / swap alleles, so that frequencies of wheat are > than those of Aegilops
# ensures we get "presences" prevailing in wheat and "absences" prevailing in Aegilops
# while preserving AFLP signal

# find which loci need swapping
freqs = aggregate(matm.clean, by = list(x.clean$TypAeg), mean)
rownames(freqs) = freqs[, 1]
freqs = freqs[, -1]
focus = which(freqs[1, ] < freqs[2, ])

#perform swapping
matm.clean =  swapfun(matm.clean, focus)


### Identify best bands
band.freqs = aggregate(matm.clean, by = list(x.clean$TypAeg), mean)
rownames(band.freqs) = band.freqs[, 1]
band.freqs = t(band.freqs[, -1])
deltafreqs = abs(band.freqs[, 1] - band.freqs[, 2])

# compute PCoA of individuals
require(vegan)
pca = rda(matm.clean)
eigen1 = scores(pca, display = "species")

# compare first eigenaxis vs deltafreqs
plot(eigen1[,1], deltafreqs)


### Select best bands and compute band frequencies in each pop (still with non-admixed specimens)
matm.clean.best = matm.clean[, deltafreqs > 0.15]

# compute freqs in pops and pool wheat specimens
pops = as.character(x.clean$Pop)
pops[x.clean$TypAeg == "Taes"] = "Wheat"

band.freqs = aggregate(matm.clean.best, by = list(pops), mean)
rownames(band.freqs) = band.freqs[, 1]
band.freqs = band.freqs[, -1]

# make sure that this table is sorted as Ini_dist_to_fields.txt
info = read.delim("../data/Ini_dist_to_fields.txt",header=T)
correctorder = c("Wheat", as.character(info$Pop))

band.freqs = band.freqs[ match(correctorder, rownames(band.freqs)),]
band.freqs[ rownames(band.freqs) == "824598",] = band.freqs[ rownames(band.freqs) == "849587",]

require(RawGeno)
ACOP(dist(band.freqs[-1,]), rownames(band.freqs[-1,]))

# Save these band frequencies
write.table(band.freqs, file = "../data/Ini_real_phen_freqs.txt", col.names = T, row.names = T, quote = F, sep = "\t")
write.table(colnames(band.freqs), file = "BandSelection.txt", row.names = F, col.names = F, quote = F)
