### This script computes the empirical summary statistics
# based on the complete AFLP dataset
# using only the subset of 458 + 33 specimens included in the study
# we compute the summary stats from a subsampling of 100 AFLP bands that are highly informative.
# To do so we:
# 1. Select 100 bands to compute the empirical summary statistics (the simulations are based on 100 loci)
# 2. From these, we further select AFLP bands that differ in frequencies among the wheat and Aegilops pools (difference >= 15%)
# 3. We compute summary statistics based on this sample of AFLP bands
# 4. We repeat the step 1 -> 3, 100 times to get averaged summary statistics
# As a result, we get summary stats that are representative of the AFLP dataset, 
# but based on the most informative bands among 100 randomly selected loci
# Note that the same procedure is followed to compute summary stats from the simulations


### Functions
# function for merging tables
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

# fuzzy c-means
require(e1071)
search.fuz = function(matm, grp, fuz, nreps){
  res = NULL
  werr = NULL
  for(i in 1:nreps){
    c.fuz = cmeans(matm, grp, m=fuz)
    res = rbind(res, data.frame(i, c.fuz$membership))  
    werr = c(werr, c.fuz$withinerror) 
    }
  best = which.min(werr)[1]
  res[res[,1] == best,-1]
  }


## Open datasets
data = read.delim("TGCTAATCBinarytest.txt",header = T, row.names = 1)
matinfo = read.delim("FocalSpecimens.txt", header = T)

# recode species column
matinfo.wh = matinfo[ matinfo$TypAeg != 'triun',]
matinfo.ae = matinfo[ matinfo$TypAeg == 'triun',]
matinfo = rbind(matinfo.wh, matinfo.ae)
matinfo[ matinfo$TypAeg != 'triun', ]$TypAeg = "Taes"

# join metadata and AFLP tables, and extract 01 matrix of it
x = merging(matinfo, data, matinfo$Tag, rownames(data))
matm = x[, (ncol(matinfo) + 1):ncol(x)]
rownames(matm)=x$Tag
pops=as.character(x$Pop)


# Recode / swap alleles, so that frequencies of wheat are > than those of Aegilops
# 1. find which loci need swapping
freqs = aggregate(matm, by = list(x$TypAeg), mean)
rownames(freqs) = freqs[, 1]
freqs = freqs[, -1]
focus = which(freqs[1, ] < freqs[2, ])

# 2. perform swapping
matm_swap =  swapfun(matm, focus)

# 3. check job
freqs_swap = aggregate(matm_swap, by = list(x$TypAeg), mean)
rownames(freqs_swap) = freqs_swap[, 1]
freqs_swap = freqs_swap[, -1]
boxplot(t(freqs_swap))

require(RawGeno)
ACOP(dist(matm), x$Aeg)
ACOP(dist(matm_swap), x$Aeg)


## Identify non-informative bands (from complete AFLP dataset): proceed as in summary stats from simulations, based on deltafreqs
# this part of the script checks that these bands are really informative
# there is no summary stats done here.
deltafreqs = abs(freqs_swap[1 ,] - freqs_swap[2 ,]) #WARNING: had bug here.
bestbands = colnames(deltafreqs[, deltafreqs > 0.15]) #deltafreq = 0.15 as in SummaryStats_V5.r

# check that job is correct: SUPER IMPORTANT
keptbands = which(!is.na(match(colnames(matm_swap), bestbands)))
matm_swap_best = matm_swap[, keptbands]

require(RawGeno)
ACOP(dist(matm_swap_best), x$Aeg) # OK, we keep a good wheat-Aegilops differentiation

# check allele frequencies of these bands, at species level
freqs_swap_best = aggregate(matm_swap_best, by = list(x$TypAeg), mean)
rownames(freqs_swap_best) = freqs_swap_best[, 1]
freqs_swap_best = freqs_swap_best[, -1]
checkdelta = abs(freqs_swap_best[1, ] - freqs_swap_best[2, ])
range(checkdelta)
boxplot(t(freqs_swap_best)) # still good

# check allele frequencies of these bands, at pop level
x$Pop = as.character(x$Pop)
x[ x$TypAeg == "Taes", 1] = "Wheat"
freqs_final_check = aggregate(matm_swap_best, by = list(x$Pop), mean)
rownames(freqs_final_check) = freqs_final_check[, 1]
freqs_final_check = freqs_final_check[, -1]
freqs_final_check = t(freqs_final_check)
boxplot(freqs_final_check)
image(t(freqs_final_check)) # lines = AFLP bands, columns = pops. 
			    # you can see that wheat accessions have high frequencies (white / yellow cells),
			    # Aegilops pops have low frequencies (red / orange cells)
rownames(freqs_final_check)


## Compute summary statistics
# Check 1: visualize allele frequencies in each pop (still using all the available AFLP bands)
x$Pop = as.character(x$Pop)
x[ x$TypAeg == "Taes", 1] = "Wheat"

freqs_final = aggregate(matm_swap, by = list(x$Pop), mean)
rownames(freqs_final) = freqs_final[, 1]
freqs_final = t(freqs_final[, -1])
boxplot(freqs_final)


# Mean and SD of wheat-diagnostic markers
# This must be done on a subset of 100 loci, in order to stick to what is done in the simulations
OUT = NULL
for(i in 1:100){
  # get 100 loci from the complete AFLP dataset
  freqs_100 = freqs_final[sample(1:nrow(freqs_final), 100, replace = F), ]

  # keep only the informative ones (identified earlier)
  freqs_it = freqs_100[ !is.na(match(rownames(freqs_100), bestbands)),]

  # compute mean and SD of allele frequencies (corresponds to prevalence of wheat diagnostic bands)
  summary_Mmean = apply(freqs_it, 2, mean)
  summary_SDmean = apply(freqs_it, 2, sd)
  summary_mean = data.frame(names(summary_Mmean), summary_Mmean, summary_SDmean, nrow(freqs_it) / 100)
  OUT = rbind(OUT, summary_mean)
  }
summary_mean = aggregate(OUT[, -1], by = list(OUT[, 1]), median)
rownames(summary_mean) = summary_mean[, 1]
summary_mean = summary_mean[, -1]


# Mean and SD of pop-level admixtures (using c-means)
# Subsample 100 loci from there
OUT = NULL
for(i in 1:100){
  matm_swap_100 = matm_swap[, sample(1:ncol(matm_swap), 100, replace = F)]
  matm_swap_sub = matm_swap_100[, !is.na(match(colnames(matm_swap_100), bestbands))]
  
  # Place it on species centers
  centers = aggregate(matm_swap_sub, by = list(x$TypAeg), mean)[, -1]
  aeg.fuzz = search.fuz(matm_swap_sub, centers, 1.25, 10) #1.25 works just fine here, same as used in simulations
  aeg.test = aeg.fuzz[ x$TypAeg == "Taes",]
  aeg.fuzz = aeg.fuzz[, which.min(colSums(aeg.test))]
  plot(aeg.fuzz, x$Aeg, main = ncol(matm_swap_sub)) #check we get consistent admixture values

  # Aggregate at pops level
  summary_Mcmeans = aggregate(aeg.fuzz, by = list(x$Pop), mean)
  summary_SDcmeans = aggregate(aeg.fuzz, by = list(x$Pop), sd)
  summary_cmeans = data.frame(summary_Mcmeans[[1]], unclass(summary_Mcmeans[[2]]), unclass(summary_SDcmeans[[2]]))

  OUT = rbind(OUT, summary_cmeans)
  }
summary_cmeans = aggregate(OUT[, -1], by = list(OUT[, 1]), mean)
mean(OUT[,2] > .45 & OUT[,2] < .55) #proportion of failed runs


## Prepare and save outputs
OUT = data.frame(summary_mean, summary_cmeans)
OUT = OUT[ rownames(OUT) != "Wheat", ]

# sort it as in Ini_dist_to_fields.t
dists.geo = read.delim("../data/Ini_dist_to_fields.txt",header=T)
OUT = merging(dists.geo, OUT, dists.geo[, 1], rownames(OUT))

# quick visual check: verifying that Aegilops pops
# that are close to cultivations have higher admixture values
par(mfrow = c(2,2), pty = 's')
boxplot(OUT[,4] ~ cut(log(OUT$Dist + 1), 2))
boxplot(OUT[,5] ~ cut(log(OUT$Dist + 1), 2))
boxplot(OUT[,8] ~ cut(log(OUT$Dist + 1), 2))
boxplot(OUT[,9] ~ cut(log(OUT$Dist + 1), 2))


summary_stats = unlist(OUT[, c(6, 4:5, 8:9)])

n_patches = nrow(OUT)

output = "SummaryStats_empirical_mean_sd_100loci_pruning0.15_R1bis.25.txt"
param_names = c("N_aegil", "s_aegil", "shape", "magnitude")
band_names = paste("P.bands", 1:n_patches, sep = "")
meanM_names = paste("freq.mean", 1:n_patches, sep = "")
meanSD_names = paste("freq.sd", 1:n_patches, sep = "")
cmeansM_names = paste("cmeans.mean", 1:n_patches, sep = "")
cmeansSD_names = paste("cmeans.sd", 1:n_patches, sep = "")
headers = c(band_names, meanM_names, meanSD_names, cmeansM_names, cmeansSD_names)
cat(headers ,"\n", sep = "\t", file = output, append = F)
cat(summary_stats, file = output, sep = "\t", append = T)



