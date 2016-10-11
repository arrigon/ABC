args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }

require(e1071)


### Script IO : TO COMMENT UNLESS FOR DEBUG
# prepare input/output paths
# input = "../Qoutputs/Run1/simulation_g175.dat.txt"
# output = "../Qstats/SummaryStats_simulation.txt"
# ini_generations = 50
# N_aegil = 320
# s_aegil = 0.95
# shape = 500
# magnitude = .5
# nloc = 100


### import dataset
# data = read.delim(input, skip = nloc + 1, sep = "\t", header = F) #used when R deals directly with QNEMO outputs
data = read.delim(input, sep = "\t", header = F)

### extract info from data
# get population column
pops = data[, 1]

# count number of available pops
nsimpop = length(table(pops))

# get loci columns
mat01 = data[, 2:(ncol(data) - 1)] 


### Subsampling (NOW TAKEN CARE BY QNEMO and perl: QNEMO takes real + XX indivs, so that we can purge out wheat immigrants and refine subsampling with perl script). 
# Here, we just make checks to ensure that the available data meets sampling requirements
# get sampling efforts from AFLP/inifile.
real = read.delim("../data/Ini_dist_to_fields.txt", header = T)
npop = nrow(real) + 1 #add wheat in counting pops

# check that nb avail pops = nb required pops
if(npop != nsimpop){ #several pops were lost, skip that simulation
  ## Prepare and save output
  summary_stats = rep(NA, 5 * (npop - 1))
  summary_params = c(N_aegil, s_aegil, shape, magnitude)
  OUT = c(summary_params, summary_stats)
  cat(OUT, "\n", file = output, sep = "\t", append = T)
  exit()
  }

# If this is OK, proceed to next check
samp.aeg = real[, 3] # sampling effort Aegilops
if(table(pops)[1] > 100){ # sampling effort Wheat
  samp.wh = 100
  } else {
  samp.wh = table(pops)[1]
  }
samp.effort = c(samp.wh, samp.aeg)

# compare what is available to samp.effort; skip estimation if discrepency > 25% effectif of a single pop
avail = table(pops)
discr_pct = (avail - samp.effort) / samp.effort

missn = discr_pct[discr_pct < 0] # find any pops where samples are missing

if(length(missn) > 0){ #if there are any missing specimens, do the following
  if(any(missn < -.25)){ #check that missing specimens are not exceeding 25% of a given pop, skip that simulation if that happens
    ## Prepare and save output
    summary_stats = rep(NA, 5 * (npop - 1)) #actually replace results with NA data
    summary_params = c(ini_generations, N_aegil, s_aegil, shape, magnitude)
    OUT = c(summary_params, summary_stats)
    cat(OUT, "\n", file = output, sep = "\t", append = T)
    exit()
    } else { # Otherwise, adjust subsampling number to what is available (prevent code crash)
    correct = which(discr_pct < 0)
    samp.effort[correct] = avail[correct]
    }
  }


# recode diploid-codominant -> dominant. NOW TAKEN CARE BY PERL SCRIPT


### Compute summary stats
# Prune non-informative bands
centers = pops
centers[centers > 1] = 0
sp.freqs = aggregate(mat01, by = list(centers), mean)
sp.freqs = sp.freqs[, -1]
deltafreqs = abs(sp.freqs[1, ] - sp.freqs[2, ])
bestbands = which(deltafreqs > 0.15) #using 0.15, as in empirical data

if(length(bestbands) < 5){ #if there are any missing specimens, do the following
  ## Prepare and save output
  summary_stats = rep(NA, 5 * (npop - 1)) #actually replace results with NA data
  summary_params = c(ini_generations, N_aegil, s_aegil, shape, magnitude)
  OUT = c(summary_params, summary_stats)
  cat(OUT, "\n", file = output, sep = "\t", append = T)
  exit()
  }

mat01 = mat01[, bestbands]
nbands = ncol(mat01) / 100


# provide centroid coordinates (enforce specimen assignment in cmeans, also based on frequencies, done here to minimize computational load)
sp.centroids = sp.freqs[, bestbands]


## Compute cmeans admixtures
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

# compute admixtures 0 = Aegilops, 1 = Wheat
aeg.fuzz = search.fuz(mat01, sp.centroids, 1.25, 10) ## WARNING using 1.25 to stick with empirical runs
aeg.test = aeg.fuzz[ pops == 1,] ## previous version had pops > 1
aeg.fuzz = aeg.fuzz[, which.min(colSums(aeg.test))] ## previous version had which.max
summary_cmeans = aggregate(aeg.fuzz, by = list(pops),  mean)[, -1] #we get the mean
summary_cmeanssd = aggregate(aeg.fuzz, by = list(pops),  sd)[, -1] #we get the SD
summary_cmeans = data.frame(summary_cmeans, summary_cmeanssd)
summary_cmeans = summary_cmeans[ -1,]  #get rid of first line (wheat patch) and first column (pop names)
summary_cmeans = unlist(summary_cmeans)


### Compute Allele frequencies by pop: mean and SD
all.freqs = aggregate(mat01, by = list(pops),  mean)[,- 1] #get rid of first column (pop names)
summary_mean = rowMeans(all.freqs)
summary_sd = apply(all.freqs, 1, sd)
summary_mean = data.frame(summary_mean, summary_sd)
summary_mean = summary_mean[-1 ,]  #get rid of first line (wheat patch) and first column (pop names)
summary_mean = unlist(summary_mean)


### Prepare and save output
### output order: freq.mean, freq10, freq50, freq90, cmeans.mean, cmeans10, cmeans50, cmeans90
summary_stats = c(rep(nbands, length(summary_mean) / 2), summary_mean, summary_cmeans)
summary_params = c(ini_generations, N_aegil, s_aegil, shape, magnitude)
OUT = c(summary_params, summary_stats)
cat(OUT, "\n", file = output, sep = "\t", append = T)

