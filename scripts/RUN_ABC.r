args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }

### Set random seed (based on system time, ensures that different runs will have different random values)
# print with possibly greater accuracy:
op <- options(digits.secs=6)
hour = as.character(Sys.time())
seed = as.numeric(strsplit(hour, "\\:", perl = T)[[1]][3]) * 1e6
set.seed(seed)



### Parallel run params:
# run = 1 #COMMENT ME!!
# niter = 5 #ABC = 1000
niter = 2000 #ABC = 1000, Defines how many iterations are assigned per CPU


### Fixed parameters (overridden by PRIORS), uncomment for debugging purposes
N_wheat = 5000 #KEEP UNCOMMENTED #ABC = 5000
ini_generations = 50  #KEEP UNCOMMENTED

# selfing rate
# N_aegil = 100 # N
# s_aegil = 0.95 # s

# dispersal params
# shape = 3000 # c
# magnitude = 0.05 # m

# number of generations



### get patch number
ini_dist_to_fields = "../data/Ini_dist_to_fields.txt"
tmp = read.delim(ini_dist_to_fields, header = T, row.names = 1)
n_patches = nrow(tmp) #without accountig for wheat patch

### IO
output = paste("../Qstats/SummaryStats_simulation", run,".txt", sep = "")
param_names = c("N_generations", "N_aegil", "s_aegil", "shape", "magnitude")
band_names = paste("P.bands", 1:n_patches, sep = "")
meanM_names = paste("freq.mean", 1:n_patches, sep = "")
meanSD_names = paste("freq.sd", 1:n_patches, sep = "")
cmeansM_names = paste("cmeans.mean", 1:n_patches, sep = "")
cmeansSD_names = paste("cmeans.sd", 1:n_patches, sep = "")
headers = c(param_names, band_names, meanM_names, meanSD_names, cmeansM_names, cmeansSD_names)
cat(headers ,"\n", sep = "\t", file = output, append = F)



### Iterate all parameters
#############################
for(iter in 1:niter){
  # draw params from priors
  ini_generations = round(300*rbeta(1, 1, 2), 0)
  shape = runif(1, min = 1, max = 15000) #ABC: 1 - 15000
  N_aegil = runif(1, min = 300, max = 5000) #ABC: 300 - 5000 (cannot go below 300 because of subsampling issues)
  magnitude = runif(1, min = .001, max = 0.5) #MAX MUST BE 0.9 !!! (cannot go above 0.5 because of subsampling issues)
  s_aegil = runif(1, min = .05, max = 1) ### WARNING: Self = UNIF[0-1] here

  ini_generations
  shape
  N_aegil
  magnitude
  s_aegil

  # run pipeline
  command = paste("R CMD BATCH \"--args output='", output,
				    "' ini_generations=", ini_generations,
				      " run=", run,
				      " N_wheat=", N_wheat,
				      " N_aegil=", N_aegil,
				      " s_aegil=", s_aegil,
				      " shape=", shape,
				      " magnitude=", magnitude,
				      " \" Pilot.r", sep = "")
  system(command) 
  }
