### Script parameters
## get arguments from bash
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }



# Uncomment ONLY for debugging this specific script.
# shape = 1325 #shape parameter (given through bash)
# magnitude = 0.5 # must be between 0 and 1; should be < than 1, otherwise would empty wheat patch within one generation.
# input = "../data/Ini_dist_to_fields.txt"
# output = "../Qparams/Ini_dispersal_matrix.txt" # (could be given through bash), remains hard-coded for now

### priors Aegilops dispersal
magn_aeg = runif(1, 0, .5)


### Dispersal function (exponential family)
expfamily = function(dist, shape){
  1 / (2 * pi * shape^2) * exp(-(dist / shape))
  }


## Get distances of pops to cultivations
data = read.delim(input, header = T, row.names = 1)

## Get distances among Aegilops pops
aeg.path = "../data/Ini_dist_among_pops.txt"
aeg.d = read.delim(aeg.path, header = T, row.names = 1)

# store them in a vector
x = rbind(data[, 1], aeg.d)


### Compute dispersal probabilities
# prepare emigration matrix
# compute the dispersal proba, using these distances
# standardize these proba, so that they sum to one (they will still scale identically)
# also multiply by proportion of outmigrants
# mynorm = function(n, mean, sd = .05){
#   tmp = rnorm(n, mean, sd)
#   tmp[tmp < 0] = 1e-3
#   tmp[tmp > 1] = 0.999
#   tmp
#   }

# plot(density(mynorm(10000, .5, sd = 0.05)), xlim = c(0, 1), ylim = c(0, 50))
# for(i in seq(0, 1, length.out = 20)){
#   lines(density(mynorm(10000, i, 0.05)))
#   }

# magns = c(magnitude, mynorm(nrow(aeg.d), magn_aeg, sd = 0.05)) #with hyperprior
magns = c(magnitude, rep(magn_aeg, nrow(aeg.d))) #with cst prior -> keep that one

disp.matrix = NULL
for(i in 1:nrow(x)){
  disp = expfamily(x[i ,], shape)
  if(i > 1) disp[i - 1] = 0
  disp = magns[i] * disp / sum(disp)
  disp.matrix = rbind(disp.matrix, disp)
  }

# prepare the first column of that matrix (disp aeg to wheat, set to zero, except first line, wheat->wheat)
disp.matrix = cbind(c(1-sum(disp.matrix[1,]), rep(0, ncol(disp.matrix))), disp.matrix)

# fill the diagonal with 1-sum(outmigrants) (except the wheat x wheat migration), so that Aegilops indivuals are not migrating anywhere
for(i in 2:nrow(disp.matrix)){
  disp.matrix[i,i] = 1 - sum(disp.matrix[i,-i])
  }


### Write outputs
# Write this table with the QuantiNEMO format
cat("{", file = output, append = F)

for(i in 1:nrow(disp.matrix)){
  line = disp.matrix[i, ]
  line.txt = paste(line, collapse = " ")
  cat("{ ", line.txt, "}\n", file = output, append = T)
  }

cat(" }\n", file = output, append = T)

