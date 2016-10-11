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


### Dispersal function (exponential family)
expfamily = function(dist, shape){
  1 / (2 * pi * shape^2) * exp(-(dist / shape))
  }


## Get distances of pops to cultivations
data = read.delim(input, header = T, row.names = 1)

# store them in a vector
x = data[, 1]


### Compute dispersal probabilities
# compute the dispersal proba, using these distances
disp.proba = expfamily(x, shape)

# standardize these proba, so that they sum to one (they will still scale identically)
disp.proba = magnitude * disp.proba / sum(disp.proba)

# prepare emigration matrix
disp.matrix = matrix(0, length(x) + 1, length(x) + 1) # make an empty (filled with zeroes) matrix of pop x pop dimensions

# prepare the first line of that matrix, which is defined to be the wheat patch
disp.line = c(1-sum(disp.proba), disp.proba)
disp.matrix[1 ,] = disp.line

# fill the diagonal with ones (except the wheat x wheat migration), so that Aegilops indivuals are not migrating anywhere
for(i in 2:nrow(disp.matrix)){
  disp.matrix[i,i] = 1
  }

# HACK: we want the NULL model, replace first line by no-migration values
disp.matrix[1, ] = c(1, rep(0, ncol(disp.matrix) - 1))

### Write outputs
# Write this table with the QuantiNEMO format
cat("{", file = output, append = F)

for(i in 1:nrow(disp.matrix)){
  line = disp.matrix[i, ]
  line.txt = paste(line, collapse = " ")
  cat("{ ", line.txt, "}\n", file = output, append = T)
  }

cat(" }\n", file = output, append = T)

