### Script parameters
## get arguments from bash
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }


# Uncomment ONLY for debugging this specific script.
# input = "../data/Ini_template_file.ini"
# output = "../Qparams/Ini_params_simulation.txt" # (could be given through bash), remains hard-coded for now
# path_to_disp = "../Qparams/Ini_dispersal_matrixXXXXX.ini"
# path_to_freq = "../Qparams/Ini_allele_freqsXXXXX.ini"
# RunFolder = "Run1"
# ini_generations = 50
# n_patches = 5
# N_wheat = 100
# N_aegil = 20
# s_aegil = 0.9
# n_loci = 100

### Function to update fields accordingly
# Find - Replace function
findrepl = function(x, targetfield, value){
  line_nr = grep(targetfield, x)
  x[line_nr] = paste(targetfield, value, sep = "\t")
  x
  }

ini_generations


### Modify Ini file
## Import template of Ini file
ini_file = readLines(input)

## Update fields
# Outputs folder
ini_file = findrepl(ini_file, "ntrl_genot_dir", RunFolder)

# N generations
ini_file = findrepl(ini_file, "generations", ini_generations)

# N generations in genotype logtimes
ini_file = findrepl(ini_file, "ntrl_genot_logtime", ini_generations)

# patch number
ini_file = findrepl(ini_file, "patch_number", n_patches)

# patch capacity Wheat is first cell and cap_Aegilops is repeated npatches-1 times
all_capacities = paste("{", N_wheat, " rep(", N_aegil, ", ", n_patches - 1, ")}", sep = "")
ini_file = findrepl(ini_file, "patch_capacity", all_capacities)

# dispersal matrix
path_to_disp = paste("$", path_to_disp, sep = "")
ini_file = findrepl(ini_file, "dispersal_rate", path_to_disp)

# Allele frequencies
ini_file = findrepl(ini_file, "ntrl_allelic_file", path_to_freqs)

# mating_proportion
ini_file = findrepl(ini_file, "mating_proportion", s_aegil)

# ntrl_loci
ini_file = findrepl(ini_file, "ntrl_loci", n_loci)


### Save the output
writeLines(ini_file, output)

