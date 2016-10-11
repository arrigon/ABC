args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }



####### Parameters of the simulation COMMENT IF NOT USED FOR DEBUG
#### Parameters to optimize (N, s, m and d <-> c)
# setwd("/home/arrigon/Dropbox/Projets/01PROJECTS/01 - Aegilops/ABC/Qnemo/scripts/") 
# # Run number
# run = 10
# output = "../Qstats/SummaryStats_simulation_test.txt"
# 
# # carrying capacities
# N_aegil = 200 # N
# N_wheat = 5000
# 
# # selfing rate
# s_aegil = 0.95 # s
# 
# # dispersal params
# shape = 100 # c
# magnitude = 0.5 # m
# 
# # number of generations
# ini_generations = 100



#### Initial conditions / templates
# Template for parameters file
ini_template = "../data/Ini_template_file.ini"

# Initial allele frequencies
ini_allele_freqs = "../data/Ini_real_phen_freqs.txt"

# number of loci: get it from frequency file
tmp = read.delim(ini_allele_freqs, header = T, row.names = 1)
# n_loci = nrow(tmp)
n_loci = 100 

# Distance of pops to wheat fields (!!! The ordering of that specific file determines the ordering of pops in the simulation)
ini_dist_to_fields = "../data/Ini_dist_to_fields.txt"

# number of patches: get it from the distances file
tmp = read.delim(ini_dist_to_fields, header = T, row.names = 1)
n_patches = nrow(tmp) + 1 #including wheat patch
# n_patches = 2 # or specify it directly; This must be in line with dispersal_matrix! Also account for wheat patch!



##### IO DO NOT COMMENT
path_to_summstats = paste(output, sep = '')
RunFolder = paste("Run", run, sep = '')
path_to_ini =  paste("../Qparams/Ini_params_simulation", run, ".txt", sep = '')
path_to_freqs = paste("../Qparams/Ini_allele_freqs", run, ".txt", sep = '')
path_to_disp = paste("../Qparams/Ini_dispersal_matrix", run, ".txt", sep = '')



####### SIMULATION PIPELINE
## Produce Dispersal Probabilities (set c - shape and m - magnitude parameters)
# Uncomment according to desired model
# generate command line
command = paste("R CMD BATCH \"--args input='", ini_dist_to_fields,
				   "' output='", path_to_disp,
				   "' shape=", shape,
				   " magnitude=", magnitude, 
 				   "\" DispersalProba_Prim.r", sep = "") ## WARNING: Model 5: Primary Disp + Exponential Wheat Gene Flow
# 				   "\" DispersalProba_Prim_NULL.r", sep = "") ## WARNING: Model 1: Primary Disp + No Wheat Gene Flow
# 				   "\" DispersalProba_Prim_Uniform.r", sep = "") ## WARNING: Model 3: Primary Disp + Uniform Wheat Gene Flow
# 				   "\" DispersalProba_complex.r", sep = "") ## WARNING: Model 6: Secondary Disp + Exponential Wheat Gene FlowSecondary
#				   "\" DispersalProba_complex_NULL.r", sep = "") ## WARNING: Model 2: Secondary Disp + No Wheat Gene Flow
# 				   "\" DispersalProba_complex_Uniform.r", sep = "") ## WARNING: Model 4: Secondary Disp + Uniform Wheat Gene FlowSecondary

command
system(command)



## Produce initial allelic frequencies (defined using external file, npop must be given here)
# NOTE: npop = pops Aegilops + 1 pop wheat. Pop Nr1 = wheat 
command = paste("R CMD BATCH \"--args input='", ini_allele_freqs,
				   "' output='", path_to_freqs, 
				   "' npop=", n_patches, 
				   "\" Inifreq2.r", sep = "")
command
system(command) 



## Produce ini file of QuantiNEMO
path_to_disp = gsub("\\.\\./", "", path_to_disp)
path_to_freqs = gsub("\\.\\./Qparams/", "", path_to_freqs, perl = T)

command = paste("R CMD BATCH \"--args input='", ini_template,
				   "' output='", path_to_ini,
				   "' ini_generations=", ini_generations,
				    " path_to_disp='", path_to_disp,
				   "' path_to_freqs='", path_to_freqs,
				   "' n_patches=", n_patches,
				    " RunFolder='", RunFolder,
				   "' N_wheat=", N_wheat,
				    " N_aegil=", N_aegil,
				    " s_aegil=", s_aegil,
				    " n_loci=", n_loci,
                                    " \" UpdateIniFile.r", sep = "")
command
system(command) 



## Run QuantiNEMO
setwd("..")
path_to_ini = gsub("\\.\\./", "", path_to_ini, perl = T)
command = paste("./quantiNemo ", path_to_ini, sep = "")
command
system(command)



## Convert outputs quickly
setwd("scripts/")
command = paste("perl convert_counts.pl ../Qoutputs/", RunFolder, "/simulation_g", ini_generations, ".dat ../data/Ini_sampling_scheme.txt", sep = "")
command
system(command);



## Produce summary statistics
command = paste("R CMD BATCH \"--args input='../Qoutputs/", RunFolder,"/simulation_g", ini_generations,".dat.txt'",
				    " output='", path_to_summstats,"'",
				    " N_aegil=", N_aegil,
				    " s_aegil=", s_aegil,
				    " shape=", shape,
				    " magnitude=", magnitude,
				    " ini_generations=", ini_generations,
				    " nloc=", n_loci,
                                    " \" SummaryStats_V5.r", sep = "")
command
system(command)



## Clean mess
command = paste("rm ../Qoutputs/", RunFolder,"/simulation_g", ini_generations,".*", sep = "")
command
system(command)
