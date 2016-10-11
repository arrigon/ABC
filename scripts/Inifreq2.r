args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }


### Uncomment ONLY for debugging this specific script.
# fake.wheat = seq(0, 1, length.out = 100)
# fake.aegil = seq(1, 0, length.out = 100)
# fake.all = rbind(fake.wheat, fake.aegil)
# input = "../data/FakeFreqs.txt"
# output = "../Qparams/Ini_allele_freqs.txt"
npop = npop - 1
mu = 1e-7 #mutation rate


### Open allelic frequency table (wheat must be in first row)
fake.all = read.delim(input, header = T, row.names = 1)


### Pick randomly 100 loci from that table
# Set random seed (based on system time, ensures that different runs will have different random values)
op <- options(digits.secs=6)
hour = as.character(Sys.time())
seed = as.numeric(strsplit(hour, "\\:", perl = T)[[1]][3]) * 1e6
set.seed(seed)

# subsample 100 loci from table
fake.all = fake.all[, sample(1:ncol(fake.all), 100, replace = F)]


### Convert dominant to codominant frequencies
# Using modification HWE for inbreeding, assuming diploidy
# where FreqPhen0 = q²s + qs
# so q²s + qs - FreqPhen0 = 0, and we solve for q
# then p = 1 - q
correction = function(freqPRESENCE, selfing = .9){
  freqABSENCE = 1 - freqPRESENCE
  a = 1 - selfing
  b = selfing
  c = -freqABSENCE
  x1 = (-b + sqrt(b^2 - 4 * a * c))/(2 * a)
  x2 = (-b - sqrt(b^2 - 4 * a * c))/(2 * a)
  outs = c(x1, x2)
  q = outs[which.max(outs)]
  p = 1 - q
  if(p < 0) p = 0
  if(p > 1) p = 1
  p
  }
fake.all.corr = matrix(NA, nrow(fake.all), ncol(fake.all))
for(i in 1:nrow(fake.all)){
  for(j in 1:ncol(fake.all)){
    fake.all.corr[ i, j] = correction(fake.all[ i, j], 0.9)
    }
  }
fake.all = fake.all.corr


### Convert this data to Qnemo file ## WARNING ASSUMES THAT FIRST POP IS WHEAT. Make sure this is correct in Inifile
cat("# Allelic file
###############
[FILE_INFO] {
  col_locus 1
  col_allele 2
  col_mut_freq 3
  col_ini_freq 4
  }\n#locus allele mut_freq ini_freq\n", file = output, append = F)

for(loc in 1:ncol(fake.all)){
  freqs.pres = fake.all[, loc]  
  freqs.abs = 1 - freqs.pres
  
  cat(loc, " 1 ",mu," {", paste(freqs.pres, collapse = " "), "}\n", sep = "", file = output, append = T)
  cat(loc, " 2 ",mu," {", paste(freqs.abs, collapse = " "), "}\n", sep = "", file = output, append = T)
  }