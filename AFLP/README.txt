This folder contains the following items:
- TGCTAATCBinary.txt = complete binary AFLP matrix, includes extra specimens, blank controls and replicates. 
		       Failed specimens have already been removed from this dataset.
- TGCTAATCBinary.info = corresponding metadata (including PCR plate informations, replicate IDs, outgroups, etc).
		        but look at FocalSpecimens.txt to get the focal set of specimens on which is based the present study

- FocalSpecimens.txt = Informations for specimens analysed in the present study:
    Tag - specimen ID
    Pop - population ID
    Comment - includes wheat cultivar name and accession number in CHangin's seed bank
    TypAeg - species status
    NDD, EDD - Geographical coordinates, in decimal degrees
    Alt - Altitude
    Taes,Tdur,Aeg - Indicative SRTUCTURE admixture values. Not used in the MS, we instead rely on InSctruct.
    Mil - Level of remoteness to cultivated fields. Not used in the MS, we rely on actual geographic distances
    Focused on specimens used in the present study.

- ComputeInitialFrequencies_variant1.r = script to compute the initial AFLP band frequencies -> produce data/Ini_real_phen_freqs.txt
- SummaryStats_empirical_mean_sd_v3.r = script to compute summary statistics from the empirical AFLP dataset
- SummaryStats_empirical_mean_sd_100loci_pruning0.15_R1.25.txt = empirical summary statistics, as used for parameter estimation and model selection.
