# ABC pipeline

This repository contains R and Perl scripts developped for Pajkovic et al. 2014. 

http://onlinelibrary.wiley.com/doi/10.1111/mec.12918/abstract
http://datadryad.org/resource/doi:10.5061/dryad.7np24

An updated version of those programs were as well used in 

http://link.springer.com/article/10.1007/s10531-016-1165-z
http://datadryad.org/resource/doi:10.5061/dryad.8f0d4/1

Related to the Pajkovic paper, you can find the following items:
- AFLP/ = empirical AFLP dataset
- data/ = ini files, hard coded
- Qoutputs/ = quantinemo outputs
- Qparams/ = quantinemo parameter files
- Qstats/ = final summary stats
- scripts/ = all needed scripts, look at README.jpg to get an idea how how the pipeline works
- quantinemo = QuantiNemo V 1.6.0 (unreleased), can be downloaded at 
http://www.unil.ch/popgen/softwares/quantinemo/blaser_et_al_evolution.zip
otherwise, contact Samuel Neuenschwander for getting a copy.

Note that the pipeline must be started from scripts/ using the following command:
perl ABC_loader.pl ncores nmax
-ncore = number of parallel instances of quantiNemo you want to run
-nmax = launch jobs until reaching nmax CPU load

Dependencies:
Linux OS, Ubuntu 10.04
Perl v5.10.1
threads;
threads::shared;
File::Basename;
List::Util qw(shuffle);
R >= 3.0.0
e1071
QuantiNemo: V 1.6.0


PS. This project is finished, but the pipeline can be easily adapted to new studies. Feel free to contact me for any collaborations or help. If you use those codes, please cite either reference.# ABC
