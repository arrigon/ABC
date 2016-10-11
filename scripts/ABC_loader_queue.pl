#!/bin/perl

# Usage:
# perl ABC_loader_queue.pl NCPU MaxUse
#
# NCPU = Number of cores to be used by each process
# MaxUse = Maximal load allowed on the computing server.
#
#

use threads;
use threads::shared;
 

use File::Basename;

my $rootdir = getcwd;

#### Get script arguments
$CPU = $ARGV[0];
$MaxUse = $ARGV[1];
$scriptname = "ABC_loader";

@listruns = (1..$CPU);


my @threads;
foreach $run (@listruns){
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0] + 1;

  print "### $scriptname : Current CPU load is $load\n\n";
  
  if($load < $MaxUse) {
    #### Automated loader
    print "### $scriptname : RUNNING run #$run\n";

    my $t = threads->new(\&rungenome, $run);
    push(@threads,$t);
    sleep(10);

    } else {
    print "### $scriptname : run $run is waiting... CPU load is maximised\n\n";
    sleep(10);
    redo;
    }
  }

foreach (@threads) {
  my $run = $_->join;
  print "done with $run\n";
  }


sub rungenome {
  ### Get CPU load
  my ($run) = ($_[0]);

  # Launch ABC pipeline
  $command = "R CMD BATCH \"--args run=$run\" RUN_ABC.r";
  print "### $scriptname : $command\n";
  system("$command");
  }
