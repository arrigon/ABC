#!/bin/perl 

use List::Util qw(shuffle);

### Script args
my $input = $ARGV[0];
my $sampling = $ARGV[1];
my $output = "$input.txt";

### IO
open(IN, "$input");
open(SMP, "$sampling");
open(OUT, ">$output");

### parse sampling file
my %samp;
while(my $line = <SMP>){
  chomp($line);

  # split pop Vs Samp effort
  $_ = $line;
  s/[{|}]//g;
  s/^\s+//;
  s/\s+$//;
  $line = $_;

  # store data in hash
  my @fields = split(/\s/, $line);
  $samp{$fields[0]} = $fields[1];
  }


### parse input file (Qnemo results), remove recent migrants, recode to dominant markers.
my $nrline = 1;
my $nloci;
my %store;

while(my $line = <IN>){
  chomp($line);
  my @fields = split(/\s/, $line);

  if($nrline == 1){
    $nloci = $fields[1];
    }
  
  if($nrline > $nloci + 2){
    # get pop and genotype info
    my $pop = $fields[0];
    my $genotype = join("\t", @fields[1..$nloci]);

    # check that focal specimen is not wheat immigrant, skip it if this is the case
    my @parents = @fields[($nloci+3)..($nloci+5)];
    if($pop !~ /01/ && $parents[1] =~ /\_1$/ && $parents[2] =~ /\_1$/){
      next();
      }

    # recode genotypes into dominant markers, WATCH OUT, must fit with coding scheme from IniFreq
    $_ = $genotype;
    s/11/1/g;
    s/12/1/g;
    s/21/1/g;
    s/22/0/g;
    $genotype = $_;

    # store data in a hash
    $store{$pop}{$nrline} = "$pop\t$genotype";
    }

  $nrline++;
  }
close(IN);


### select samples from pops, save to outputs
# Open sampling scheme file
my $popitr = 1;

foreach my $pop (sort(keys(%store))){
  # see how many samples we need from each pop
  my $needed = $samp{$popitr} - 1;
 
  # see what samples are available in pop and sample from there
  my @randkeys;
  @randkeys = shuffle(keys(%{$store{$pop}}));
  @randkeys = @randkeys[0..$needed];
  @randkeys = grep defined, @randkeys;
  foreach my $ind (@randkeys){
    my $line = $store{$pop}{$ind};
    print OUT "$line\n";
    }

  $popitr++;
  }
close(OUT);
