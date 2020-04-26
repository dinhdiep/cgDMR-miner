#!/usr/bin/perl -w

use strict;

#use Statistics::Distributions qw(chisqrprob);
my %cpgTable;
my %samplesID;
my $minDepth = $ARGV[0]; #15

my $line = <STDIN>;
chomp($line);
my @header = split "\t", $line;
for(my $i = 1; $i < scalar(@header); $i++){
  $samplesID{$header[$i]}= $i;
}


print "index\ttotal_depth\taverage_methylation\tweighted_average_methylation\tabsolute_difference\n";

while($line = <STDIN>){
  chomp($line);
  my @f = split "\t", $line;
  my $value = $f[0];
  my ($o_m, $o_d) = (0,0);
  my $sum_ave = 0;
  my $sum_depth = 0;
  my $num_samples_mindepth=0;
  my @n_values;
  my @m_values;
  my $max_mf = 0;
  my $min_mf = 1;
  my $denom = 0;
  foreach my $id (keys %samplesID){
    my $i = $samplesID{$id};
    my ($m, $d) = split ":", $f[$i];
    next if($d < $minDepth); # ignore samples with low coverage
    $sum_depth += $d;
    $o_m+=$m;
    $o_d+=$d;
    push(@m_values, $m);
    push(@n_values, $d);
    if($m/$d > $max_mf){
      $max_mf = $m/$d;
    }
    if($m/$d < $min_mf){
      $min_mf = $m/$d;
    }
    $sum_ave+= $m/$d;
    $num_samples_mindepth++;
    $denom += $d * ($d - 1);
  }
  #my $p = $o_m/$o_d;
  my $p = $sum_ave/$num_samples_mindepth; # use unweighted average to get estimation of binomial p 
  print $value, "\t", $sum_depth, "\t", $p, "\t", $o_m/$o_d, "\t", $max_mf-$min_mf, "\n";
}
