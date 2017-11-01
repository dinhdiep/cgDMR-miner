#!/usr/bin/perl -w
use strict;

my $chr = $ARGV[0];
my $last_state = "NA";
my $last_pos = "NA";
my $start_pos = "NA";
my $cpg_cnt = 0;
my $sum_stats = 0;
my $maxgap = 10000; # 10kbp
my $line = <STDIN>; # skip header
while($line = <STDIN>){
  chomp($line);
  my @f = split "\t", $line;
  my ($pos, $stat, $state) = ($f[0], $f[1], $f[2]);
  if($last_state eq "NA"){
    $start_pos = $pos;
    $cpg_cnt = 1;
    $sum_stats = $stat;
    $last_state = $state;
    $last_pos = $pos;
    next;
  }
  if($last_state ne "NA" and ( $last_state ne $state || $pos - $last_pos > $maxgap ) ){
    # output the current segment and start new
    print $chr, "\t", $start_pos - 1, "\t", $last_pos, "\t", $cpg_cnt, "\t", sprintf("%4.3f", $sum_stats/$cpg_cnt), "\n";
    $start_pos = $pos;
    $cpg_cnt = 1;
    $sum_stats = $stat;
    $last_state = $state;
    $last_pos = $pos;
  }else{
    $cpg_cnt++;
    $sum_stats += $stat;
    $last_state = $state;
    $last_pos = $pos;
  }
}

print $chr, "\t", $start_pos - 1, "\t", $last_pos, "\t", $cpg_cnt, "\t", sprintf("%4.3f", $sum_stats/$cpg_cnt), "\n";

exit 0;
