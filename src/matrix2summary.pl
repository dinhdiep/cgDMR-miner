#!/usr/bin/perl  -w
#Should disable buffering of STDOUT

use strict;
use Statistics::Descriptive;
use POSIX qw(floor);

if(!$ARGV[0]){
  print "Usage:\n";
  print "$0 [matrix file] [minSamples] [chr:start-end] [outPrefix]\n";
  exit 0;
}

$| = 1;


### Functions ###

sub log2{
  my $n = shift;
  return log($n)/log(2);
}

sub methyl_transform{
  my $x = shift;
  my $delta = 0.001;
  return ( $delta + $x * (1 - 2*$delta) );
}

sub entropy_log2{
  my @x = @_;
  my $entropy = 0;
  foreach my $value(@x) {
    $entropy += $value * log2($value);
  }
  return ( -1 * $entropy );
}

### END OF FUNCTIONS ###

my $matrixFile = $ARGV[0];

my $minSamples = 1;
$minSamples = $ARGV[1] if($ARGV[1]);
print "Number of samples for including window: $minSamples\n";

my ($chr, $start, $end) = ("NA", 0, 3e9);
($chr, $start, $end) = split /[:-]/, $ARGV[2] if($ARGV[2] and $ARGV[2] ne "NA");
print "Region subset: $chr,$start,$end\n";

my $outPrefix = $ARGV[3];
my $outfile = "MethylSummary.$outPrefix";
print "Output file name: $outfile\n";

my %allSampleName;
my $print_header = 1;

my %methylData;

open(IN, "$matrixFile") || die("Error opening $matrixFile\n");
my $line = <IN>;
my @header = split /[\t\ ]/, $line;
for(my $i = 2; $i < scalar(@header); $i++){
  $allSampleName{$header[$i]} = 1;
}
my @samples = keys %allSampleName;

open(OUT, ">$outfile") || die("Error writing $outfile.\n");
if($print_header){  
  print OUT "first_pos\tmiddle_pos\tlast_pos\tchsq\tjsd\ttotal_depth\ttotal_ave\n";
}

while($line = <IN>){
  chomp($line);
  my @fields = split /[\t\ ]/, $line;
  my ($qchr, $qpos) = ($fields[0], $fields[1]);
  next if($qchr ne $chr and $chr ne "NA");
  next if($qpos > $end || $qpos < $start);
  for(my $i = 2; $i < scalar(@fields); $i++){
    $methylData{$qpos}->{$header[$i]} = $fields[$i];
  }
}

close(IN);

my @sorted_pos = sort {$a <=> $b} keys %methylData;
my $numSites = scalar(@sorted_pos);
print "Total number of sites to consider: $numSites\n";

for( my $i = 0; $i < $numSites; $i++){

  my $pos = $sorted_pos[$i];

  #ave/chsq/jsd/depth/cnt
  my ($sum_m, $sum_d, $chsq) = (0, 0, 0); # for ave, cnt, depth, chsq
  my (@x , @y, @p1, @p2, @mid1, @mid2, @uni); # for jsd
  my (@m, @d); # for chsq
  my ($sum_x, $sum_y) = (0, 0); # for jsd

  foreach my $sId (@samples){

    next if($methylData{$pos}->{$sId} eq "NA");
    my $ave_x = $methylData{$pos}->{$sId};

    $sum_m += $ave_x;
    push(@m, $ave_x);
    
    my $transformed_x = methyl_transform($ave_x);
    push(@x, $transformed_x);
    push(@y, 1 - $transformed_x);
    $sum_x += $transformed_x;
    $sum_y += ( 1 - $transformed_x );
  }

  next if(scalar(@m) < $minSamples || $sum_m == 0);
  
  for( my $k = 0; $k < scalar(@m); $k++){
    my $observed_mf = $m[$k];
    my $expected_mf = $sum_m/scalar(@m);
    $chsq += ($observed_mf  - $expected_mf)*($observed_mf - $expected_mf)/$expected_mf;    

    my $prob_x = $x[$k]/$sum_x;
    my $prob_y = $y[$k]/$sum_y;

    push(@mid1, (0.5*$prob_x + 0.5/scalar(@m)));
    push(@mid2, (0.5*$prob_y + 0.5/scalar(@m)));
    push(@p1, $prob_x);
    push(@p2, $prob_y);
    push(@uni, 1/scalar(@m));
  }

  my $jsd1 = entropy_log2(@mid1) - 0.5*entropy_log2(@p1) - 0.5*entropy_log2(@uni);
  my $jsd2 = entropy_log2(@mid2) - 0.5*entropy_log2(@p2) - 0.5*entropy_log2(@uni);
  my $jsd = $jsd1;
  $jsd = $jsd2 if($jsd2 < $jsd1);
 
  print OUT "$pos\t$pos\t$pos\t$chsq\t$jsd\t$sum_m\t", scalar(@m), "\n";
  
}  
close(OUT);

undef(%methylData);
undef(%allSampleName);
exit 0;
