#!/usr/bin/perl  -w
#Should disable buffering of STDOUT
#Usage:
#mf2slidingWindow SAMPLES_LIST [numCpGs] [minDepth] [minSamples] [chsq/jsd/depth/cnt/freq] [outPrefix]

use strict;
use POSIX qw(floor);

if(!$ARGV[0]){
  print "Usage:\n";
  print "$0 SAMPLES_LIST [numCpGs] [minDepth] [minSamples] [chr:start-end] [outPrefix]\n";
  print "this script is for the new style methylFreq files, each site should have context either CG, CHG, or CHH\n";
  print "sliding window across adjacent CpGs, not considering base pair distances\n";
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

my $list = $ARGV[0];

my $numCpGs = 1;
$numCpGs = $ARGV[1] if($ARGV[1]);
print "Number of CpGs per sliding window: $numCpGs\n";

my $minDepth = 10;
$minDepth = $ARGV[2] if($ARGV[2]);
print "Depth threshold for including sample: $minDepth\n";

my $minSamples = 1;
$minSamples = $ARGV[3] if($ARGV[3]);
print "Number of samples for including window: $minSamples\n";

my ($chr, $start, $end) = ("NA", 0, 3e9);
($chr, $start, $end) = split /[:-]/, $ARGV[4] if($ARGV[4]);
#print "Region subset: $chr,$start,$end\n";

my $outPrefix = $ARGV[5];
my $outfile = "MethylSummary.$outPrefix";
print "Output file name: $outfile\n";

open(IN, "$list") || die("Error opening list file.");
my @fileList = <IN>;
close(IN);

my %allSampleName;
my $print_header = 1;

my %methylData;

foreach my $file (@fileList){
  chomp($file);
  my @f = split /\t/, $file;
  my ($sampleID, $fileName) = ($f[0], $f[1]);
  next if(!$fileName);
  open(INFILE, "$fileName") || die("Error in opening file $fileName.\n");  
  while(my $line = <INFILE>){
    #next if($line !~ /^chr/);
    chomp($line);
    my @fields = split(/\t/, $line);

    next if($fields[3] == 0);
    next if($chr ne "NA" and $chr ne $fields[0]); # Cpgs subset
    next if($chr ne "NA" and $fields[1] > $end || $fields[1] < $start); # Cpgs subset

    my ($c, $t) = (0,0);
    for(my $i = 5; $i < scalar(@fields); $i+=2){
      $c = $fields[$i+1] if($fields[$i] eq 'C');
      $t = $fields[$i+1] if($fields[$i] eq 'T');
    }
    next if( ($c+$t)/$fields[3] < 0.9 );
    $methylData{$fields[1]}->{$sampleID}->{'C'}+=$c;
    $methylData{$fields[1]}->{$sampleID}->{'T'}+=$t;
  }
  close(INFILE);  
  print "Finish adding $file.\n";
  $allSampleName{$sampleID} = 1;
    
}

my @samples = keys %allSampleName;

open(OUT, ">$outfile") || die("Error writing $outfile.\n");
if($print_header){  
  print OUT "first_pos\tmiddle_pos\tlast_pos\tchsq\tjsd\ttotal_depth\ttotal_ave\n";
}

my @sorted_pos = sort {$a <=> $b} keys %methylData;
my $numSites = scalar(@sorted_pos);
print "Total number of sites to consider: $numSites\n";

for( my $i = 0; $i < $numSites - $numCpGs + 1; $i++){

  #ave/chsq/jsd/depth/cnt
  my ($sum_m, $sum_d, $chsq) = (0, 0, 0); # for ave, cnt, depth, chsq
  my (@x , @y, @p1, @p2, @mid1, @mid2, @uni); # for jsd
  my (@m, @d); # for chsq
  my ($sum_x, $sum_y) = (0, 0); # for jsd

  foreach my $sId (@samples){

    my ($t_m, $t_d) = (0, 0);

    for( my $j = 0; $j < $numCpGs; $j++){
      #last if($i + $j >= scalar(@sorted_pos));
      my $pos = $sorted_pos[$i+$j];
      if($methylData{$pos}->{$sId}->{'C'}){
        $t_m += $methylData{$pos}->{$sId}->{'C'};
        $t_d += $methylData{$pos}->{$sId}->{'C'};
      }
      if($methylData{$pos}->{$sId}->{'T'}){
        $t_d += $methylData{$pos}->{$sId}->{'T'};
      }
      
    }
   
    next if($t_d < $minDepth);

    $sum_m += $t_m;
    $sum_d += $t_d;
    push(@m, $t_m);
    push(@d, $t_d);
    
    my $ave_x = $t_m/$t_d;
    my $transformed_x = methyl_transform($ave_x);
    push(@x, $transformed_x);
    push(@y, 1 - $transformed_x);
    $sum_x += $transformed_x;
    $sum_y += ( 1 - $transformed_x );
  }

  next if(scalar(@m) < $minSamples || $sum_m == 0);
  
  for( my $k = 0; $k < scalar(@m); $k++){
    my $observed_mf = $m[$k];
    my $expected_mf = $d[$k] * ($sum_m/$sum_d);
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
  
  print OUT "$sorted_pos[$i]\t$sorted_pos[$i + int($numCpGs/2)]\t$sorted_pos[$i+$numCpGs-1]\t$chsq\t$jsd\t$sum_d\t", $sum_m/$sum_d, "\n";
  
}  
close(OUT);

undef(%methylData);
undef(%allSampleName);
exit 0;
