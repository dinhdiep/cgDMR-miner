#!/usr/bin/perl -w

use strict;
use strict;
use Statistics::LSNoHistory;
use Statistics::Basic qw(median);
use Cwd qw(abs_path);

my $cur_script = $0;
my @tmp = split /\//, abs_path($cur_script);
pop(@tmp);
my $script_dir = join("\/", @tmp);

my $out_dir = $ARGV[0];
my $cur_sample = $ARGV[1];
my $mf_list = $ARGV[2];
my $sites = $ARGV[3]; 
my $minDepth = 1; 
$minDepth = $ARGV[4] if($ARGV[4]);
my $cur_chrom = "NA";
$cur_chrom = $ARGV[5] if($ARGV[5]);
my $nsites=14;
my $hsize=500;
my $maxgap=50000;
my $cmd;

my %methylFreq;

sub main{
  
  my $mycnt = "MethylMatrix.$cur_sample.cnt";
  
  my $start_time = time();
  
  # Reads the sites file
  open(SITES, "$sites") || die("[methylFreq2SmoothedValues] Sites file not provided!\n");
  while(my $line = <SITES>){
    chomp($line);
    my @tmp = split "\t", $line;
    $methylFreq{$tmp[1]}->{"C"} = 0;
    $methylFreq{$tmp[1]}->{"T"} = 0;
  }
  close(SITES);

  # Reads the methylation frequency file
  open(LIST, "$mf_list") || die("[methylFreq2SmoothedValues] $mf_list is not readable\n");
  while(my $line = <LIST>){
    chomp($line);
    my @tmp = split "\t", $line;
    open(MF, "$tmp[1]") ||  die("[methylFreq2SmoothedValues] $tmp[1] is not readable\n");
    while(my $mf_line = <MF>){
       chomp($mf_line);
       my @fields = split "\t", $mf_line;
       my $pos = $fields[1] - 1;
       next if(!$methylFreq{$pos});
       my ($c, $t) = (0,0);
       for(my $i = 5; $i < scalar(@fields); $i+=2){
         $c = $fields[$i+1] if($fields[$i] eq 'C');
         $t = $fields[$i+1] if($fields[$i] eq 'T');
       }
       if( ($c+$t)/$fields[3] < 0.9 ){
         $c = 0;
         $t = 0;
       }
       next if(!$methylFreq{$pos});
       $methylFreq{$pos}->{"C"} += $c;
       $methylFreq{$pos}->{"T"} += $t;
    }
    close(MF);
  }
  close(LIST); 

  my $out_name = $mycnt;
  $out_name =~ s/cnt/smoothed/g;


  # Split chromosome into subsets of regions
  my $sub_i = 1;
  my $test_i = 0;
  my $cur_i_size = 0;
  my %split_matrix;
  my $last_pos = 0;

  my @sorted_pos = sort {$a <=> $b} keys %methylFreq;
  for(my $i = 0; $i < scalar(@sorted_pos); $i++){
    my $pos = $sorted_pos[$i];
    my $methyl = $methylFreq{$pos}->{"C"};
    my $depth = $methylFreq{$pos}->{"C"} + $methylFreq{$pos}->{"T"};
    my $out_line = "$cur_chrom\t$pos\t$pos\t$methyl\t$depth\n";

    if($last_pos > 0 and $pos - $last_pos > $maxgap){
      $sub_i++;
      $cur_i_size=0;
      $last_pos = $pos;
    }else{
      $last_pos = $pos;
    }
    push(@{$split_matrix{"sub.$sub_i"}}, $out_line);
    $cur_i_size++;
    if($test_i eq 0 and $cur_i_size > 100000){
      $test_i = $sub_i;
    }
  }
  
  # Cross validation of smoothing parameters  
  open(SUBSORTED, ">$mycnt.tmp") || die("Error writing tmp file\n");
  print SUBSORTED @{$split_matrix{"sub.$test_i"}};
  close(SUBSORTED);

  # Generate test sets
  if(!generateTestSets("$mycnt.tmp", $out_dir)){
     die "[methylFreq2SmoothedValues] Failed to generate test data\n";
  }

  my $best_cor = 0;
  for(my $i = 14; $i <= 50; $i+=4){
    my $cur_cor = ValidateNsHsize($i, $hsize, "$out_dir/test_set");
    #print "$i\t$hsize\t$cur_cor\n";
    if($cur_cor > $best_cor){
      $nsites = $i;
      $best_cor = $cur_cor;
    }
  }

  $best_cor = 0;

  for(my $j = 500; $j<=2000; $j+=200){
    my $cur_cor = ValidateNsHsize($nsites, $j, "$out_dir/test_set");
    #print "$nsites\t$j\t$cur_cor\n";
    if($cur_cor > $best_cor){
      $hsize = $j;
      $best_cor = $cur_cor;
    }
  }

  print "Choosen smoothing params: ns $nsites h $hsize for $cur_sample with $best_cor\n";

  # perform smoothing
  
  unlink("$out_dir/$out_name");
  for(my $i = 1; $i <= $sub_i; $i++){
    open(SUBSORTED, ">$mycnt.tmp") || die("Error writing tmp file\n");
    print SUBSORTED @{$split_matrix{"sub.$i"}};
    close(SUBSORTED);
    $cmd = "Rscript $script_dir/runSmoothing_latest.R $mycnt.tmp $nsites $hsize $out_dir/$out_name.$i.tmp";
    #print $cmd, "\n";
    system($cmd) == 0 or die "[methylFreq2SmoothedValues] Failed to run smoothing with bsseq files (exit $?): $!\nCommand used:\n\t$cmd\n";
    $cmd = "cat $out_dir/$out_name.$i.tmp | grep -v V2  >> $out_dir/$out_name";
    system($cmd) == 0 or die "[methylFreq2SmoothedValues] Failed to merge bsseq files (exit $?): $!\nCommand used:\n\t$cmd\n";
    unlink("$out_dir/$out_name.$i.tmp");
    unlink("$mycnt.tmp");
    @{$split_matrix{"sub.$i"}} = ();
    undef(@{$split_matrix{"sub.$i"}});
  }
  my $end_time = time();
  print "Total time to process $cur_sample: ", $end_time - $start_time, "\n";
  unlink($mycnt);
  unlink("$mycnt.sorted");
  unlink("$out_dir/cross_validate_out");
  unlink("$out_dir/test_set");
}

sub ValidateNsHsize{
  my $i = shift;
  my $j = shift;
  my $validate_set = shift;
  my $c = "ns$i"."_"."h$j";
  $cmd = "Rscript $script_dir/runSmoothing_latest.R $out_dir/cross_validate_out $i $j $out_dir/$c 2>&1 $out_dir/Rout.txt\n";
  system($cmd) == 0 or die "[methylFreq2SmoothedValues] Failed to generate test ns$i h$j smoothed methylation file (exit $?): $!\nCommand used:\n\t$cmd\n";  
  my $corrStat = Statistics::LSNoHistory->new;
  my %table;
  my @diff;
  open(IN1, "$out_dir/$c") || die("error reading $c\n");
  open(IN2, "$validate_set") || die("error reading $validate_set\n");
  while(my $line = <IN2>){
    chomp($line);
    my @tmp = split "\t", $line;
    my ($c, $depth) = ($tmp[2], $tmp[3]);
    next if($depth == 0);
    $tmp[1]++;
    $table{$tmp[0]}->{$tmp[1]} = sprintf("%4.3f", $c/$depth);
  }
  close(IN2);

  my $header = <IN1>;
  my $N=0;
  my $DiffSum=0;
  while(my $line = <IN1>){
    chomp($line);
    my @tmp = split "\t", $line;
    for(my $i = 2; $i < scalar(@tmp); $i++){
      next if(!defined($table{$tmp[0].":".$tmp[1]}->{$i}));
      my $m_validate = $table{$tmp[0].":".$tmp[1]}->{$i};
      my $m_Smoothed = sprintf("%4.3f", $tmp[$i]);
      $corrStat->append_point($m_Smoothed, $m_validate);
      push(@diff, abs($m_Smoothed - $m_validate));
      $DiffSum+=abs($m_Smoothed - $m_validate);
      $N++;
    }
  }
  close(IN1);
  unlink("$out_dir/$c");
  return($corrStat->pearson_r);
}

sub generateTestSets {

  # sampling rate is 10%

  my $matrix_in = shift;
  my $out_dir = shift;

  my $matrix_out = "$out_dir/cross_validate_out";
  my $list_out = "$out_dir/test_set";
  my $shell_out = "$out_dir/test_smoothing.sh";

  open(IN, "$matrix_in") || die("error opening $matrix_in\n");
  open(OUT1, ">$matrix_out") || die("error writing $matrix_out\n");
  open(OUT2, ">$list_out") || die("error writing $list_out\n");

  my $num_sites = 0;
  while(defined(my $line = <IN>) and $num_sites < 100000){
    chomp($line);
    my @tmp = split "\t", $line;
    print OUT1 $tmp[0], "\t", $tmp[1], "\t", $tmp[2];
    my $throwOut = int(rand(10)); #integer from 0-9, so 0 for 1/10 chance
    if($throwOut <= 0.5){
      print OUT1 "\t0\t0\n";
      print OUT2 $tmp[0], ":", $tmp[1], "\t1\t", $tmp[3], "\t", $tmp[4], "\n";
    }else{
      print OUT1 "\t", $tmp[3], "\t", $tmp[4], "\n";
    }
    $num_sites++;
  }
  close(IN);
  close(OUT1);
  close(OUT2);

  return 1;

}
main();
