#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Path qw(make_path);
use File::Basename;
use Cwd 'abs_path';

my $script_dir = abs_path(dirname(__FILE__));

my %chrs_mf;
my %chrs_bed;

my $cmd;

sub main{
  my %opts = ();
  getopts('d:i:o:m:n:s:c:p:', \%opts);

  # the default values
  my $min_depth = 10;
  my $mf_list = "NA";
  my $out_dir = "cgDMR_results";
  my $mode = "jsd";
  my $num_cgs = 5;
  my $ifsmooth = "no";
  my $sites = "NA";
  my $pval_cutoff = 0.01;

  $min_depth = $opts{d} if($opts{d});
  $mf_list = $opts{i} if($opts{i});
  $out_dir = $opts{o} if($opts{o});
  $mode = $opts{m} if($opts{m});
  $num_cgs = int($opts{n}) if($opts{n});
  $ifsmooth = $opts{s} if($opts{s});
  $sites = $opts{c} if($opts{c});
  $pval_cutoff = $opts{p} if($opts{p});

  if($mf_list eq "NA"){
    print "\n[cgDMR-miner] Input list is not provided\n";
    printUsage();
  }
  if($ifsmooth =~ m/yes/i and $sites eq "NA" ){
    print "\n[cgDMR-miner] Smoothing is required but CpG bed list file is not provided\n";
    printUsage();
  }

  # generating the output directory
  make_path($out_dir);

  # split the cpg bed file into chromosomes
  if($sites ne "NA"){
    open(SITES_PATHS, "$sites") || die("Error opening $sites\n");
    while(my $line = <SITES_PATHS>){
      chomp($line);
      my @fields = split "\t", $line;
      $chrs_bed{$fields[0]} = $fields[1];
    }
    close(SITES_PATHS);
  }

  # split methylation frequency files into chromosomes
  open(IN, "$mf_list") || die("error reading $mf_list\n");
  while(my $line = <IN>){
    chomp($line);
    my @tmp = split "\t", $line;
    my ($sample, $mf_loc, $chr) = ($tmp[0], $tmp[1], $tmp[2]);
    push(@{$chrs_mf{$chr}}, "$sample\t$mf_loc\t$chr\n");
  }	
  close(IN);

  # change to the output directory
  chdir $out_dir;

  # Processing one chromosome at a time
  foreach my $chr (keys %chrs_mf){
 
    # Write mf_list for current chromosome  
    my @mf_locs = @{$chrs_mf{$chr}};
    open(OUT, ">mf_list_$chr") || die("error writing mf_list_$chr");
    print OUT join("", @mf_locs);
    close(OUT); 

    # Generate methylation summary file
    if($ifsmooth =~ m/no/i){
      $cmd = "perl $script_dir/allMethyl2SlidingWindow.pl mf_list_$chr $num_cgs $min_depth 2 NA $chr";
      system($cmd);
    }else{
     if(generateSmoothedMatrix($chr)){
        $cmd = "perl $script_dir/matrix2summary.pl $chr.matrix 2 NA $chr";
        system($cmd);
     }else{
       print "\n[cgDMR-miner] Unable to finish smoothing $chr. Skipping\n";
       next;
     }
   }
  

    # Segmentation with mhsmm
    $cmd = "Rscript $script_dir/hmm_4states.R MethylSummary.$chr $chr $mode";
    system($cmd);

    # Generate segment bed file from HMM states
    $cmd = "perl $script_dir/sites2segments.pl $chr < $chr.HMM4.sites.txt > $chr.bed";
    system($cmd);
 
    # Generate AMF across the segments
    $cmd = "perl $script_dir/allMethyl2Matrix.pl mf_list_$chr $min_depth 1 $chr.bed cnt dmr $chr";
    system($cmd);

    # Use MethylPy's approach for differential methylation testing within segments
    $cmd = "perl $script_dir/Ctest.pl $min_depth < MethylMatrix.$chr.cnt > $chr.HMM4.p_info";
    system($cmd);

    # Generate methylation level matrix for DMRs
    if(!generateDmrsMatrix("$chr.HMM4.p_info", $pval_cutoff, "MethylMatrix.$chr.cnt")){
      print "[cgDMR-miner] failed to generate weighted average methylation levels for DMRs on $chr. Skipping...\n";
      next;
    }
   
    # cleaning up
    unlink("mf_list_$chr");
    system("gzip $chr.matrix");
  }


}


sub generateDmrsMatrix{

  my $ctest_res = shift;
  my $cutoff = shift;
  my $cnt_matrix_file = shift;
  my $out_matrix_file = $cnt_matrix_file;
  $out_matrix_file =~ s/cnt/levels.txt/g;
  
  my %keep;
  open(CTEST_RES, "$ctest_res") || die("[cgDMR-miner] Error reading $ctest_res\n");
  my $line = <CTEST_RES>;
  while($line = <CTEST_RES>){
    chomp($line);
    my @fields = split "\t", $line;
    if($fields[5] < $cutoff){
      $keep{$fields[0]} = 1;
    }
  }
  close(CTEST_RES);

  open(OUT, ">$out_matrix_file") || die("[cgDMR-miner] Error writing $out_matrix_file\n");
  open(MF_CNT, "$cnt_matrix_file") || die("[cgDMR-miner] Error reading $cnt_matrix_file\n");
  $line = <MF_CNT>;
  chomp($line);
  print OUT $line, "\n";

  while($line = <MF_CNT>){
    chomp($line);
    my @fields = split "\t", $line;
    next if(!$keep{$fields[0]});
    my $cur_methyl = 0;
    my $cur_depth = 0;
    for(my $i = 1; $i < scalar(@fields); $i++){
      my ($methyl, $depth) = split ":", $fields[$i];
      $cur_methyl+=$methyl;
      $cur_depth+=$depth;
    }
    next if($cur_depth == 0);
    print OUT  $fields[0];

    for(my $i = 1; $i < scalar(@fields); $i++){
      my ($methyl, $depth) = split ":", $fields[$i];
      if($depth == 0){
        print OUT "\tNA";
        next;
      }
      print OUT "\t", sprintf("%4.3f", $methyl/$depth);
      next;
    }
    print OUT "\n";  
  }
  close(MF_CNT);

  return 1;

}


sub generateSmoothedMatrix{

  my $chr = shift;
  my $cmd;
  my $directorySmoothed = "MethylSmoothingFiles";
  my $smf_list = "$directorySmoothed/matrices.$chr.list";
  my $cur_list = "$directorySmoothed/cur_mf_list";

  # generating the smoothed output directory
  make_path($directorySmoothed); 
  print "Procssing $chr .. \n";
  if(!$chrs_bed{$chr}){
    print "Error no CpG bed provided for $chr!\n";
    return 0;
  }
  my $cpg_bed = $chrs_bed{$chr};

  open(SMF_LIST_OUT, ">$smf_list") || die("Error writing to $smf_list\n");
  open(MF_LIST_IN, "mf_list_$chr") || die("Error reading mf_list_$chr\n");
  while(my $line = <MF_LIST_IN>){
    chomp($line);
    my ($sample, $mf_file, $mf_chr) = split "\t", $line;
    open(MF_LIST_OUT, ">$cur_list") || die("Error writing to $cur_list\n");
    print MF_LIST_OUT $sample, "\t", $mf_file, "\n";
    close(MF_LIST_OUT);
    $cmd = "$script_dir/goBsmooth_v4.pl $directorySmoothed $sample.$chr $cur_list $cpg_bed 1 $chr";
    system($cmd);
    print SMF_LIST_OUT "$sample\t$directorySmoothed/MethylMatrix.$sample.$chr.smoothed\n"; 
  }
  close(SMF_LIST_OUT);

  $cmd = "$script_dir/table2Matrix.pl $smf_list 3 > $chr.matrix";
  system($cmd);
 
  # cleaning up
  $cmd = "rm -r $directorySmoothed";
  system($cmd);
 
  return 1;

}

sub printUsage{
  print "USAGE: $0 \n";
  print " Note: Rscript must be installed and this script needs the full path to the scripts folder\n";
  print "     -d Default = 10. Minimum read depth per sample for analysis \n";
  print "     -i Required. Input file, a list of methylation frequency files, one line per sample: id<tab>full path to file<tab>chromosome_id \n";
  print "     -o Required. Output directory \n";
  print "     -m Default = jsd. Methylation variation mode, to calculate the chi-square, use \"chsq\", to calculate the Jensen-Shannon divergence, use \"jsd\" \n";
  print "     -n Default = 3. Number of CpGs per sliding window, default is 3, use only odd numbers for best performance \n";
  print "     -s Default = yes [no]. Instead of sliding window, do smoothing with Bsmooth [takes more time], yes[no] \n";
  print "     -c Required for smoothing. A CpG bed files list: chromID<tab>full path to file \n";
  exit 1; 
}

main();

