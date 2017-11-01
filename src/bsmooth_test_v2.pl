#!/usr/bin/perl -w

# sampling rate is 10%
my $matrix_in = $ARGV[0];
my $out_dir = $ARGV[1];

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

exit 0;
