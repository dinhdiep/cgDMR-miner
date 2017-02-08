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
my $minDepth = $ARGV[4];
$minDepth = 1 if(!$minDepth);
my $nsites=14;
my $hsize=500;
my $maxgap=50000;
my $cmd;

sub main{
	die("MethylFreq file not provided!\n") if(!$mf_list);
	die("Sample name not provided!\n") if(!$cur_sample);
	my $start_time = time();
	$cmd = "$script_dir/allMethyl2Matrix.pl $mf_list 1 1 $sites cnt 1 $cur_sample";
	my $mycnt = "MethylMatrix.$cur_sample.cnt";
	system($cmd) == 0 or die "[methylFreq2SmoothedValues] Failed to process methylfreq files (exit $?): $!\nCommand used:\n\t$cmd\n";
	$cmd = "sed 's/[,:]/\\t/g' $mycnt | grep -v position | sort -k2,2n > $mycnt.sorted";
	system($cmd) == 0 or die "[methylFreq2SmoothedValues] Failed to generate test cnt file (exit $?): $!\nCommand used:\n\t$cmd\n";

	my $out_name = $mycnt;
	$out_name =~ s/cnt/smoothed/g;
	my $sub_i = 1;
	my $test_i = 0;
	my $cur_i_size = 0;
	my %split_matrix;
	my $last_pos = 0;

	open(SORTED, "$mycnt.sorted") || die("Error reading $mycnt.sorted\n");
	while(my $line = <SORTED>){
		my @f = split "\t", $line;
		if($last_pos > 0 and $f[1] - $last_pos > $maxgap){
			$sub_i++;
			$cur_i_size=0;
			$last_pos = $f[1];
		}else{
			$last_pos = $f[1];
		}
		push(@{$split_matrix{"sub.$sub_i"}}, $line);
		$cur_i_size++;
		if($test_i eq 0 and $cur_i_size > 100000){
			$test_i = $sub_i;
		}
	}
	close(SORTED);
	
	# cross validation of smoothing parameters	
	open(SUBSORTED, ">$mycnt.tmp") || die("Error writing tmp file\n");
	print SUBSORTED @{$split_matrix{"sub.$test_i"}};
	close(SUBSORTED);
        $cmd = "$script_dir/bsmooth_test_v2.pl $mycnt.tmp $out_dir";
        system($cmd) == 0 or die "[methylFreq2SmoothedValues] Failed to generate test data (exit $?): $!\nCommand used:\n\t$cmd\n";
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
		$cmd = "cat $out_dir/$out_name.$i.tmp | grep -v V2 >> $out_dir/$out_name";
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
	system($cmd) == 0 or die "[ methylFreq2SmoothedValues] Failed to generate test ns$i h$j smoothed methylation file (exit $?): $!\nCommand used:\n\t$cmd\n";	
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
main();
