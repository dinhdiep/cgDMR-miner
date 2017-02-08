#!/usr/bin/perl -w

use strict;
use Cwd qw(abs_path);

my $cur_script = $0;
my @tmp = split /\//, abs_path($cur_script);
pop(@tmp);
my $script_dir = join("\/", @tmp);
my $result_dir = "segmentation_results";
my $temp_dir = "TempFiles";
mkdir $temp_dir;
my $do_remove_temp_dir = 0;
my $samples_info = $ARGV[0];
my $sites = $ARGV[1];
my $minDepth = $ARGV[2];
my $minEffectSize = 0.2;
$minEffectSize = $ARGV[4] if($ARGV[4]);
my $pCutoff = 0.001;
$pCutoff = $ARGV[3] if($ARGV[3]);


my $cmd;

my %samples;
my %chrs;
my %ignoreChrs;

sub main{
        if(!$minDepth || !$sites || !$samples_info){
		printUsage();
	}
	open(SITES_PATHS, "$sites") || die("Error opening $sites\n");
	while(my $line = <SITES_PATHS>){
		chomp($line);
		my @fields = split "\t", $line;
		$chrs{$fields[0]}->{"sites"} = $fields[1];
	}
	close(SITES_PATHS);
	open(SAMPLE_INFO, "$samples_info") || die("Error opening $samples_info\n");
	open(LABEL_GROUPS, ">label_groups") || die("Error writing label_groups\n");
	while(my $line = <SAMPLE_INFO>){
		chomp($line);
		my @fields = split "\t", $line;	
		$samples{$fields[0]}->{"id"} = $fields[1];
		$samples{$fields[0]}->{"mf_list"} = $fields[2];
		print LABEL_GROUPS $fields[0], "\t", $fields[1], "\n";
		open(MF_LIST, "$fields[2]") || die("Error opening $fields[2]\n");
		while(my $mf_line = <MF_LIST>){
			chomp($mf_line);
			my @mf_fields = split "\t", $mf_line;
			if(!$chrs{$mf_fields[2]}->{"sites"}){
				$ignoreChrs{$mf_fields[2]} = 1;
				next;
			}
			push(@{$chrs{$mf_fields[2]}->{$fields[0]}}, $mf_fields[1]);
			
		}
		close(MF_LIST);			
	}
	close(SAMPLE_INFO);
	close(LABEL_GROUPS);
	
	foreach my $chr(keys %ignoreChrs){
		print "[Master] Ignoring $chr because CpG sites list not provided\n";
	}
	##  perform smoothing, by chromosomes, make matrix and run segmentation
	mkdir "MethylMatrixSmoothed";
	my $cmd;
	my %processed_chr;
	foreach my $chr(keys %chrs){
		next if(!$chrs{$chr}->{"sites"});
		print "Processing $chr..\n";
		my $smf_list="$temp_dir/matrices.$chr.list";
		my $amf_list="$temp_dir/amf.$chr.list";
		my $cur_list="$temp_dir/cur_mf_list";
		my $cpg_bed = $chrs{$chr}->{"sites"};
		my $directorySmoothed = "MethylMatrixSmoothed/mat_$chr";
		mkdir $directorySmoothed;
		open(SMF_LIST_OUT, ">$smf_list") || die("Error writing to $smf_list\n");
		open(AMF_LIST_OUT, ">$amf_list") || die("Error writing to $amf_list\n");
		foreach my $sample(keys %samples){
			next if(!$chrs{$chr}->{$sample});
			open(MF_LIST_OUT, ">$cur_list") || die("Error writing to $cur_list\n");
			foreach my $mf_file(@{$chrs{$chr}->{$sample}}){
				print MF_LIST_OUT $sample, "\t", $mf_file, "\n";
				print AMF_LIST_OUT $sample, "\t", $mf_file, "\n";
			}
			close(MF_LIST_OUT);
			$cmd = "$script_dir/goBsmooth_v3.pl $directorySmoothed $sample.$chr $cur_list $cpg_bed 1";
			system($cmd) == 0 or die "[Master] Failed to generate smoothed file (exit $?): $!\nCommand used:\n\t$cmd\n";
			print SMF_LIST_OUT "$sample\t$directorySmoothed/MethylMatrix.$sample.$chr.smoothed\n";
		}
		close(AMF_LIST_OUT);
		close(SMF_LIST_OUT);
		
		$cmd = "$script_dir/table2Matrix.pl $smf_list 3 > $chr.matrix";
		system($cmd) == 0 or warn "[Master] Failed to generate $chr matrix file (exit $?): $!\nCommand used:\n\t$cmd\n";
		$cmd = "Rscript $script_dir/analyzeMatrix.R $chr.matrix 1 label_groups $chr";
		system($cmd) == 0 or warn "[Master] Failed to generate statistics file (exit $?): $!\nCommand used:\n\t$cmd\n";
		$cmd = "Rscript $script_dir/segmentOnly.R $chr.jsd.txt $chr";
		system($cmd) == 0 or warn "[Master] Failed to generate bed file (exit $?): $!\nCommand used:\n\t$cmd\n";
		$cmd = "$script_dir/allMethyl2Matrix.pl $amf_list $minDepth 1 $chr.HMM.seg.txt cnt dmr $chr";
		system($cmd) == 0 or warn "[Master] Failed to generate methylFreq file (exit $?): $!\nCommand used:\n\t$cmd\n";
		$cmd = "$script_dir/Ctest_shuffle.pl $minDepth $pCutoff $minEffectSize < MethylMatrix.$chr.cnt > $chr.HMM.pval";
		system($cmd) == 0 or warn "[Master] Failed to generate statistics file (exit $?): $!\nCommand used:\n\t$cmd\n"; 		
		
		$processed_chr{$chr} = 1;
	}
	system("rm -r $temp_dir") if($do_remove_temp_dir);	
}

sub printUsage{
        print "cgDMR-miner usage\n";
        print "\t ./cgDMR-miner.pl <samples_info> <cpg sites bed list> <minDepth> <p-value cutoff> <minimum effect size>\n";
        print "\t\tsamples_info:  <sample name><tab><sample ID, ie 'case_XX', 'control_XX', 'sample_XX'><tab><path to methylFreq file table>\n";
        print "\t\t\tmethylFreq paths table: <sample name><tab><path to methylFreq><tab><chromosome name>.\n";
        print "\t\t\tnote: for paired samples, indicate the pairing, ie. case1, control1, case2, control2 etc.\n";
        print "\t\tcpg sites paths table: <chr name><tab><bed path> \n\t\t\tbed files containing the CpG positions to be considered, one file per chromosome.\n";
        print "\t\tminDepth: minimum total depth required in each sample for DMR summarization\n";
        print "\t\tp-value cutoff: [0-1] maximum p value for DMRs. Default 0.001 \n\n";
        print "\t\tminimum efect size: [0-1] effect size cutoff for DMRs. Default 0.2 \n\n";
        exit 0;
}

main();
