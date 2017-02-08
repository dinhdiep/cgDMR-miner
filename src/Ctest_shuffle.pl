#!/usr/bin/perl -w

use strict;
use Math::Random qw(random_binomial random_poisson);

my %cpgTable;
my %samplesID;
my $minDepth = $ARGV[0]; 
my $pCutoff = $ARGV[1]; 
my $minEffectSize = $ARGV[2]; 
my $num_perm = 3000;
my $num_r = 36;

my $line = <STDIN>;
chomp($line);
my @header = split "\t", $line;
for(my $i = 1; $i < scalar(@header); $i++){
	$samplesID{$header[$i]}= $i;
}

print "index\ttotal_cpg_coverage\taverage_methyl\ttotal_methyl_freq\tchsq_statistisc_df1\tpermutation_p_value\tmaximum_difference\tdmr_status\n";
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
	next if($num_samples_mindepth < 2);
	my $chsq = 0;
	#my $p = $o_m/$o_d;
	my $p = $sum_ave/$num_samples_mindepth; # use unweighted average to get estimation of binomial p 
	if($p == 0 || $p == 1){
		print $value, "\t", $o_d, "\t", $p, "\t", $o_m/$o_d, "\t", 0, "\t", 1, "\t", 0, "\n";
		next;
	}
	for(my $i = 0; $i < $num_samples_mindepth; $i++){	
		my $observed_count = $m_values[$i];
		my $expected_count = $n_values[$i]*$p;
		$chsq += ($observed_count  - $expected_count)*($observed_count - $expected_count);
	}
	$chsq = $chsq / ($p*(1-$p));
	my $chsq_final = ($chsq - $o_d)*($chsq - $o_d)/(2*$denom);
	my $cur_perm = 0;
	my $cur_R = 0;
	while($cur_perm < $num_perm){
		my $k_sum_ave = 0;
		my $k_chsq = 0;
		my @k_boot;
		for(my $i = 0; $i < $num_samples_mindepth; $i++){
			my $cur_k = random_binomial(1, $n_values[$i], $p);
			push(@k_boot,$cur_k);
			$k_sum_ave+= $cur_k/$n_values[$i];
		}
		my $k_p = $k_sum_ave / $num_samples_mindepth;
		if($k_p == 0 || $k_p == 1){
			$cur_perm++;
			next;
		}
		for(my $i = 0; $i < $num_samples_mindepth; $i++){
			my $observed_count = $k_boot[$i];
			my $expected_count = $n_values[$i]*$k_p;
			$k_chsq += ($observed_count  - $expected_count)*($observed_count - $expected_count);
		}
		$k_chsq = $k_chsq / ($k_p*(1-$k_p));
		$cur_perm++;
		$cur_R++ if($k_chsq >= $chsq);
		last if($cur_R == $num_r);
	}
	my $est_p = 1;
	if($cur_R == $num_r){
		$est_p = $num_r/$cur_perm;
	}else{
		$est_p = ($cur_R+1)/($num_perm+1);
	}
	my $dmr_status = "NOPASS";
	if($est_p < $pCutoff and sqrt($chsq_final/$sum_depth) > $minEffectSize){
		$dmr_status = "PASS";
	}
	print $value, "\t", $sum_depth, "\t", $p, "\t", $o_m/$o_d, "\t", $chsq_final, "\t", $est_p, "\t", $max_mf-$min_mf, "\t", $dmr_status, "\n";
}
