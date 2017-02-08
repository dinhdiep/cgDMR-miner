#!/usr/bin/perl -w

use strict;

my $line = <STDIN>;
chomp($line);
print $line, "\n";

while(my $line = <STDIN>){
	chomp($line);
	my @fields = split "\t", $line;
	my $cur_methyl = 0;
	my $cur_depth = 0;
	for(my $i = 1; $i < scalar(@fields); $i++){
		my ($methyl, $depth) = split ":", $fields[$i];
		$cur_methyl+=$methyl;
		$cur_depth+=$depth;
	}
	next if($cur_depth == 0);
	print $fields[0];
	my $average = int(1000*($cur_methyl/$cur_depth));
	for(my $i = 1; $i < scalar(@fields); $i++){
		my ($methyl, $depth) = split ":", $fields[$i];
		if($depth == 0){
			print "\tNA";
			next;
		}
		if($ARGV[0] eq "convertOnly"){
			print "\t", sprintf("%4.3f", $methyl/$depth);
			next;
		}
		$depth = 1000 if($depth > 1000);
		$methyl = 0;
		for(my $j = 0; $j < $depth; $j++){
			my $random = 1 + int(rand(1000));
			$methyl++ if($random <= $average);
		}
		print "\t", sprintf("%4.3f", $methyl/$depth);
	}
	print "\n";
}
