#!/usr/bin/perl -w

use strict;

my %table;
my %samples;

while(my $line = <STDIN>){
	chomp($line);
	open(IN, "$line") || die("Error opening $line\n");
	my $mline = <IN>;
	chomp($mline);
	my @header = split "[\t\ ]", $mline;
	for(my $i = 1; $i < scalar(@header); $i++){
		$samples{$header[$i]} = 1;
	}
	while(my $mline = <IN>){
		chomp($mline);
		my @fields = split "[\t\ ]", $mline;
		for(my $i = 1; $i< scalar(@fields); $i++){
			$table{$fields[0]}->{$header[$i]} = $fields[$i];
		}
	}
	close(IN);
}
my @all_samples = sort keys %samples;
print "chr_position\t", join("\t", @all_samples), "\n";
foreach my $index(keys %table){
	print $index;
	for(my $i = 0; $i < scalar(@all_samples); $i++){
		print "\t", $table{$index}->{$all_samples[$i]};
	}
	print "\n";
}
