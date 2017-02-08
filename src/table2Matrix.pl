#!/usr/bin/perl -w
use strict;

# make a matrix from specified column from all the files in a list

my $list_files = $ARGV[0];
my $column_id = $ARGV[1];

if(!$column_id){
	print "Missing input files\n";
	print "./makeMatrix.pl <files_list> <column_id>\n";
	print "list is a two columns table: column 1 - sample_id, column 2 - path to file\n";
	print "file should be format: <chr><tab><position><tab><value1><tab><value2><tab><valueN>\n";
}

my %table;
my @samples;

open(IN, "$list_files") || die("Error reading $list_files\n");
while(my $line = <IN>){
	chomp($line);
	my @tmp = split "\t", $line;
	push(@samples, $tmp[0]);
	open(DATA, "$tmp[1]") || die("Error opening $tmp[1]\n");
	while(my $data_in = <DATA>){
		chomp($data_in);
		my @fields = split "\t", $data_in;
		next if($fields[0] eq "chr");
		$table{$fields[0]."\t".$fields[1]}->{$tmp[0]} = $fields[$column_id-1];
	}
}
print "chr\tstart\t", join("\t", @samples), "\n";
foreach my $value (keys %table){
	print $value;
	for(my $i = 0; $i < scalar(@samples); $i++){
		if($table{$value}->{$samples[$i]}){
			print "\t", $table{$value}->{$samples[$i]};
		}else{
			print "\tNA";
		}
	}
	print "\n";
}
