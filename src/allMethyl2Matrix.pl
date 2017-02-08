#!/usr/bin/perl  -w
#Should disable buffering of STDOUT
#Usage:
#allMethyl2Matrix.pl SAMPLES_LIST [minDepth] [minSamples] [sitesList/NA] [cnt/freq/depth]

#this script is more the new style methylFreq files

use strict;
use Statistics::Descriptive;
use POSIX qw(floor);

if(!$ARGV[0]){
	print "Usage:\n";
	print "allMethyl2Matrix.pl [samplesList] [minDepth] [minSamples] [regionsBED] [cnt/freq/depth] [windowSize] [outPrefix]\n";
	print "this script is more the new style methylFreq files, each site should have context either CG, CHG, or CHH\n";
	print "window size is number of bp for a sliding window or dmr for over a region\n";
	print "minSamples is not considered when there is a regionsBED file\n";
	exit 0;
}


$| = 1;

my $list = $ARGV[0];
my $minSamples = 1;
$minSamples = $ARGV[2] if($ARGV[2]);
my $regionsBED = "NA";
$regionsBED = $ARGV[3] if($ARGV[3]);
my $style = "freq";
$style = $ARGV[4];
my $windowSize = 1;
$windowSize = $ARGV[5] if($ARGV[5]);
$windowSize = 0 if($ARGV[5] eq "dmr");
my $outPrefix = $ARGV[6];

open(IN, "$list") || die("Error opening list file.");
my @fileList = <IN>;
close(IN);

my $minDepth=10;
$minDepth = $ARGV[1];
my %allSampleName;
my $print_header =1;
my $binSize = 1000;

my %methylData;
my %methylCount;
my %regionsBins;

if($regionsBED eq "NA"){
	foreach my $file(@fileList){
		chomp($file);
		my @f = split /\t/, $file;
		my ($sampleID, $fileName) = ($f[0], $f[1]);
		next if(!$fileName);
		open(INFILE, "$fileName") || warn("Error in opening file $file.\n");
		while(my $line = <INFILE>){
			next if($line !~ /^chr/);
			chomp($line);
			my @fields = split(/\t/, $line);
			next if($fields[3] == 0);
			my ($chr, $start, $end) = ($fields[0], $fields[1]-1, $fields[1]);
			if($windowSize > 1){
				$fields[1] = floor($fields[1]/$windowSize);
				($start, $end) = ($fields[1]*$windowSize, $fields[1]*$windowSize+$windowSize-1);
			}elsif($windowSize == 0){
				die("Window size set to DMR but no regions file specified\n");
			}
			my $index = "$chr,$start,$end";
			$methylCount{$index}++;
		}
		close(INFILE);
		print "Finish scanning $file.\n";
	}
}else{
	print "Scanning $regionsBED\n";
	open(INFILE, "$regionsBED") || die("Error opening $regionsBED\n");
	while(my $line = <INFILE>){
		chomp($line);
		my @fields  = split "\t", $line;
		my ($chr, $start, $end) = ($fields[0], $fields[1], $fields[2]);
		my $index = "$chr,$start,$end";
		$methylCount{$index} = $minSamples;
		# as long as windowSize is 1 or DMR, make sure to print all of the regions regardless of coverage
		if($end - $start == 1){
			# region is a siteList
			my $i = floor($start/$binSize);
			my $key = "$chr:$i";
			push(@{$regionsBins{$key}}, $index);
		}
		for(my $i = floor($start/$binSize)-3; $i le floor($start/$binSize)+3; $i++){
			my $key = "$chr:$i";
			push(@{$regionsBins{$key}}, $index);
		}
	}
	close(INFILE);
}

foreach my $file (@fileList){
	chomp($file);
	my @f = split /\t/, $file;
	my ($sampleID, $fileName) = ($f[0], $f[1]);
	next if(!$fileName);
	open(INFILE, "$fileName") || die("Error in opening file $fileName.\n");	
	while(my $line = <INFILE>){
		next if($line !~ /^chr/);
		chomp($line);
		my @fields = split(/\t/, $line);
		next if($fields[3] == 0);
		my ($chr, $start, $end) = ($fields[0], $fields[1]-1, $fields[1]);
		my $index = "$chr,$start,$end";
		if($regionsBED ne "NA"){
			my $keepsite = 0;				
			my $bin = floor($fields[1]/$binSize);
			if($regionsBins{$chr.":".$bin}){
				my @candidates = @{$regionsBins{$fields[0].":".$bin}};
				#print scalar(@candidates), "\n";
				foreach my $can (@candidates){
					my ($chr, $start, $end) = split ",", $can;
					next if($fields[0] ne $chr || $fields[1] > $end || $fields[1] < $start);
					if($windowSize == 0){
						$index = $can;
					}elsif($windowSize > 1){
						$start = floor($start/$windowSize);
						($start, $end) = ($start*$windowSize, $start*$windowSize+$windowSize-1);
						$index = "$chr,$start,$end";
					}
					$keepsite = 1;
					last;
				}
			}elsif($methylCount{$index}){
				# check if this is single site list
				$keepsite = 1;
			}
			#print "$index\n" if($keepsite == 1);
			next if($keepsite == 0);
		}else{	
			#skip site if not minimum samples met
			next if(!$methylCount{$index} || $methylCount{$index} < $minSamples);
		}
		my ($c, $t) = (0,0);
		for(my $i = 5; $i < scalar(@fields); $i+=2){
			$c = $fields[$i+1] if($fields[$i] eq 'C');
			$t = $fields[$i+1] if($fields[$i] eq 'T');
		}
		next if( ($c+$t)/$fields[3] < 0.9 );
		$methylData{$index}->{$sampleID}->{'C'}+=$c;
		$methylData{$index}->{$sampleID}->{'T'}+=$t;
	}
	close(INFILE);	
	print "Finish adding $file.\n";
	$allSampleName{$sampleID} = 1;
		
}

my %fileCoverage;
my @samples = keys %allSampleName;

# fillin 0 coverage sites
if($regionsBED ne "NA" && $windowSize <= 1){
	foreach my $index (keys %methylCount){
		next if($methylData{$index});
		$methylData{$index}->{"nullSample"}->{'T'} = 0;
	}
}

my $outfile = "MethylMatrix.$outPrefix.$style";
open(OUT, ">$outfile") || die("Error writing $outfile.\n");
if($print_header){	
	print OUT "1chr_position\t", join("\t", @samples), "\n";
}

foreach my $index(keys(%methylData)){	
		my @m;
		my $cnt =0;
		for(my $i=0; $i<scalar(@samples); $i++){
			if($methylData{$index}->{$samples[$i]}){
				my ($c,$t) = ($methylData{$index}->{$samples[$i]}->{'C'}, $methylData{$index}->{$samples[$i]}->{'T'});
				#$c = 0 if($t == 0);
				if($c+$t < $minDepth){
					push(@m, "NA") if($style ne "cnt");
					push(@m, "0:0") if($style eq "cnt");
					next;
				}
				my $n = $c+$t;
				push(@m, sprintf("%4.3f", $c/($c+$t))) if($style eq "freq");
				push(@m, $c.":".$n) if($style eq "cnt");
				push(@m, $c+$t) if($style eq "depth");
				$cnt++;
				#$fileCoverage{$allFiles[$i]}->{"depth 1"}++ if($c + $t > 1);
				#$fileCoverage{$allFiles[$i]}->{"depth 5"}++ if($c + $t > 5);
				#$fileCoverage{$allFiles[$i]}->{"depth 10"}++ if($c + $t > 10);
			}else{
				push(@m, "NA") if($style ne "cnt");
				push(@m, "0:0") if($style eq "cnt");
			}
		}
		print OUT "$index\t", join("\t", @m), "\n" if($cnt >= $minSamples || $regionsBED ne "NA");
}	
close(OUT);

undef(%methylData);
undef(%fileCoverage);
undef(%methylCount);
undef(%regionsBins);
undef(%allSampleName);
exit 0;