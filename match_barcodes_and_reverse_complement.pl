#!/usr/bin/perl
use strict;
use warnings;
use Bio::Perl;
use Bio::SeqIO;

# This script stores the barcode from each header in file1.  It then searches
# file2 for a sequence with a matching header.  If it finds one, it appends
# the reverse complement of the stored barcode to the header of file2 and 
# then outputs the reverse complement of the sequence in file2

my ($file1, $file2) = @ARGV;
my $output;

my @temp = split (/\//, $file2);
if (scalar (@temp)){
	if ($temp[-1] =~ m/(.*)(.fasta)/){
		my ($file2_prefix) = split (/\./, $file2);
		$output = $file2_prefix. '.cims.rc.fasta'; 
	}
}
else{
	if ($file2 =~ m/(.*)(.fasta)/){
		$output = $1.'.BC'.$2; 	
	}
}

# Assign filehandles
open (FILE1, "$file1") || die "Cant open file: $!";
open (FILE2, "$file2") || die "Cant open file: $!";
open (OUT, ">$output") || die "Cant open file: $!";


# Read in file 1
print "Reading the first file into memory...\n";
my %seen_headers = ();
while (my $line = <FILE1>){
	next if $line !~ m/\>/;
	chomp ($line);
	#my ($header, $count, $barcode) = split (/\#/, $line);
	my ($header, $blah1) =  split (/\//, $line);
	my ($blah2, $count, $barcode) =  split (/\#/, $blah1);
	#my ($header, $blah, $count, $barcode) = split (/[\#\/]/, $line);
	$seen_headers{ $header } = [$barcode, $count];
}

#Read in file 2
print "Iterating through second file...\n";
my $header_found = 'F';
while (my $line = <FILE2>){
	chomp ($line);

	if ($line =~ m/(\>.*)\/(.*)/){
		my ($header, $blah) = ($1, $2);
		if (exists $seen_headers{ $header } ){
			my ($barcode_tmp, $count) = @{$seen_headers {$header}};
			my $seq_obj = Bio::Seq->new(-seq => $barcode_tmp, -alphabet => 'dna');
			my $barcode = ($seq_obj->revcom)->seq;
			print OUT $line, '#', $count, '#', $barcode, "\n";
		}
		else{
			print 'Error: header line from second file not found in first file', "\n";
			print '--- ', $header, "\n";
		}
		$header_found = 'T'
	}
	else{
		if ($header_found eq 'T'){
			my $seq_obj = Bio::Seq->new(-seq => $line, -alphabet => 'dna');
			my $sequence = ($seq_obj->revcom)->seq;
			print OUT $sequence, "\n"
		}
		$header_found = 'F';
	}
}

print 'Done !', "\n";