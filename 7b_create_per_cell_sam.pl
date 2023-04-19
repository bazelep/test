#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;


my $sample = $ARGV[0];
my $barcode_file = $ARGV[1];
my $sam_file = $ARGV[2];
my $sam_header_file = $ARGV[3];
my $BC_DIR = $ARGV[4];
my %sam;
my %bar;

die "\n\n\t:Usage: perl $0 <sample> <barcode file> <sam file> <sam header file> <barcode dir>\n\n" 
unless (defined($sample));


open BAR, $barcode_file or die $!;
while (my $line = <BAR>) {
    chomp $line;
    $bar{$line} = 1;
}
close BAR;
#print Dumper(\%bar);
#exit;

open SAM, $sam_file or die $!;
while (my $line = <SAM>) {
    chomp $line;

    my @cols = split(/\t/, $line);

    foreach my $col (@cols) {
	if ($col =~ /CB:Z:/) {
	    
	    $col =~ s/^CB:Z://;
	    
	    last unless exists($bar{$col});
	    
	    push @{$sam{$col}}, $line;
	}
    }
}
close SAM;

open HEADER, $sam_header_file or die $!;
my @header = <HEADER>;
close HEADER;


#print Dumper(\%sam);

print "there are " . (scalar keys %sam) . " barcodes\n";
#exit;

foreach my $barcode (keys %sam) {

    my $outfile = "$BC_DIR/$sample.$barcode.sam";
    open OUT, ('>' . $outfile) or die $!;

    foreach my $header_line (@header) {
	    print OUT $header_line;
    }
    
    foreach my $line (@{$sam{$barcode}}) {
	    print OUT "$line\n";
    }
    close OUT;
}


exit;
