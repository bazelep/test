#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $file = $ARGV[0];


my $barcode = (split(/\./, basename($file)))[1];
print "$file\t$barcode\n";



open FILE, $file or die $!;
open OUT, ('>' . $file . ".updated") or die $!;

while (my $line = <FILE>) {
    chomp $line;

    if ($line =~ /^#CHROM/) {
	    my @cols = split(/\t/, $line);
	
        ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  K0
	    $cols[9] = $barcode;

	    print OUT (join("\t", @cols) . "\n");
    }
    else {
	    print OUT "$line\n";
    }
}
close FILE;

exit;
