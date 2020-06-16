#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

sub btest{
    foreach my $read (1..2){
        my $inFile = "/home/wilsonte/vamp/data/D2a/$read/D2a_$read.fa.p";
        my $outFile = "/home/wilsonte/vamp/data/btest/$read/btest_$read.fa.p";
        open my $inFileH, "<", $inFile or die "could not open $inFile";
        open my $outFileH, ">", $outFile or die "could not open $outFile";
        my $i = 0;
        while ($i < 2E6 and my $line = <$inFileH>){
            print $outFileH $line;
            $i++;
        }
        close $inFileH;
        close $outFileH;
    }
}

1;

