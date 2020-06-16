#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));

my (%inReadFiles, %inMapFiles);

sub runL1{
    createDatabase('chr1L1');
    my %readLengths = ('090a'=>36, '090b'=>35, A1A1a=>36, A1A1b=>35, A3A2ab=>35, A3A2c=>36, D2a=>39, D2b=>39);
    foreach my $sample(qw(090a 090b A1A1a A1A1b A3A2ab A3A2c D2a D2b)){
        getReadFiles($sample, 'purgedFasta', \%inReadFiles);
        getMapFiles($sample, \%inMapFiles);   
        $inMapFiles{read1} .= ".chr1L1";
        $inMapFiles{read2} .= ".chr1L1";
        $param{readLength} = $readLengths{$sample};
        mapReads_PASS_Rooted($sample, 'read1', \%inReadFiles, \%inMapFiles, 'chr1L1', 1);
        mapReads_PASS_Rooted($sample, 'read2', \%inReadFiles, \%inMapFiles, 'chr1L1', 1);
        
    }

}

1;


