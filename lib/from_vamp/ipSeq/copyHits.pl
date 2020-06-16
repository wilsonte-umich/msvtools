#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{copyHits} or $command{copyHits} = ['singleThread', '1:00:00', 1000, 0];

sub copyHits {
    my ($sampleOut, $sampleIn) = @_;
    copyTable('Hits', $sampleOut, $sampleIn);
    copyTable('HMap', $sampleOut, $sampleIn);
}

1;

