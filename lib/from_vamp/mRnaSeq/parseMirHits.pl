#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{parseMirHits} or $command{parseMirHits} = ['multiThread', '48:00:00', 2000, 0];

sub parseMirHits {
    my ($sample) = @_;
    my $featureName = "miR";
    my $mapTable = newTable('EMap', "$featureName\_$sample");
    calculateFeatureDensity($sample, $mapTable, $featureName, \&getInMirHitsSql);
    calculateNormalizedFeatureDensity($sample, $mapTable, $featureName, \&getGenomeMirSize);
    my $histogramTable = newTable('h', $mapTable);
    createHitMapHistogram($histogramTable, $mapTable);
}
sub getInMirHitsSql {
    my ($sample, $chrom) = @_;
    my $mirTable = "miR_$param{refSeq}";
    return getHitJoinSQL($sample, $chrom, $mirTable, "MIRID", "chromosome", "start_ - $param{padding}", "end_ + $param{padding}", "strand");
}
sub getGenomeMirSize {
    my $mirTable = "miR_$param{refSeq}";
    runSQL("SELECT sum((end_ + $param{padding}) - (start_ - $param{padding}) + 1) N from $mirTable", \my$size);
    fetchRow();
    status("    $param{refSeq} has $size bp in miRs\n");
    return $size; 
}

1;

