#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command));

#scaleSample.pl provides wholly generic internal functions
#for placing two hit count samples onto the same numeric scale
#the sample with the lower hit count is scaled up to match
#the sample with the higher hit count

sub scaleSamples {
    my ($sample1, $sample2, $statsType) = @_;
    if($statsType){
        $sample1 = "$statsType\_$sample1";
        $sample2 = "$statsType\_$sample2"
    }
    my $statsTable1 = getTableName('Stats', $sample1);
    my $statsTable2 = getTableName('Stats', $sample2);
    getStatistics($statsTable1,\my%stats1);
    getStatistics($statsTable2,\my%stats2);
    my $count = "posCount";
    $param{keepDups} and $count = "hitCount";
    my $sample1HitCount = $stats1{$count};
    my $sample2HitCount = $stats2{$count};
    my $scalarNum = $sample1HitCount;
    $scalarNum < $sample2HitCount and $scalarNum = $sample2HitCount;
    my $sample1Scalar = ($scalarNum/$sample1HitCount);
    my $sample2Scalar = ($scalarNum/$sample2HitCount);
    return ($sample1Scalar, $sample2Scalar, $sample1HitCount, $sample2HitCount);
}

1;

