#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

sub pairSizeHistogram{
    my ($sample) = @_;
    my $pairsTable = getTableName('Pairs', $sample);
    my $histogramTable = getTableName('h', $pairsTable);
    my $outFile = "$histogramTable.csv";
    open my $outFileH, ">", $outFile;
    runSQL("SELECT x bp,
                   sum(decode(series,0,y,0)) convergent,
                   sum(decode(series,1,y,0)) divergent
            FROM (SELECT series, round(x/10)*10 x, sum(y) y
                  FROM $histogramTable
                  GROUP BY series, x)
            WHERE x < 5000 AND series = 0
            GROUP BY x
            ORDER BY x", \my($bp,$conv,$div));
    while(fetchRow()){ print $outFileH "$bp,$conv,$div\n" }
    close $outFileH;
}

sub forcePairStats{
    my ($sample, $minNormal, $modeNormal, $maxNormal, $stDevNormal, $minFragSize, $maxRevNormal) = @_;
    print "$sample, $minNormal, $modeNormal, $maxNormal, $stDevNormal, $minFragSize, $maxRevNormal\n";
    my $statsTable = getTableName('Stats', $sample);
    updateStat($statsTable, 'meanNormal', $modeNormal);
    updateStat($statsTable, 'stDevNormal', $stDevNormal);
    updateStat($statsTable, 'minNormal', $minNormal);
    updateStat($statsTable, 'maxNormal', $maxNormal);
    updateStat($statsTable, 'modeNormal', $modeNormal);
    updateStat($statsTable, 'minFragSize', $minFragSize);
    updateStat($statsTable, 'maxReverseNormal', $maxRevNormal);
}



1;

