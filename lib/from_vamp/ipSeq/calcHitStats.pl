#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{calcHitStats} or $command{calcHitStats} = ['multiThread', '4:00:00', 1000, 0];

sub calcHitStats { 
    my ($sample) = @_;
    my $hitsTable = getTableName('Hits', $sample);
    status("counting unique and duplicate hit counts for $hitsTable...\n");
    my $statsTable = newTable('Stats', $hitsTable);
    runSQL("SELECT count(*) posCount, sum(count_) hitCount
            FROM $hitsTable", \my($posCount, $hitCount));
    fetchRow();
    $hitCount or die "retrieved no hits from $hitsTable\n";
    updateStat($statsTable, 'posCount', $posCount);
    updateStat($statsTable, 'hitCount', $hitCount);
    status("  unique hit count = $posCount\n");
    status("  duplicate hit count = $hitCount\n");
}

1;
