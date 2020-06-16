#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

sub fix090 {
    my $inTable = getTableName('Frags','090');
    my $outTable = newTable('Frags', '090_new');
    runSQL("INSERT INTO $outTable 
    SELECT FRAGMENTID, FRAGMENTTYPE,
                      CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1, 0 NHITS1,
                      CHROMOSOME2, POSITION2, LENGTH2, STRAND2, DISCREPANCIES2, 0 NHITS2,
                      FRAGMENTSIZE, PAIRID, NFRAGS,
                      EVENTSIZE, STDEVNORMAL, ENDTOLERANCE, 
                      NSETSFRAG, NSETSPAIR 
    FROM $inTable");
}

1;

