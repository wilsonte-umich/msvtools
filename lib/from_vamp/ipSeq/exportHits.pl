#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{exportHits} or $command{exportHits} = ['multiThread', '4:00:00', 2000, 0];

sub exportHits { 
    my ($sample) = @_;
    my $hitsTable = getTableName('Hits', $sample);
    my $hitsFile = "$hitsTable.csv";
    open my $hitsFileH, ">", $hitsFile or die "could not open $hitsFile\n";
    print $hitsFileH "CHROMOSOME,POSITION,STRAND,COUNT\n";
    runSQL("SELECT CHROMOSOME, POSITION, STRAND, COUNT_
            FROM $hitsTable
            ORDER BY CHROMOSOME, POSITION, STRAND",
            \my($chrom,$pos,$strand,$count));
    while(fetchRow()){ print $hitsFileH "$chrom,$pos,$strand,$count\n" }
    close $hitsFileH;
}

1;
