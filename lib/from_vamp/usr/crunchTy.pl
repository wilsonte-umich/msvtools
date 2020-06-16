#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));

sub crunchTy{
    my $sample = "vac22wt";
    open my $outputH, ">", "$sample"."TyPairs.xls";
    foreach my $chrom(1..nChrom()){
        runSQL("select f.spanstart start_, r.overlapstart end_
                from
                (select * from sets_$sample"."_ty where strand1 = 1) f,
                (select * from sets_$sample"."_ty where strand1 = 2) r
                where f.chromosome1 = $chrom
                 and f.chromosome1 = r.chromosome1
                 and r.overlapstart - f.spanstart > 0
                 and r.overlapstart - f.spanstart < 20000
                 and r.spanstart - f.overlapstart > -500
                 and r.spanstart - f.overlapstart < 2000
                order by f.chromosome1, f.spanstart, r.overlapstart", \my($start, $end));
        my ($minStart, $maxHigh) = (0, 0);
        while(fetchRow()){
            if($minStart and $start > $maxHigh){
                print $outputH "$chrom\t$minStart\t$maxHigh\n";
                ($minStart, $maxHigh) = (0, 0);
            }
            $minStart or $minStart = $start;
            $maxHigh >= $end or $maxHigh = $end;
        }
        print $outputH "$chrom\t$minStart\t$maxHigh\n";      
    }      
    close  $outputH;
}

1;



