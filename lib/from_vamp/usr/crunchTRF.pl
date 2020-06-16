#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs));

my($vampID, $chrom_, $start, $end, $spanStart, $spanEnd, $outFileH);

sub crunchTRF{
    my $outTable = "TRF_HG18";
    my $outFile = "$outTable.csv";
    foreach my $chrom(1..nChrom()){    
        open $outFileH, ">", $outFile;
        $chrom_ = $chrom;
        runSQL("SELECT START_, END_ FROM SIMPLEREPEAT_HG18 WHERE CHROMOSOME = $chrom ORDER BY START_", \($start, $end));
        ($spanStart, $spanEnd) = (undef, 0);
        while(fetchRow()){
            if(defined $spanStart and $start > $spanEnd){
                commitTRFSpan();
                ($spanStart, $spanEnd) = (undef, 0);
            }
            defined $spanStart or $spanStart = $start;
            $spanEnd >= $end or $spanEnd = $end;
        }   
        commitTRFSpan();
        close $outFileH;
        loadData($outFile, $outTable, ",", "VAMPID, CHROMOSOME, START_, END_");
    }
}
sub commitTRFSpan{
    $vampID++;
    print $outFileH "$vampID,$chrom_,$spanStart,$spanEnd\n";
}

1;












