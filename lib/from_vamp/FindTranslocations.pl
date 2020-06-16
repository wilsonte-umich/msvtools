#!/usr/bin/perl -w
use strict;
use warnings;

#########################################################################
#FindTranslocations.pl is very similar to FindSets.pl, except modified
#to find sets among pairs mapped to different chromosomes.
#Differences include that Fragment overlap has no meaning for DiffChrom Sets.
#########################################################################

use vars(qw(%param %types %fields %refSeqs $uid));

sub findTranslocations1{  #used by find, takes only one sample as the comparison input
    initializeFindTranslocations();
    status("  Set type $param{setType}:\n    Chr: ");
    #pair each chromosome with all other non-identical chromosomes, in turn
    #considering all possible strand orientations independently
    foreach my $chrom1 (1..$refSeqs{$param{refSeqBase}}{nChrom}){ 
        $param{chrom1} = $chrom1;
        status("$chrom1 ");
        foreach my $strand1(1..2){
            $param{strand1} = $strand1;
            foreach my $chrom2 (1..$refSeqs{$param{refSeqBase}}{nChrom}){  
                $param{chrom2} = $chrom2;
                foreach my $strand2(1..2){ 
                    $param{strand2} = $strand2;
                    my $sql1 = findTranslocationsSQL(1, $param{fragsTable1});
                    threadP1(" $sql1 ", \&threadP2);           
                }
            }        
        }
    }
    status("\n");
}

sub findTranslocations2{ #used by compare, takes two samples as the comparison input
    initializeFindTranslocations();
    status("  Set type $param{setType}:\n    Chr: ");
    #temporarily merge sample1 and sample2 into a single search table for set finding AND
    #pair each chromosome with all other non-identical chromosomes, in turn   
    #considering all possible strand orientations independently  
    foreach my $chrom1 (1..$refSeqs{$param{refSeqBase}}{nChrom}){
        $param{chrom1} = $chrom1;
        status("$chrom1 ");
        foreach my $strand1(1..2){
            $param{strand1} = $strand1;
            foreach my $chrom2 (1..$refSeqs{$param{refSeqBase}}{nChrom}){
                $param{chrom2} = $chrom2;
                foreach my $strand2(1..2){ 
                    $param{strand2} = $strand2;
                    my $sql1 = findTranslocationsSQL(1, $param{fragsTable1});
                    my $sql2 = findTranslocationsSQL(2, $param{fragsTable2});
                    threadP1(" ($sql1 UNION ALL $sql2) ", \&threadP2);
                }
            }         
        }
    }
    status("\n");
}

sub initializeFindTranslocations{
    $param{setType} = $types{Sets}{DiffChrom};
    $param{enforceOverlap} = 0; #enforceOverlap set to false for DiffChrom 
}

sub findTranslocationsSQL{
    my ($sampleN, $fragsTable) = @_; 
    return "SELECT FRAGMENTID, PAIRID, FRAGMENTTYPE,
                    CHROMOSOME1, POSITION1, POSITION2, FRAGMENTSIZE, 
                    EVENTSIZE, STDEVNORMAL, ENDTOLERANCE, $sampleN AS SAMPLE
            FROM $fragsTable
            WHERE FRAGMENTTYPE=$param{setType}
                AND FRAGMENTID > 0
                AND CHROMOSOME1=$param{chrom1}
                AND CHROMOSOME2=$param{chrom2}
                AND STRAND1=$param{strand1}
                AND STRAND2=$param{strand2} 
                AND NHITS1 <= $param{maxHits} AND NHITS2 <= $param{maxHits} 
                $param{strataFilter}";
                #FRAGMENTID > 0 accounts for Fragment masking by reduceNoise
}

1;
