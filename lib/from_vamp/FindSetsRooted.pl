#!/usr/bin/perl -w
use strict;
use warnings;

#########################################################################
#FindSetsRooted.pl finds one-ended Sets within the partner reads corresponding
#to the reads rooted in a repeat or other element when parameter rooted is used.
#One-ended Sets are then subjected to neighbor analysis to find forward-reverse
#combination Sets within a specified maxNeighborDistance.
#########################################################################

use vars(qw(%param %types %fields %refSeqs $uid));

sub findSetsRooted1{ #used by find, takes only one sample as the comparison input
    $param{setType} = $types{Sets}{Rooted};
    my @strands = initializeFindSets();    
    status("  Set type $param{setType}:\n    Chr: ");    
    foreach my $chrom1 (1..$refSeqs{$param{refSeqBase}}{nChrom}){    
        $param{chrom1} = $chrom1;
        status("$chrom1 ");
        foreach my $strandsRef(@strands){
            $param{strand1} = $$strandsRef[0];
            $param{strand2} = $$strandsRef[1];
            my $sql1 = findSetsSQL(1, $param{fragsTable1});      
            threadP1(" $sql1 ", \&processRootedGroup);
        }
    }
    status("\n");
}

sub findSetsRooted2{ #used by compare, takes two samples as the comparison input
    $param{setType} = $types{Sets}{Rooted};
    my @strands = initializeFindSets();
    status("  Set type $param{setType}:\n    Chr: ");
    foreach my $chrom1 (1..$refSeqs{$param{refSeqBase}}{nChrom}){ 
        $param{chrom1} = $chrom1;
        status("$chrom1 ");
        foreach my $strandsRef(@strands){
            $param{strand1} = $$strandsRef[0];
            $param{strand2} = $$strandsRef[1];
            #temporarily merge sample1 and sample2 into a single search table for set finding
            my $sql1 = findSetsSQL(1, $param{fragsTable1});
            my $sql2 = findSetsSQL(2, $param{fragsTable2});
            threadP1(" ($sql1 UNION ALL $sql2) ", \&processRootedGroup);
        }
    }
    status("\n");
}

sub processRootedGroup{ 
    my ($fragsRef) = @_;
    $param{setSigs} = {};
    my $setsRef = checkSetPositions($fragsRef, 1);
    foreach my $setRef(@$setsRef){ processSet($setRef) }
}

sub calcRootedFracOverlap{
    my ($sample) = @_;
    my $setsTable = getTableName('Sets', $sample); 
    my $statsTable = getTableName('Stats', $sample);    
    getStatistics($statsTable, \my %stats);
    my $midNormal = (($stats{maxNormal} - $stats{minNormal}) / 2) + $stats{minNormal};
    runSQL("UPDATE $setsTable 
            SET FRACTIONOVERLAPLEFT = ($midNormal - (OVERLAPSTART - SPANSTART)) / ($midNormal + (OVERLAPSTART - SPANSTART))");         
}

sub findRootedNeighbors1{ #used by find, takes only one sample as the comparison input
    my ($sample) = @_;
    $param{setType} = $types{Sets}{Rooted};
    my ($setsTable) = getTableName("Sets", $sample);    
    my $neighborsTable = newTable("Sets", "$sample\_N"); 
    my $neighborsFile = "$neighborsTable.csv";  
    open my $neighborsFileH, ">", $neighborsFile;    
    $param{setID} = 0;
    foreach my $chrom1 (1..$refSeqs{$param{refSeqBase}}{nChrom}){
        status("$chrom1 "); 
        my $sql = "SELECT *
                FROM $setsTable
                WHERE CHROMOSOME1 = $chrom1
                  AND FRACBADFRAGS <= 0.4
                  AND FRACTIONOVERLAPLEFT <= 0.9
                ORDER BY SPANSTART";
        runSQL($sql); #many chromosomes won't return any sets; check whether sql has any return
        my @test = fetchAllHashRef();
        scalar(@test) or next;   
        runSQL($sql); 
        my @fRefs; #get the first forward Set end
        while(!(scalar @fRefs) and my $__ = fetchRowHashRef()){ $$__{STRAND1} == 1 and push @fRefs, $__ }
        while (my $__ = fetchRowHashRef()){ 
            if ($$__{STRAND1} == 1){ #hold forward Set ends for checking later
                push @fRefs, $__;
            } else { #check reverse Set ends against all persistent forward Set ends
                my @fRefsHold; #to hold the forward set ends still withing the distance tolerance
                foreach my $fRef(@fRefs){
                    if($$__{OVERLAPSTART} - $$fRef{SPANSTART} < $param{maxNeighborDistance}){ #overlapStart - spanStart is the span of the candidate combined Set
                        createFullRootedSet($neighborsFileH, $fRef, $__);
                        push @fRefsHold, $fRef;  #reference forward Set end has still passed the distance tolerance, keep it
                    }
                }
                @fRefs = @fRefsHold;
            }
        }        
    }
    close $neighborsFileH; 
    loadData($neighborsFile, $neighborsTable, ",", "SETID, SETTYPE, CHROMOSOME1, CHROMOSOME2,
                                                    SPANSTART, SPANEND, OVERLAPSTART, OVERLAPEND,
                                                    FRAGMENTMEAN, EVENTMEAN, STDEV, 
                                                    MEANDELTA, FRACTIONDELTA2, FRACTIONOVERLAPLEFT, FRACTIONOVERLAPRIGHT,
                                                    STRAND1, STRAND2,
                                                    NFRAGSTOTAL, NFRAGSSAMPLE, FRACBADFRAGS,
                                                    MARK");  
    status("\n");
}

sub createFullRootedSet{
    my ($neighborsFileH, $fRef, $rRef) = @_;
    $param{setID}++;
    my $fracPromiscuous = getCombinedFraction($fRef, $rRef, 'FRACBADFRAGS');
    my $fracWithNormal = getCombinedFraction($fRef, $rRef, 'FRACTIONDELTA2'); #FRACTIONDELTA2 holds FRACTIONWITHNORMAL for Rooted Sets
    print $neighborsFileH join(",", $param{setID}, $$fRef{SETTYPE}, $$fRef{CHROMOSOME1}, 0, 
                                    $$fRef{SPANSTART}, $$rRef{OVERLAPSTART}, $$fRef{OVERLAPSTART}, $$rRef{SPANSTART},
                                    0, $$rRef{SPANSTART} - $$fRef{OVERLAPSTART}, 0, 
                                    0, $fracWithNormal, $$fRef{FRACTIONOVERLAPLEFT}, $$rRef{FRACTIONOVERLAPLEFT},
                                    $$fRef{STRAND1}, $$rRef{STRAND1},
                                    $$fRef{NFRAGSTOTAL} + $$rRef{NFRAGSTOTAL}, $$fRef{NFRAGSSAMPLE} + $$rRef{NFRAGSSAMPLE}, $fracPromiscuous,
                                    0)."\n";
}

sub getCombinedFraction{
    my ($fRef, $rRef, $fractionField) = @_;
    my $fN = $$fRef{$fractionField} * $$fRef{NFRAGSSAMPLE};
    my $rN = $$rRef{$fractionField} * $$rRef{NFRAGSSAMPLE};
    return ($fN + $rN)/($$fRef{NFRAGSSAMPLE} + $$rRef{NFRAGSSAMPLE});
}

##obsolete now that mapRooted tries to map every root read into target genome region
# sub markKnownRoots{
    # my ($sample) = @_;
    # my $rootTable = getTableName('ROOT', "$param{refSeq}\_$param{rooted}");
    # tableExists($rootTable) or (print "could not find table $rootTable" and return);     
    # my $neighborsTable = getTableName("Sets", "$sample\_N");  
    # runSQL("UPDATE $neighborsTable SET Mark = 0");     
    # my $setsCrossingRootSQL = getSetsCrossingRootSQL($neighborsTable, $rootTable, 'SPANSTART', 'SPANEND');     
    # runSQL("UPDATE $neighborsTable SET Mark = 1 WHERE SETID IN ($setsCrossingRootSQL)");
    # $setsCrossingRootSQL = getSetsCrossingRootSQL($neighborsTable, $rootTable, 'OVERLAPSTART', 'OVERLAPEND');           
    # runSQL("UPDATE $neighborsTable SET Mark = 2 WHERE SETID IN ($setsCrossingRootSQL)");
# }
# 
# sub getSetsCrossingRootSQL{
    # my ($neighborsTable, $rootTable, $start, $end) = @_;
    # return "SELECT n.SETID
            # FROM $neighborsTable n,
                 # $rootTable r
            # WHERE n.CHROMOSOME1 = r.CHROMOSOME
              # AND n.$start < r.END_
              # AND n.$end > r.START_
            # GROUP BY SETID";
# }

1;

