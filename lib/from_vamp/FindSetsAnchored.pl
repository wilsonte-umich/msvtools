#!/usr/bin/perl -w
use strict;
use warnings;

#########################################################################
#This contains code for an approach which was never completed...
#########################################################################

use vars(qw(%param %types %fields %refSeqs $uid));

my ($setsTable, $setsFileH);

sub findAnchored{
    my ($sample) = @_;
    
    status("getting statistics...\n");
        my $statsTable = getTableName('Stats', $sample);
        getStatistics($statsTable, \my %stats);
        $param{minFragsSet1} = int($stats{minCoverageFMap} / 2); #/2 to account for mono-allelic events
        $param{minFragsSet1} >= 2 or $param{minFragsSet1} = 2;
         $param{maxFragsSet1} = $stats{maxCoverageFMap};
        ($param{minFragsSet2}, $param{maxFragsSet2}, $param{minFragsSet}, $param{maxFragsSet}) =
            (999999, 0, $param{minFragsSet1}, $param{maxFragsSet1}); 

    status("finding sets...\n");
        $param{fragsTable1} = getTableName('Frags', $sample);
        $setsTable = getTableName('Sets', $sample);
        if(tableExists($setsTable)){
            runSQL("SELECT Max(SETID) FROM $setsTable", \my ($maxSetID));
            fetchRow();
            $param{setID} = $maxSetID + 1;
        } else {
            $setsTable = newTable('Sets', $sample);
            $param{setID} = 1;
        }
        my $setsFile = "$setsTable.csv";        
        open $setsFileH, ">", $setsFile;        
        findSetsAnchored1();
        close $setsFileH;
        
    status("loading sets...\n");        
        loadData($setsFile, $setsTable, ",", "SETID, SETTYPE, CHROMOSOME1, CHROMOSOME2,
                                                SPANSTART, SPANEND, OVERLAPSTART, OVERLAPEND,
                                                STRAND1, STRAND2,
                                                NFRAGSTOTAL, NFRAGSSAMPLE");

    status("estimating fractional overlap...\n");    
        calcRootedFracOverlap($sample);
        
        
        
#find the projections of anchors in sets    

#establish neighbors?


}  

sub findSetsAnchored1{ #used by find, takes only one sample as the comparison input
    $param{setType} = $types{Sets}{Anchored};
    status("  Set type $param{setType}:\n    Chr: ");  
    foreach my $chrom1 (1..nChrom()){    
        $param{chrom1} = $chrom1;
        status("$chrom1 ");
        foreach my $strand1(1,2){        
            $param{strand1} = $strand1;
            my $sql1 = findAnchoredSetsSQL(1, $param{fragsTable1});              
            threadP1(" $sql1 ", \&processAnchoredGroup);
        }
    }
    status("\n");
}

sub findAnchoredSetsSQL{
    my ($sampleN, $fragsTable) = @_;  
    #returns ALL anchored reads of ANY fragment anomaly type (except currently blocking insertion...)
    #NOT limited by NHITS, i.e. does not require that all partner maps be available
    #DOES enforce same strand orientation of anchored reads
    my $anchoredFragsSQL = anchoredFragsSQL($fragsTable);
    return "SELECT FRAGMENTID, PAIRID, FRAGMENTTYPE,
                    CHROMOSOME1, POSITION1, POSITION1 + 10000 POSITION2, FRAGMENTSIZE, 
                    EVENTSIZE, STDEVNORMAL, ENDTOLERANCE, $sampleN AS SAMPLE
            FROM ($anchoredFragsSQL)        
            WHERE ( FRAGMENTTYPE = $types{Sets}{Deletion} OR 
                    FRAGMENTTYPE = $types{Sets}{Inversion} OR 
                    FRAGMENTTYPE = $types{Sets}{Duplication} OR 
                    FRAGMENTTYPE = $types{Sets}{DiffChrom} )
              AND CHROMOSOME1 = $param{chrom1}
              AND STRAND1 = $param{strand1} "; 
#                    FRAGMENTTYPE = $types{Sets}{Insertion} OR 
}

sub anchoredFragsSQL{
    my ($fragsTable) = @_;
    #returns frags for which at least one read had only one mapping (i.e. anchored frags)
    #enforcing orientation of the reads such that read 1 is always an anchored read
    #otherwise equivalent to input Frags table
    return "
        SELECT FRAGMENTID, FRAGMENTTYPE, 
               (CASE NHITS1 WHEN 1 THEN CHROMOSOME1 ELSE (CASE CHROMOSOME2 WHEN 0 THEN CHROMOSOME1 ELSE CHROMOSOME2 END) END) CHROMOSOME1,
               (CASE NHITS1 WHEN 1 THEN POSITION1 ELSE POSITION2 END) POSITION1,
               (CASE NHITS1 WHEN 1 THEN LENGTH1 ELSE LENGTH2 END) LENGTH1,
               (CASE NHITS1 WHEN 1 THEN STRAND1 ELSE STRAND2 END) STRAND1,
               (CASE NHITS1 WHEN 1 THEN DISCREPANCIES1 ELSE DISCREPANCIES2 END) DISCREPANCIES1,
               (CASE NHITS1 WHEN 1 THEN NHITS1 ELSE NHITS2 END) NHITS1,
               (CASE NHITS1 WHEN 1 THEN CHROMOSOME2 ELSE (CASE CHROMOSOME2 WHEN 0 THEN 0 ELSE CHROMOSOME1 END) END) CHROMOSOME2,
               (CASE NHITS1 WHEN 1 THEN POSITION2 ELSE POSITION1 END) POSITION2,
               (CASE NHITS1 WHEN 1 THEN LENGTH2 ELSE LENGTH1 END) LENGTH2,
               (CASE NHITS1 WHEN 1 THEN STRAND2 ELSE STRAND1 END) STRAND2,
               (CASE NHITS1 WHEN 1 THEN DISCREPANCIES2 ELSE DISCREPANCIES1 END) DISCREPANCIES2,
               (CASE NHITS1 WHEN 1 THEN NHITS2 ELSE NHITS1 END) NHITS2,
               FRAGMENTSIZE, PAIRID, NFRAGS, EVENTSIZE, STDEVNORMAL, ENDTOLERANCE, NSETSFRAG, NSETSPAIR
        FROM $fragsTable
        WHERE (NHITS1 = 1 OR NHITS2 = 1)    
    "
}

sub processAnchoredGroup{ 
    my ($fragsRef) = @_;
    $param{setSigs} = {};
    my $setsRef = checkSetPositions($fragsRef, 1);
    foreach my $setRef(@$setsRef){ processAnchoredSet($setRef) }
}
    
    
sub processAnchoredSet{ 
    my ($fragsRef) = @_;
    my $nFragsSet1 = scalar(@$fragsRef);
    ($nFragsSet1 >= $param{minFragsSet1} and $nFragsSet1 <= $param{maxFragsSet1}) or return;
    $param{setID}++;
    printAnchoredSet($fragsRef, $nFragsSet1, 1); 
}

sub printAnchoredSet{ #commit set to files for uploading to db later on
    my ($setRef, $nFragsSet, $sampleN) = @_;
    my @parsedSet = getSet($setRef, $nFragsSet);
    my $setSig = join(":", ($sampleN, @parsedSet[0..4]));  
    unless($param{setSigs}{$setSig}){ #prevent set duplicates      
        print $setsFileH join(",", ($param{setID}, $param{setType}, $param{chrom1}, 0,
                                    @parsedSet[0..3],
                                    $param{strand1}, 0,
                                    $nFragsSet, $nFragsSet))."\n";                                                        
    }
    $param{setSigs}{$setSig} = 1;
}


sub xxxx{
    my $anchoredSets = "SELECT * FROM $setsTable WHERE SETTYPE = $types{Sets}{Anchored}";
    my $anchoredFrags = anchoredFragsSQL($param{fragsTable1});
    my $partners = "SELECT s.SETID, 
                            Decode(f.CHROMOSOME2,0,f.CHROMOSOME1,f.CHROMOSOME2) partnerChrom,
                            Trunc(f.POSITION2 + (f.LENGTH2 / 2)) partnerPos,
                            f.STRAND2 partnerStrand
                     FROM ($anchoredSets) s, ($anchoredFrags) f
                     WHERE f.CHROMOSOME1 = s.CHROMOSOME1
                       AND f.POSITION1 >= s.SPANSTART
                       AND f.POSITION1 <= s.OVERLAPSTART
                       AND f.STRAND1 = s.STRAND1";
    my $rmskTable = "";
    my $rmskSQL = "SELECT p.SETID, r.CLASS_ XXXXETC FROM RMSK
                   FROM ($partners) p, $rmskTable r
                   WHERE p.partnerChrom = r.CHROMOSOME
                     AND p.partnerPos >= r.START_
                     AND p.partnerPos <= r.END_ ";
                     
    #ANY REASON TO do THE SAME for CDS OR NCDS??
    #THOSE ELEMENTS WOULD NOT BE EXPECTED TO BE REPETITIVE

    #query for count of class, or whatever, during use

}

1;

