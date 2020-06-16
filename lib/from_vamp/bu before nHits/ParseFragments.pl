#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#ParseFragments.pl parses Pairs into Fragments of defined types.
#Output = Fragments tables.
#########################################################################

use vars(qw(%param %types %fields %fieldNames %refSeqs));

my $pairsTable;
my $fragsTable;
my $statsTable;
my %stats;

sub parseFragments{
    my ($sample) = @_;

    $pairsTable = getTableName('Pairs', $sample);
     
    status("getting statistics...\n");
        $statsTable = getTableName('Stats', $sample);
        getStatistics($statsTable, \%stats);

    status("parsing fragments...\n");
        $fragsTable = newTable('Frags', $sample);
        executeParse();

    status("cleaning up fragments...\n");
        cleanUpFragments($sample);          

    status("counting normal fragments...\n");
        countNormals();
        
    status("creating fragment size histogram...\n");
        my $histogramTable = newTable('h', $fragsTable);
        createFragSizeHistogram($histogramTable, $fragsTable);

    status("creating coverage maps...\n");
        my ($fragMapTable, $readMapTable) = createMaps($sample, $fragsTable, $statsTable);

    status("calculating coverages...\n");
        my ($fragCoverage, $readCoverage) = calculateCoverages($fragsTable, $statsTable);
        $fragMapTable and $fragCoverage and calculateNormalizedCoverage($fragMapTable, 'FMap', $fragCoverage);
        $readMapTable and $readCoverage and calculateNormalizedCoverage($readMapTable, 'RMap', $readCoverage);
        
    if ($param{dropInput}){dropTable($pairsTable)} 
}

sub callExecuteParse{ #for remote calls
    my $statsRef;
    ($pairsTable, $fragsTable, $statsRef) = @_;
    %stats = %{$statsRef};
    executeParse();
}

sub executeParse{
    #copies PAIRS into FRAGS, creating new FRAGMENTTYPE based on PAIRTYPE and FRAGMENTSIZE
    my $timeStamp = getTimeStamp();
    my $diffChromCase = " WHEN $types{Pairs}{DiffChrom} THEN $types{Frags}{DiffChrom} ";
    $param{noDiffChrom} and $diffChromCase = ' '; 
    my $reverseNormalCase = ' ';
    $param{isCircles} and $reverseNormalCase = " WHEN FRAGMENTSIZE < $stats{maxReverseNormal} THEN $types{Frags}{ReverseNormal} ";
    runSQL("INSERT INTO $fragsTable
                SELECT (CHROMOSOME1||POSITION1||LENGTH1||STRAND1||CHROMOSOME2||POSITION2||LENGTH2||STRAND2||$timeStamp) FRAGMENTID,
                        (CASE PAIRTYPE
                            WHEN $types{Pairs}{Convergent} THEN
                                (CASE
                                    WHEN FRAGMENTSIZE < $stats{minFragSize} THEN $types{Frags}{TooSmall}
                                    WHEN FRAGMENTSIZE < $stats{minNormal} THEN $types{Frags}{Insertion}
                                    WHEN FRAGMENTSIZE < $stats{maxNormal} THEN $types{Frags}{Normal}
                                    ELSE $types{Frags}{Deletion}
                                END)
                            WHEN $types{Pairs}{Colinear} THEN
                                (CASE
                                    WHEN FRAGMENTSIZE < $stats{minFragSize} THEN $types{Frags}{TooSmall}
                                    ELSE $types{Frags}{Inversion}
                                END)
                            WHEN $types{Pairs}{Divergent} THEN
                                (CASE
                                    WHEN FRAGMENTSIZE < $stats{minFragSize} THEN $types{Frags}{TooSmall}
                                    $reverseNormalCase
                                    ELSE $types{Frags}{Duplication}
                                END)
                            WHEN $types{Pairs}{SingleRead} THEN $types{Frags}{SingleRead}
                            $diffChromCase
                            ELSE -1
                        END) FRAGMENTTYPE,
                        CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1,
                        CHROMOSOME2, POSITION2, LENGTH2, STRAND2, DISCREPANCIES2,
                        FRAGMENTSIZE, PAIRID, NFRAGS,
                        ($stats{modeNormal} - FRAGMENTSIZE) EVENTSIZE,
                        $stats{stDevNormal} STDEVNORMAL,
                        (CASE PAIRTYPE
                            WHEN $types{Pairs}{Convergent} THEN
                                (CASE
                                    WHEN FRAGMENTSIZE < $stats{minNormal} THEN FRAGMENTSIZE
                                    ELSE $stats{maxNormal}
                                END)
                            ELSE $stats{maxNormal}
                        END) ENDTOLERANCE,
                        1 NSETSFRAG, 0 NSETSPAIR
                FROM $pairsTable  ");               
}

sub cleanUpFragments{
    my ($sample) = @_;
    #Logical flow for purging of false and duplicate mappings outlined below
    
    #STEP 1a: Purge Fragments with Discrepancies (i.e. imperfect) derived from Pairs 
    #that also gave rise to Fragments without Discrepancies (i.e. perfect)
    #Parameter preferPerfect determines whether this executes before or after
    #purge to Normal, or not at all.
    purgeImperfect(1);
    
    #STEP 2: Purge false anomalies, i.e. anomalies which also have
    #        a corresponding Normal, ReverseNormal or TooSmall
    #        Priority order for keeping = Normal, then ReverseNormal, then TooSmall, then anomaly
    #1a: purge anomalies, TooSmall and ReverseNormal when Normal is also present (Normal persists)
    purgeFalseWithNormal();
    #1b: purge anomalies and TooSmall when ReverseNormal or TooSmall is also present (Normal and ReverseNormal persist)
    purgeFalseWithReverseNormal();
    
    #STEP 1b: as above
    purgeImperfect(2);

    #At this point, any pair will only have Normal or ReverseNormal or anomaly fragments.
    #However, there may be more than one Fragment of a type per Pair, and different Pairs may be duplicates

    #STEP 3: Find the best Fragment when a Pair gave rise to more than one Normal or ReverseNormal
    #        Priority for keeping = least number of discrepancies, closest to mean Fragment size (Normal only)
    #        If a best Fragment cannot be found, discard the entire Pair as non-informative
    #        All anomalies not purged above persist prior to Set finding
    pickBestPair($sample, 'Normal');
    pickBestPair($sample, 'ReverseNormal');
    pickBestSingleRead($sample);

    #At this point, a Pair will have either:
    #        1 (and only 1) Normal, or
    #        1 (and only 1) ReverseNormal, or
    #        potentially multiple anomalies, or
    #        no remaining fragments.
    #However, there may still be duplicates, i.e. the different pairs mapping to the same coordinates
    #(such duplicate must have different discrepancies - identical read pairs are purged during prepareReads)

    #STEP 4: Find the best Fragment when more than one Pair mapped to the same coordinates
    #        Priority for keeping = least number of discrepancies
    purgeDuplicateFragments();
}

sub purgeImperfect{
    my ($preferPerfect) = @_;
    ($param{preferPerfect} and $param{preferPerfect} == $preferPerfect) or return;
    status("  purging imperfect Fragments with perfect Fragments from same Pair\n");    
    runSQL("DELETE FROM $fragsTable    
            WHERE (DISCREPANCIES1 > 0 OR DISCREPANCIES2 > 0)
            AND PAIRID IN
                (SELECT PAIRID FROM $fragsTable
                    WHERE (DISCREPANCIES1 = 0 AND DISCREPANCIES2 = 0)
                    GROUP BY PAIRID)");
}

sub purgeFalseWithNormal{
    status("  purging Fragments with Normal from same Pair\n");
    runSQL("DELETE FROM $fragsTable
            WHERE FRAGMENTTYPE != $types{Frags}{Normal}
            AND PAIRID IN
                (SELECT PAIRID FROM $fragsTable
                    WHERE FRAGMENTTYPE = $types{Frags}{Normal}
                    GROUP BY PAIRID)");
}

sub purgeFalseWithReverseNormal{
    status("  purging Fragments with ReverseNormal from same Pair\n");
    #TooSmall are deliberately deleted along with all anomalies that also have a TooSmall
    #TooSmall will never be used again, were only progagated to allow this anomaly purging
    runSQL("DELETE FROM $fragsTable
            WHERE (FRAGMENTTYPE != $types{Frags}{Normal}
                AND FRAGMENTTYPE != $types{Frags}{ReverseNormal})
            AND PAIRID IN
                (SELECT PAIRID FROM $fragsTable
                    WHERE FRAGMENTTYPE = $types{Frags}{ReverseNormal}
                       OR FRAGMENTTYPE = $types{Frags}{TooSmall}
                    GROUP BY PAIRID)");
}

sub pickBestPair{
    my ($sample, $frType) = @_;
    status("  picking best Fragment from Pairs with multiple $frType\n");
    my $fragmentType = $types{Frags}{$frType};
    my $fragsTableTMP = newTable('Frags', "$sample\_TMP");
    my $fragsFileTMP = "$fragsTableTMP.csv";
    open my $fragsFileTMPH, ">", $fragsFileTMP;
    #DBI note: pickBestPair uses a hash reference to track Frags in Pair group
    #therefore, SQL column binding to scalars is not possible, i.e must use fetchRowHashRef() 
    runSQL("SELECT *
            FROM $fragsTable
            WHERE FRAGMENTTYPE=$fragmentType
            AND PAIRID IN
                (SELECT PAIRID
                    FROM (SELECT PAIRID, Count(PAIRID) as N
                            FROM $fragsTable
                            WHERE FRAGMENTTYPE=$fragmentType
                            GROUP BY PAIRID)
                    WHERE N>1)
            ORDER BY PAIRID");
    my @currentRowRefs;
    my $prevPairID = 0;
    while (my $__ = fetchRowHashRef()){
        if ($$__{PAIRID} > $prevPairID and $prevPairID){
            checkCurrentRowRefs($fragsFileTMPH, $fragmentType, @currentRowRefs);
            @currentRowRefs = ();
        }
        push @currentRowRefs, $__;
        $prevPairID = $$__{PAIRID};
    }
    checkCurrentRowRefs($fragsFileTMPH, $fragmentType, @currentRowRefs);
    close $fragsFileTMPH;
    loadData($fragsFileTMP, $fragsTableTMP, ",", "FRAGMENTID, FRAGMENTTYPE,
            CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1,
            CHROMOSOME2, POSITION2, LENGTH2, STRAND2, DISCREPANCIES2,
            FRAGMENTSIZE, PAIRID, NFRAGS, 
            EVENTSIZE, STDEVNORMAL, ENDTOLERANCE, 
            NSETSFRAG, NSETSPAIR");
    runSQL("DELETE 
            FROM $fragsTable
            WHERE FRAGMENTTYPE=$fragmentType
            AND PAIRID IN
                (SELECT PAIRID
                    FROM (SELECT PAIRID, Count(PAIRID) as N
                            FROM $fragsTable
                            WHERE FRAGMENTTYPE=$fragmentType
                            GROUP BY PAIRID)
                    WHERE N>1) ");      
    runSQL("INSERT INTO $fragsTable SELECT * FROM $fragsTableTMP");
    dropTable($fragsTableTMP);
}

sub checkCurrentRowRefs{
    my ($fragsFileH, $fragmentType, @rowRefs) = @_;  #filter by sequence matches first, then fragment size
    if (scalar @rowRefs > 1){@rowRefs = checkByNDiscs(@rowRefs)}
    if ($fragmentType == $types{Frags}{Normal} and scalar @rowRefs > 1){@rowRefs = checkBySizeDelta(@rowRefs)}
    if (scalar @rowRefs == 1){commitRow($fragsFileH, @rowRefs)} #if still can't find best Fragment, drop the Pair
}

sub checkByNDiscs{
    my (@rowRefs) = @_;
    my @currentRowRefs;    
    my $minNDisc = 999999;
    foreach my $rowRef (@rowRefs){
        my $nDisc = getNDisc($$rowRef{DISCREPANCIES1}) + getNDisc($$rowRef{DISCREPANCIES2});
        if ($nDisc < $minNDisc){
            $minNDisc = $nDisc;
            @currentRowRefs = ($rowRef);
        } elsif ($nDisc = $minNDisc){
            push @currentRowRefs, $rowRef; 
        }
    }
    return @currentRowRefs;
}

sub getNDisc{
    my ($discs) = @_;
    if ($discs){
        my $l = length($discs);
        my $n = int($l / 4);
        if ($l % 4){$n++}
        return $n;
    } else {
        return $discs 
    }   
}

sub checkBySizeDelta{
    my (@rowRefs) = @_;
    my @currentRowRefs;    
    my $minSizeDelta = 999999;
    foreach my $rowRef (@rowRefs){
        my $sizeDelta = abs($$rowRef{EVENTSIZE});
        if ($sizeDelta < $minSizeDelta){
            $minSizeDelta = $sizeDelta;
            @currentRowRefs = ($rowRef);
        } elsif ($sizeDelta = $minSizeDelta){
            push @currentRowRefs, $rowRef; 
        }
    }
    return @currentRowRefs;
}

sub commitRow{
    my ($fragsFileH, $rowRef) = @_;
    print $fragsFileH join(",",
        ($$rowRef{FRAGMENTID}, $$rowRef{FRAGMENTTYPE},
        $$rowRef{CHROMOSOME1}, $$rowRef{POSITION1}, $$rowRef{LENGTH1}, $$rowRef{STRAND1}, $$rowRef{DISCREPANCIES1},
        $$rowRef{CHROMOSOME2}, $$rowRef{POSITION2}, $$rowRef{LENGTH2}, $$rowRef{STRAND2}, $$rowRef{DISCREPANCIES2},
        $$rowRef{FRAGMENTSIZE}, $$rowRef{PAIRID}, $$rowRef{NFRAGS},
        $$rowRef{EVENTSIZE}, $$rowRef{STDEVNORMAL}, $$rowRef{ENDTOLERANCE},
        $$rowRef{NSETSFRAG}, $$rowRef{NSETSPAIR}
        ))."\n" ;  
}

sub pickBestSingleRead{
    my ($sample) = @_;
    status("  picking best SingleReads\n");    
    my $fragsTableTMP = newTable('Frags', "$sample\_TMP");     
    my $fragsFileTMP = "$fragsTableTMP.csv";   
    open my $fragsFileTMPH, ">", $fragsFileTMP; 
    runSQL("SELECT * 
            FROM $fragsTable
            WHERE FRAGMENTTYPE = $types{Frags}{SingleRead}
            ORDER BY PAIRID"); 
    my @srRefs;
    my $prevPairID = 0;      
    while (my $__ = fetchRowHashRef()){     
        if ($$__{PAIRID} > $prevPairID and $prevPairID){
             processSingleReads($fragsFileTMPH, \@srRefs);           
            @srRefs = ();
        }                     
        push @srRefs, $__;
        $prevPairID = $$__{PAIRID};
    } 
    processSingleReads($fragsFileTMPH, \@srRefs);  
    close $fragsFileTMPH;
    loadData($fragsFileTMP, $fragsTableTMP, ",", "FRAGMENTID, FRAGMENTTYPE,
            CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1,
            CHROMOSOME2, POSITION2, LENGTH2, STRAND2, DISCREPANCIES2,
            FRAGMENTSIZE, PAIRID, NFRAGS, 
            EVENTSIZE, STDEVNORMAL, ENDTOLERANCE, 
            NSETSFRAG, NSETSPAIR");
    runSQL("DELETE FROM $fragsTable
            WHERE FRAGMENTTYPE = $types{Frags}{SingleRead} ");  
    runSQL("INSERT INTO $fragsTable SELECT * FROM $fragsTableTMP");      
    dropTable($fragsTableTMP);    
} 

sub processSingleReads{
    my ($fragsFileH, $srRefsRef) = @_;
    scalar @$srRefsRef == 1 and return printBestReadsSR($fragsFileH, $srRefsRef);
    $srRefsRef = stratifyReadsSR(\&getNDisc, $srRefsRef);
    scalar @$srRefsRef == 1 and return printBestReadsSR($fragsFileH, $srRefsRef);
    $srRefsRef = stratifyReadsSR(\&getNInternalDiscs, $srRefsRef);     
    scalar @$srRefsRef == 1 and return printBestReadsSR($fragsFileH, $srRefsRef);
    #single reads which still map to more than genome location are dropped in current usage  
}

sub printBestReadsSR{
    my ($fragsFileH, $srRefsRef) = @_;
    foreach my $srRef(@$srRefsRef){ commitRow($fragsFileH, $srRef) }  
}

sub stratifyReadsSR{
    my ($nDiscsSubRef, $srRefsRef) = @_;
    my %strata; 
    foreach my $srRef(@$srRefsRef){ 
        my $nDiscs = &$nDiscsSubRef($$srRef{DISCREPANCIES1}, $$srRef{LENGTH1});
        push @{$strata{$nDiscs}}, $srRef;
    }
    my @strata = sort {$a <=> $b} keys %strata;
    my $minDiscs = $strata[0]; 
    return $strata{$minDiscs};
}

sub purgeDuplicateFragments{
    my ($inputFragsTable) = @_;
    $inputFragsTable or $inputFragsTable = $fragsTable;
    #The GROUP BY in this SQL purges duplicate Fragments.
    #Since all Pairs were re-ordered so that Position1 < Position2
    #this also means that reverse-orientation duplicates are also purged.
    #Unlike Pair purging during mapReads, this purging will purge mismatched Fragments
    #as long as the map positions and orientations match.
    #The Fragment with the least number of discrepancies is kept.
    status("  purging any remaining duplicate Fragments\n");
    my $sql = "SELECT Min(FRAGMENTID) FRAGMENTID, FRAGMENTTYPE,
                  CHROMOSOME1, POSITION1, LENGTH1, STRAND1, Min(DISCREPANCIES1) DISCREPANCIES1,
                  CHROMOSOME2, POSITION2, LENGTH2, STRAND2, Min(DISCREPANCIES2) DISCREPANCIES2,
                  FRAGMENTSIZE, Min(PAIRID) PAIRID, Max(NFRAGS) NFRAGS,
                  EVENTSIZE, STDEVNORMAL, ENDTOLERANCE, 
                  1 AS NSETSFRAG, 0 AS NSETSPAIR
                FROM $inputFragsTable
                GROUP BY  FRAGMENTTYPE,
                  CHROMOSOME1, POSITION1, LENGTH1, STRAND1,
                  CHROMOSOME2, POSITION2, LENGTH2, STRAND2,
                  FRAGMENTSIZE, EVENTSIZE, STDEVNORMAL, ENDTOLERANCE";
    updateTable('Frags', $inputFragsTable, $sql);
}

sub countNormals{
    runSQL("SELECT Count(FRAGMENTID)
            FROM $fragsTable
            WHERE FRAGMENTTYPE=$types{Frags}{Normal}");
    my ($normalCount) = fetchRowArray();
    updateStat($statsTable, 'normalCount', $normalCount);    
    status("  $normalCount Normal Fragments\n");    
    if ($param{isCircles}){
        runSQL("SELECT Count(FRAGMENTID)
                FROM $fragsTable
                WHERE FRAGMENTTYPE=$types{Frags}{ReverseNormal}");
        my ($reverseNormalCount) = fetchRowArray();
        updateStat($statsTable, 'reverseNormalCount', $reverseNormalCount);
        status("  $reverseNormalCount ReverseNormal Fragments\n");
        my $percentReverse = int(($reverseNormalCount/($reverseNormalCount + $normalCount))*100);
        status("  $percentReverse% of library is reversed (i.e. bead sticking) artifact\n")    
    }
}

sub createFragSizeHistogram{
    my ($histogramTable, $fragsTable) = @_;
    runSQL("INSERT INTO $histogramTable
            (SELECT FRAGMENTTYPE SERIES, FRAGMENTSIZE X, Count(FRAGMENTID) Y
                FROM $fragsTable
                GROUP BY FRAGMENTTYPE, FRAGMENTSIZE)  ");
}

sub createMaps {
    my ($sample, $fragsTable, $statsTable) = @_;  #pass these since also called by mergeFragments
    my $fragMapTable = newTable('FMap', $sample);
    createFragMap($fragsTable, $fragMapTable);
    my $histogramTable = newTable('h', $fragMapTable);
    createMapHistogram($histogramTable, $fragMapTable);
    calculateCoverageStats($histogramTable, 'FMap', $statsTable);
    my $readMapTable;
    unless ($param{noDisc}){
        $readMapTable = newTable('RMap', $sample);
        createReadMap($fragsTable, $readMapTable);
        $histogramTable = newTable('h', $readMapTable);
        createMapHistogram($histogramTable, $readMapTable);
        calculateCoverageStats($histogramTable, 'RMap', $statsTable);
    }
    return ($fragMapTable, $readMapTable);
}

sub createFragMap{
    my ($fragsTable, $fragMapTable) = @_;
    foreach my $chrom (1..$refSeqs{$param{refSeqBase}}{nChrom}){
        my %coverage = ();
        runSQL("SELECT ((trunc(POSITION1/$param{binSize})*$param{binSize})+$param{binSize}) AS LOWBIN,
                       (trunc((POSITION2 + LENGTH2)/$param{binSize})*$param{binSize}) AS HIGHBIN
                FROM $fragsTable
                WHERE CHROMOSOME1 = $chrom
                AND (FRAGMENTTYPE=$types{Frags}{Normal}
                     OR FRAGMENTTYPE=$types{Frags}{ReverseNormal})",
                \my($lowBin, $highBin));
        while (fetchRow()){
            for (my $bin = $lowBin; $bin<= $highBin; $bin += $param{binSize}){
                $coverage{$bin}++;
            }
        }
        my $fragMapFile = "$fragMapTable.csv";  #load the results into map table
        open my $fragMapFileH, ">", $fragMapFile;
        foreach my $bin (keys %coverage){
            print $fragMapFileH join(",", ($chrom, $bin, $coverage{$bin}))."\n";
        }
        close $fragMapFileH;
        loadData($fragMapFile, $fragMapTable, ",", "CHROMOSOME, POSITION, COVERAGE");
    }
}

sub createReadMap{
    my ($fragsTable, $readMapTable) = @_;
    #this will calculate read coverages at the same spacing as frag coverage
    #adjust binSize if tighter read coverage map is desired
    foreach my $chrom (1..$refSeqs{$param{refSeqBase}}{nChrom}){
        my %coverage = ();
        runSQL("SELECT ((trunc(POSITION1/$param{binSize})*$param{binSize})+$param{binSize}) AS LOWBIN1,
                        (trunc((POSITION1 + LENGTH1 - 1)/$param{binSize})*$param{binSize}) AS HIGHBIN1,
                       ((trunc(POSITION2/$param{binSize})*$param{binSize})+$param{binSize}) AS LOWBIN2,
                        (trunc((POSITION2 + LENGTH2 - 1)/$param{binSize})*$param{binSize}) AS HIGHBIN2,
                        FRAGMENTTYPE
                FROM $fragsTable
                WHERE CHROMOSOME1 = $chrom
                AND (FRAGMENTTYPE = $types{Frags}{Normal}
                     OR FRAGMENTTYPE = $types{Frags}{ReverseNormal})
                     OR FRAGMENTTYPE = $types{Frags}{SingleRead} ",
                \my($lowBin1, $highBin1, $lowBin2, $highBin2, $fragType));          
        while (fetchRow()){
            fillRead(\%coverage, $lowBin1, $highBin1);
            $fragType == $types{Frags}{SingleRead} or fillRead(\%coverage, $lowBin2, $highBin2);
        }
        my $readMapFile = "$readMapTable.csv";
        open my $readMapFileH, ">", $readMapFile;
        foreach my $pos (keys %coverage){
            print $readMapFileH join(",", ($chrom, $pos, $coverage{$pos}))."\n";
        }
        close $readMapFileH;
        loadData($readMapFile, $readMapTable, ",", "CHROMOSOME, POSITION, COVERAGE");
    }
}

sub fillRead{
    my ($coverageRef, $lowBin, $highBin) = @_;
    for (my $bin = $lowBin; $bin<= $highBin; $bin += $param{binSize}){
        $$coverageRef{$bin}++;
    }
}

sub createMapHistogram{
    my ($histogramTable, $mapTable) = @_;
    runSQL("INSERT INTO $histogramTable
            (SELECT 0 SERIES, COVERAGE X, Count(COVERAGE) Y
                FROM $mapTable
                GROUP BY COVERAGE)  ");
}

sub calculateCoverages{
    my ($fragsTable, $statsTable) = @_;
    runSQL("SELECT Sum(CHROMLENGTH) GENOMESIZE
                FROM
                (SELECT CHROMOSOME1, Max(POSITION1) CHROMLENGTH
                FROM $fragsTable
                WHERE CHROMOSOME1 <= $refSeqs{$param{refSeqBase}}{nChrom}
                GROUP BY CHROMOSOME1)");
    my ($genomeSize) = fetchRowArray();
    runSQL("SELECT Sum(FRAGMENTSIZE)
            FROM $fragsTable
            WHERE (FRAGMENTTYPE=$types{Frags}{Normal}
                   OR FRAGMENTTYPE=$types{Frags}{ReverseNormal})");
    my ($fragCoverage) = fetchRowArray();
    $fragCoverage /= $genomeSize;
    runSQL("SELECT Sum(LENGTH1 + LENGTH2)
            FROM $fragsTable
            WHERE (FRAGMENTTYPE=$types{Frags}{Normal}
                   OR FRAGMENTTYPE=$types{Frags}{ReverseNormal})");
    my ($readCoverage) = fetchRowArray();
    $readCoverage /= $genomeSize;
    updateStat($statsTable, 'fragCoverage', $fragCoverage);
    updateStat($statsTable, 'readCoverage', $readCoverage);
    status("  Fragment coverage = $fragCoverage\n");
    status("  Read coverage = $readCoverage\n");
    return ($fragCoverage, $readCoverage);
}

sub calculateNormalizedCoverage{
    my ($mapTable, $mapType, $coverageMean) = @_;
    #calculate normalized coverage per bin, used during visualization
    updateTable($mapType, $mapTable,
                "SELECT CHROMOSOME, POSITION, COVERAGE,
                    (COVERAGE/$coverageMean) AS NORMALIZEDCOVERAGE
                FROM $mapTable
                GROUP BY CHROMOSOME, POSITION, COVERAGE
                ORDER BY CHROMOSOME, POSITION");
}

1;



