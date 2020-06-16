#!/usr/bin/perl -w
use strict;
use warnings;

#########################################################################
#FindSet.pl finds Sets of overlapping Fragments whose ends are within
#distances consistent with them being different real fragments crossing
#the same anomalous junction.
#Information about such sets of fragments is stored in the Sets table, and
#the Fragments table is updated to include counts of Sets per Fragment and Pair.
#########################################################################

use vars(qw(%param %types %fields %fieldNames %refSeqs $uid));

sub findSets1{ #used by find, takes only one sample as the comparison input
    ($param{setType}) = @_;
    my @strands = initializeFindSets();
    status("  Set type $param{setType}:\n    Chr: ");
    foreach my $chrom1 (1..$refSeqs{$param{refSeqBase}}{nChrom}){
        $param{chrom1} = $chrom1;
        status("$chrom1 ");        
        foreach my $strandsRef(@strands){
            $param{strand1} = $$strandsRef[0];
            $param{strand2} = $$strandsRef[1];
            my $sql1 = findSetsSQL(1, $param{fragsTable1});
            threadP1(" $sql1 ", \&threadP2);
        }
    }
    status("\n");
}

sub findSets2{ #used by compare, takes two samples as the comparison input
    ($param{setType}) = @_;
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
            threadP1(" ($sql1 UNION ALL $sql2) ", \&threadP2);
        }
    }
    status("\n");
}

sub initializeFindSets{
    $param{chrom2} = 0;
    $param{enforceOverlap} = 1;
    my %strands = ( #allowable strand combinations for each Frag/Set type
        $types{Sets}{Deletion} => [[1,2]],
        $types{Sets}{Insertion} => [[1,2]],
        $types{Sets}{Inversion} => [[1,1], [2,2]],
        $types{Sets}{Duplication} => [[2,1]],
        $types{Sets}{Rooted} => [[1,1], [2,2]],
    );
    return @{$strands{$param{setType}}};
}

sub findSetsSQL{
    my ($sampleN, $fragsTable) = @_;  
    return "SELECT FRAGMENTID, PAIRID, FRAGMENTTYPE,
                    CHROMOSOME1, POSITION1, POSITION2, FRAGMENTSIZE, 
                    EVENTSIZE, STDEVNORMAL, ENDTOLERANCE, $sampleN AS SAMPLE
            FROM $fragsTable
            WHERE FRAGMENTTYPE=$param{setType}
                AND FRAGMENTID > 0
                AND CHROMOSOME1=$param{chrom1}
                AND STRAND1=$param{strand1} 
                AND NHITS1 <= $param{maxHits} AND NHITS2 <= $param{maxHits} 
                $param{strataFilter}"; #on 2nd find iteration, only consider "best" frags flagged on first iteration
                #FRAGMENTID > 0 accounts for Fragment masking by reduceNoise
}

#####################################################################################
#This block defines the core logic for set finding.  The overall logic is to:
#1) order all fragments by position1
#2) step through until a "break" occurs that identifies the current fragment
#    as NOT being in the current position1 group (i.e. P1s violate max allowed as-sequenced overlap)
#3) reorder the fragments in the current position1 group by position2
#4) step through until a position2 break occurs
#5) split position2 group into as-mapped overlap groups (i.e. sets must be overlapping in as-mapped as well as as-sequenced views)
#6) re-examine the position1 and position2 groups to determine end-to-end span and therefore singularity
#7) for compare, separate fragments back into sample groups
#8) if an allowable number of fragments persist after all breaks, save them as a set
#9) continue looping back interatively until all sets and fragments are processed

sub threadP1{ #thread through all fragments of type in chromosome, looking for Position1 breaks between fragments
    my ($sql, $nextSubRef) = @_;
    #Primary sort by position1.  To account for reads with identical position1, also use 
    #secondary sort by end tolerance with order designed to be most permissive, i.e. least likely to break at this stage.
    #End tolerance sort depends on the strand of read1 since strand1 break test uses endtolerance of lower position1, and vice versa.   
    my $etSortOrder = " ASC ";
    $param{strand1} == 2 and $etSortOrder = " DESC ";    
    runSQL(" $sql ORDER BY POSITION1, ENDTOLERANCE $etSortOrder ");
    my $fragRef = fetchRowHashRef();
    $fragRef or return;    
    my @frags;    
    push @frags, $fragRef;               
    my $prevP1 = $$fragRef{POSITION1}; 
    my $prevET = $$fragRef{ENDTOLERANCE}; 
    while (my $fragRef = fetchRowHashRef()){ #at this stage, just check Position1 against _previous_ fragment
        if (breakValue($prevP1, $$fragRef{POSITION1}, $prevET, $$fragRef{ENDTOLERANCE}, $param{strand1})){
            &$nextSubRef(\@frags);
            @frags = ();
        }
        push @frags, $fragRef;        
        $prevP1 = $$fragRef{POSITION1};
        $prevET = $$fragRef{ENDTOLERANCE};
    }
    &$nextSubRef(\@frags);
}

sub threadP2{ #thread through all fragments in Position1 group, looking for Position2 breaks between fragments
    my ($fragsRef) = @_;
    my $nFrags = scalar(@$fragsRef);
    $nFrags >= $param{minFragsSet} or return; #not good to enforce max yet, P1 alone is not restrictive enough!
    $param{setSigs} = {}; #track sets to ensure uniquness of all committed sets
    $fragsRef = reorderFrags($fragsRef, 2);
    my @frags;     
    push @frags, $$fragsRef[0];      
    my $prevP2 = $$fragsRef[0]{POSITION2};
    my $prevET = $$fragsRef[0]{ENDTOLERANCE};
    foreach my $i(1..($nFrags - 1)){ #at this stage, just check Position2 against _previous_ fragment  
        if (breakValue($prevP2, $$fragsRef[$i]{POSITION2}, $prevET, $$fragsRef[$i]{ENDTOLERANCE}, $param{strand2})){
            processP2Group(\@frags);
            @frags = ();
        }    
        push @frags, $$fragsRef[$i];        
        $prevP2 = $$fragsRef[$i]{POSITION2};
        $prevET = $$fragsRef[$i]{ENDTOLERANCE};
    }
    processP2Group(\@frags);
}

sub reorderFrags{ #used to establish the proper sort (based on read and strand) of a candidate set for the next threading/checking step
    my ($fragsRef, $read) = @_;
    my ($posField, $strand) = ("POSITION$read", $param{"strand$read"});
    my %fragsByPos;
    my @frags;      
    foreach my $fragRef(@$fragsRef){ push @{$fragsByPos{$$fragRef{$posField}}}, $fragRef }      
    my @pos = sort {$a <=> $b} keys %fragsByPos; #primary sort by position
    foreach my $pos(@pos){      
        if(scalar @{$fragsByPos{$pos}} > 1){        
            my %fragsByET;          
            foreach my $fragRef(@{$fragsByPos{$pos}}){ push @{$fragsByET{$$fragRef{ENDTOLERANCE}}}, $fragRef }    
            my @ETs;   
            if ($strand == 1){  #secondary sort by endTolerance, order determined by strand
                @ETs = sort {$a <=> $b} keys %fragsByET; 
            } else {
                @ETs = sort {$b <=> $a} keys %fragsByET; 
            }      
            foreach my $ET(@ETs){
                push @frags, @{$fragsByET{$ET}}; 
            }      
        } else {
            push @frags, @{$fragsByPos{$pos}}
        }
    }
    return \@frags;
}

sub processP2Group{ #break P1/P2 groups into non-overlapping as-mapped groups of fragments (required for adjacent very small fragments)
    my ($fragsRef) = @_;
    my $nFrags = scalar(@$fragsRef);   
    ($nFrags >= $param{minFragsSet} and $nFrags <= (2 * $param{maxFragsSet})) or return;     
    my @overlapGroups;    
    if ($param{enforceOverlap}){ 
        #perform initial break into as-mapped overlap groups (as-sequenced overlap was tested by the end tolerance test performed above)
        my (@minP1s, @maxP2s);
        push @overlapGroups, []; 
        my $oGI = 0;  
        push @{$overlapGroups[$oGI]}, $$fragsRef[0];      
        @minP1s = ($$fragsRef[0]{POSITION1});
        @maxP2s = ($$fragsRef[0]{POSITION2});
        foreach my $i(1..($nFrags - 1)){ #fragments supplied by threadP2 in P2 sort order, P1s don't necessarily come in order   
            if ($$fragsRef[$i]{POSITION1} >= $maxP2s[$oGI]){
                push @overlapGroups, [];
                $oGI = scalar @overlapGroups - 1;
                push @minP1s, 1E9;     
                push @maxP2s, 0;
            }
            push @{$overlapGroups[$oGI]}, $$fragsRef[$i];           
            $minP1s[$oGI] <= $$fragsRef[$i]{POSITION1} or $minP1s[$oGI] = $$fragsRef[$i]{POSITION1};  
            $maxP2s[$oGI] = $$fragsRef[$i]{POSITION2};                                
        }
        #check to ensure that overlap groups themselves do not overlap (still possible given lack of P1 ordering in above loop)
        $oGI = 1; 
        while ($overlapGroups[$oGI]){
            if($minP1s[$oGI] <= $maxP2s[$oGI - 1]){
                push @{$overlapGroups[$oGI - 1]}, @{$overlapGroups[$oGI]}; #merge overlapping overlap groups          
                splice(@overlapGroups, $oGI, 1); #removes array element, shifts all remaining elements down one      
            } else {
                $oGI++;
            } 
        }          
    } else { #as-mapped fragment overlap has no meaning for DiffChrom Sets
        @overlapGroups = ($fragsRef);
    }
    #finally check overlap groups for self consistency, breaking into smaller groups if needed
    foreach my $oGroupRef (@overlapGroups){ 
        my $p1SetsRef = checkSetPositions($oGroupRef, 1);
        foreach my $p1SetRef(@$p1SetsRef){
            my $p2SetsRef = checkSetPositions($p1SetRef, 2);
            foreach my $p2SetRef(@$p2SetsRef){
                processSet($p2SetRef)
            }
        }
    }    
}

sub checkSetPositions{ #check if the fragments remaining in a group are legitimately contained within a single set
    my ($fragsRef, $read) = @_;
    my @sets;     
    my $nFrags = scalar(@$fragsRef);
    if ($nFrags >= $param{minFragsSet}){  
        $fragsRef = reorderFrags($fragsRef, $read);      
        my ($posField, $strand) = ("POSITION$read", $param{"strand$read"});
        if (breakValue($$fragsRef[0]{$posField},    $$fragsRef[$nFrags - 1]{$posField}, 
                       $$fragsRef[0]{ENDTOLERANCE}, $$fragsRef[$nFrags - 1]{ENDTOLERANCE}, $strand)){
            breakSet(\@sets, $fragsRef, $posField, $strand);
        } else {
            push @sets, $fragsRef;
        }        
    }
    return \@sets;
}

sub breakSet{ #if P1 or P2's are too widely separated, break into multiple, potentially overlapping, sets
    my ($setsRef, $fragsRef, $posField, $strand) = @_;
    my $maxI = scalar(@$fragsRef) - 1; 
    my $iOffset = $param{minFragsSet} - 1;
    my %goodSpans;
    #find all non-breaking fragment groups within the collection of fragments
    foreach my $lowI(0..($maxI - $iOffset)){
        my $lastI;
        foreach my $highI(($lowI + $iOffset)..$maxI){
            if (breakValue($$fragsRef[$lowI]{$posField},    $$fragsRef[$highI]{$posField}, 
                           $$fragsRef[$lowI]{ENDTOLERANCE}, $$fragsRef[$highI]{ENDTOLERANCE}, $strand)){
                last;            
            } else {
                $lastI = $highI;
            }
        }
        defined $lastI and $goodSpans{"$lowI:$lastI"}++;
    }
    #mask non-breaking fragment groups entirely subsumed by other fragment groups
    my @goodSpans = keys %goodSpans;
    foreach my $testSpan(@goodSpans){
        my $keepTestSpan = 1;
        my ($testLow, $testHigh) = split(":", $testSpan);
        foreach my $refSpan(@goodSpans){
            my ($refLow, $refHigh) = split(":", $refSpan);
            if(($testLow > $refLow and $testHigh <= $refHigh) or ($testLow >= $refLow and $testHigh < $refHigh)) {
                $keepTestSpan = 0;
                last;
            }
        }
        $goodSpans{$testSpan} = $keepTestSpan;
    }
    #save the final distinct groups
    foreach my $span(@goodSpans){
        $goodSpans{$span} or next;
        my ($lowI, $highI) = split(":", $span); 
        push @$setsRef, [@$fragsRef[$lowI..$highI]]; 
    }
}

###########################################
##these commented checkSetPositions  and breakSet carry important updates to vamp 2.1
##mostly corrected a bug that was causing all breaking uber groups to always be called as a set
##If don't like the consequences of the (probably more correct) version above
##would likely want to revert to using the ones below rather than 2.1
##Not sure which is faster, probably this one.
###########################################
##sub checkSetPositions{ #check if the fragments remaining in a group are legitimately contained within a single set
##    my ($setRef, $position) = @_;
##    my @sets;     
##    if (scalar(@$setRef) >= $param{minFragsSet}){    
##        my $pField = "POSITION$position"; 
##        my $strand = $param{"strand$position"};
##        my (%fragsByP, %endTolerances);                
##        foreach my $fragRef(@$setRef){
##            push @{$fragsByP{$$fragRef{$pField}}}, $fragRef;
##            ($endTolerances{$$fragRef{$pField}} and 
##                $endTolerances{$$fragRef{$pField}} >= $$fragRef{ENDTOLERANCE}) or
##                    $endTolerances{$$fragRef{$pField}} = $$fragRef{ENDTOLERANCE};
##        }
##        my @ps = sort {$a <=> $b} keys %fragsByP; 
##        if(breakIndex(\%endTolerances, \@ps, 0, scalar(@ps) - 1, $strand)){ #here, check highest and lowest positions in set, i.e. the set span
##            while (scalar @ps){ breakSet(\@sets, \@ps, \%fragsByP, \%endTolerances, $strand) }          
##        } else {
##            push @sets, $setRef;
##        }
##    }
##    return \@sets;
##}

##sub breakSet{ #if P1 or P2's are too widely separated, break into multiple, potentially overlapping, sets
##    my ($setsRef, $psRef, $fragsByPRef, $endTolsRef, $strand) = @_;
##    my $maxPsI = scalar(@$psRef) - 1;  
##    if (breakIndex($endTolsRef, $psRef, 0, $maxPsI, $strand)){
##        my (@lowSet, @highSet, $lowI, $highI);
##        my ($lowSetIndex, $highSetIndex, $lowSetIs, $highSetIs) = (0, $maxPsI, '', ''); 
##        while (scalar(@lowSet) < $param{minFragsSet} and $lowSetIndex < $maxPsI){
##            @lowSet = ();
##            $highI = $maxPsI; 
##            while (breakIndex($endTolsRef, $psRef, $lowSetIndex, $highI, $strand)){$highI--}  
##            foreach my $pI ($lowSetIndex..$highI){ push @lowSet, @{$$fragsByPRef{$$psRef[$pI]}} }  
##            $lowSetIndex++;
##        }
##        $lowSetIndex--;
##        $lowSetIs = "$lowSetIndex:$highI";
##        while (scalar(@highSet) < $param{minFragsSet} and $highSetIndex > 0){
##            @highSet = ();
##            $lowI = 0;
##            while (breakIndex($endTolsRef, $psRef, $lowI, $highSetIndex, $strand)){$lowI++}
##            foreach my $pI ($lowI..$highSetIndex){ push @highSet, @{$$fragsByPRef{$$psRef[$pI]}} }
##            $highSetIndex--;  
##        }
##        $highSetIndex++;
##        $highSetIs = "$lowI:$highSetIndex";
##        push @$setsRef, \@lowSet;         
##        !($lowSetIs eq $highSetIs) and push @$setsRef, \@highSet;
##        my $remaining = $lowI - $highI - 1;
##        if ($remaining > 0){
##            @$psRef  = splice(@$psRef, $highI + 1, $remaining);     
##        } else {
##            @$psRef = ();
##        }       
##    } else {
##        my @passSet;
##        foreach my $pI(0..$maxPsI){ push @passSet, @{$$fragsByPRef{$$psRef[$pI]}} } 
##        push @$setsRef, \@passSet;
##        @$psRef = ();
##    }
##}

sub breakValue{
    my ($lowPos, $highPos, $lowET, $highET, $strand) = @_;
    my $endTolerance = $strand == 1 && $lowET || $highET;
    return ($highPos - $lowPos) > $endTolerance;
}

#sub breakIndex{
#    my ($endTolsRef, $psRef, $lowPsI, $highPsI, $strand) = @_;
#    my $endTolerance = $strand == 1 && $$endTolsRef{$$psRef[$lowPsI]} || $$endTolsRef{$$psRef[$highPsI]};
#    return ($$psRef[$highPsI] - $$psRef[$lowPsI]) > $endTolerance;
#}

sub processSet{ #separate the set into its samples, check one last time for min and max nFrags
    my ($fragsRef) = @_;
    (scalar(@$fragsRef) >= $param{minFragsSet}) or return;
    my (@set1, @set2);
    foreach my $fragRef(@$fragsRef){
        if ($$fragRef{SAMPLE} == 1){
            push @set1, $fragRef;
        } else {
            push @set2, $fragRef; 
        }          
    }
    my $nFragsSet1 = scalar @set1; 
    my $nFragsSet2 = scalar @set2;
    my $set1Pass = ($nFragsSet1 >= $param{minFragsSet1} and $nFragsSet1 <= $param{maxFragsSet1});
    my $set2Pass = ($nFragsSet2 >= $param{minFragsSet2} and $nFragsSet2 <= $param{maxFragsSet2});
    ($set1Pass or $set2Pass) and $param{setID}++;  
    if (($set1Pass and $nFragsSet2) or ($set2Pass and $nFragsSet1)){ #if one set is good, keep them both regardless
        my $nFragsTotal = $nFragsSet1 + $nFragsSet2;
        printSet(\@set1, $nFragsSet1, $nFragsTotal, 1);
        printSet(\@set2, $nFragsSet2, $nFragsTotal, 2);
    } elsif ($set1Pass) {
        printSet(\@set1, $nFragsSet1, $nFragsSet1, 1);
    } elsif ($set2Pass) {
        printSet(\@set2, $nFragsSet2, $nFragsSet2, 2);
    } 
}

sub printSet{ #commit set to files for uploading to db later on
    my ($setRef, $nFragsSet, $nFragsTotal, $sampleN) = @_;
    my @parsedSet = getSet($setRef, $nFragsSet);
    my $setSig = join(":", ($sampleN, @parsedSet[0..4])); 
    unless($param{setSigs}{$setSig}){ #prevent set duplicates
        my ($setsFileH, $fIDsFileH, $pIDsFileH) = getFSFileHs($sampleN);
        $param{strataFilter} and 
            print $setsFileH join(",", ($param{setID}, $param{setType}, $param{chrom1}, $param{chrom2},
                                        @parsedSet, $param{strand1}, $param{strand2},
                                        $nFragsTotal, $nFragsSet))."\n"; 
        foreach my $fragRef (@$setRef){
            print $fIDsFileH "$param{setID},$$fragRef{FRAGMENTID}\n";
            $param{strataFilter} and print $pIDsFileH "$param{setID},$$fragRef{PAIRID}\n";
        }        
    }
    $param{setSigs}{$setSig} = 1;
}

sub getSet{  #determine span, overlap region, and mean and stDev of size for Fragments in a Set
    my ($setRef, $nFragsSet) = @_;
    my ($spanStart, $spanEnd, $overlapStart, $overlapEnd, $fragSizeSum, $eventSizeSum) = calculateSetLimits($setRef, $nFragsSet);
    $param{strataFilter} or return ($spanStart, $spanEnd, $overlapStart, $overlapEnd, $fragSizeSum);
    my ($fragmentMean, $eventMean, $stDev) = calculateSetSize($setRef, $nFragsSet, $fragSizeSum, $eventSizeSum);
    my ($meanDelta, $fractionDelta2, $fractionOverlapLeft, $fractionOverlapRight) = calculateSetDelta($setRef, $nFragsSet, $spanStart, $spanEnd, $overlapStart, $overlapEnd);
    return ($spanStart, $spanEnd, $overlapStart, $overlapEnd, $fragmentMean, $eventMean, $stDev, 
            $meanDelta, $fractionDelta2, $fractionOverlapLeft, $fractionOverlapRight);
}

sub calculateSetLimits{  #determine set span and overlap limits
    my ($setRef, $nFragsSet) = @_;
    my ($spanStart, $spanEnd, $overlapStart, $overlapEnd, $fragSizeSum, $eventSizeSum) = (0, 0, 0, 0, 0, 0);
    foreach my $fragRef (@$setRef){
        my %frag = %$fragRef;        
        if (!$spanStart or ($frag{'POSITION1'}<$spanStart)){$spanStart=$frag{'POSITION1'}} #smallest P1
        if (!$spanEnd or ($frag{'POSITION2'}>$spanEnd)){$spanEnd=$frag{'POSITION2'}} #largest P2
        if (!$overlapStart or ($frag{'POSITION1'}>$overlapStart)){$overlapStart=$frag{'POSITION1'}} #largest P1
        if (!$overlapEnd or ($frag{'POSITION2'}<$overlapEnd)){$overlapEnd=$frag{'POSITION2'}} #smallest P2
        $fragSizeSum += $frag{'FRAGMENTSIZE'};
        $eventSizeSum += $frag{'EVENTSIZE'};
    }
    return ($spanStart, $spanEnd, $overlapStart, $overlapEnd, $fragSizeSum, $eventSizeSum);
}

sub calculateSetSize{  #determine event size and stdev of size
    my ($setRef, $nFragsSet, $fragSizeSum, $eventSizeSum) = @_;
    my ($fragmentMean, $eventMean, $stDev) = (0, 0, 0);
    unless ($param{setType} == $types{Sets}{DiffChrom}) { #size not meaningful for DiffChrom
        $fragmentMean = int($fragSizeSum/$nFragsSet);
        $eventMean = int($eventSizeSum/$nFragsSet);
        if ($nFragsSet > 1){
            my $dev = 0;
            my $dev2 = 0;
            my $sumDev2 = 0;
            foreach my $fragRef (@$setRef){
                my %frag = %$fragRef; 
                $dev = $frag{'FRAGMENTSIZE'} - $fragmentMean;
                $dev2 = $dev**2;
                $sumDev2 += $dev2;
            }
            $stDev = int(sqrt ($sumDev2/($nFragsSet -1)));
        }    
    }
    return ($fragmentMean, $eventMean, $stDev);
}

sub calculateSetDelta{  #infer physical fragment size (i.e. as sequenced) and delta of each fragment from mode, expressed as # of SDs
    my ($setRef, $nFragsSet, $spanStart, $spanEnd, $overlapStart, $overlapEnd) = @_;
    my ($leftRef, $rightRef) =  ($overlapStart, $overlapEnd);
    $param{strand1} == 2 and $leftRef = $spanStart; 
    $param{strand2} == 1 and $rightRef = $spanEnd;    
    my $overlapSum;
    foreach my $fragRef (@$setRef){ 
        $$fragRef{modeNormal} = $$fragRef{FRAGMENTSIZE} + $$fragRef{EVENTSIZE};      
        $$fragRef{excess} = abs($leftRef - $$fragRef{POSITION1}) +  abs($$fragRef{POSITION2} - $rightRef);
        $overlapSum += $$fragRef{modeNormal} - $$fragRef{excess};             
    }      
    my $overlap = int(($overlapSum/$nFragsSet) + 0.5); 
    my $overlapTMP = $overlap;
    $overlapTMP < 1 and $overlapTMP = 1; #overlap can be negative! avoid the confusion this creates by disallowing
    my $excessLeft =  abs($overlapStart - $spanStart) * 2;
    my $excessRight = abs($spanEnd - $overlapEnd) * 2; 
    my $fractionOverlapLeft =  1 - ($excessLeft /  ($overlapTMP + $excessLeft)); 
    my $fractionOverlapRight = 1 - ($excessRight / ($overlapTMP + $excessRight)); 
    my ($deltaSum, $delta2Count, $withNormalSum) = (0, 0, 0);
    foreach my $fragRef(@$setRef){ 
        my $delta = abs(($$fragRef{excess} + $overlap) - $$fragRef{modeNormal}) / $$fragRef{STDEVNORMAL}; 
        $delta > 2 and $delta2Count++;   
        $deltaSum += $delta;        
        $param{rooted} and $withNormalSum += abs($$fragRef{FRAGMENTSIZE});
    } 
    my $meanDelta = $deltaSum / $nFragsSet;
    my $fractionDelta2 = $delta2Count / $nFragsSet;   
    $param{rooted} and $fractionDelta2 = $withNormalSum / $nFragsSet;   #FRACTIONDELTA2 holds FRACTIONWITHNORMAL for Rooted Sets
    return ($meanDelta, $fractionDelta2, $fractionOverlapLeft, $fractionOverlapRight);
}

sub establishFragmentStrata{ #on 1st find iteration, pick "best" frag from among promiscuous sister frags based on lowest nDisc
    my ($fragsTable, $sampleN) = @_;
    my ($setsFile, $fIDsFile, $pIDsFile) =  getFSFiles($sampleN);
    my ($setsTable, $fIDsTable, $pIDsTable) =  getFSTables($sampleN);
    loadIDs($fIDsFile, $fIDsTable);  
    runSQL("SELECT f.FRAGMENTID, f.PAIRID, 
                   f.DISCREPANCIES1 || 'x' DISCS1, f.DISCREPANCIES2 || 'x' DISCS2,
                   f.LENGTH1, f.LENGTH2
            FROM $fragsTable f, $fIDsTable fis
            WHERE f.FRAGMENTID = fis.ID_",
           \my($fragID, $pairID, $discs1, $discs2, $length1, $length2) ); #gather discrepancy information on fragments in sets (fis)    
    my %strata;
    while (fetchRow()){ 
        my $nDiscs = getNInternalDiscs($discs1, $length1) + getNInternalDiscs($discs2, $length2);
        push @{$strata{$pairID}{$nDiscs}}, $fragID;
    }
    my $priorityFIDsTable = newTable('IDs', $fIDsTable);
    my $priorityFIDsFile = "$priorityFIDsTable.csv";
    open my $priorityFIDsFileH, ">", $priorityFIDsFile;      
    foreach my $pairID(keys %strata){
        my ($lowestStrata) = sort {$a <=> $b} keys %{$strata{$pairID}};    
        foreach my $fragID(@{$strata{$pairID}{$lowestStrata}}){ print $priorityFIDsFileH "0,$fragID\n" }
    }    
    close $priorityFIDsFileH;  
    loadIDs($priorityFIDsFile, $priorityFIDsTable);
    runSQL("UPDATE $fragsTable 
            SET NSETSFRAG = -1 
            WHERE FRAGMENTID IN 
            (SELECT ID_ FRAGMENTID FROM $priorityFIDsTable)"); #flag best frags by temporarily setting NSETSFRAG to -1
    dropTable($priorityFIDsTable);
}

sub getNInternalDiscs{  #rank frags among promiscusous sisters by counting discs, ignoring unreliabe discs at the ends of reads
    my ($discs, $length) = @_;
    $discs =~ m/(.*)x$/; #strip the 'x' appended to force Perl to use $discs as string    
    $1 or return 0;
    length($1) < 3 and return $1;
    my @discPos;    
    while ($1){ 
        #$1 = discreps remaining; $2 = relPos, $3 = discType, $4 = extra    
        $1 =~ m/(.*)(..)(.)(.)$/ or $1 =~ m/(.*)(.)(.)(.)$/;   
        $3 and push @discPos, $2;
    }   
    @discPos = sort {$a <=> $b} @discPos; 
    my $exclude = 1;
    while (scalar @discPos and $discPos[0] == $exclude){shift @discPos; $exclude++}
    @discPos = sort {$b <=> $a} @discPos; 
    $exclude = $length;    
    while (scalar @discPos and $discPos[0] == $exclude){shift @discPos; $exclude--}   
    return scalar @discPos;
}

sub loadFSData{
    my ($sampleN) = @_;
    my ($setsFile, $fIDsFile, $pIDsFile) =  getFSFiles($sampleN);
    my ($setsTable, $fIDsTable, $pIDsTable) =  getFSTables($sampleN);
    loadSets($setsFile, $setsTable);
    loadIDs($fIDsFile, $fIDsTable);
    loadIDs($pIDsFile, $pIDsTable);
}

sub loadSets{
    my ($setsFile, $setsTable) = @_;   
    loadData($setsFile, $setsTable, ",", "SETID, SETTYPE, CHROMOSOME1, CHROMOSOME2,
                                                SPANSTART, SPANEND, OVERLAPSTART, OVERLAPEND,
                                                FRAGMENTMEAN, EVENTMEAN, STDEV, 
                                                MEANDELTA, FRACTIONDELTA2, FRACTIONOVERLAPLEFT, FRACTIONOVERLAPRIGHT,
                                                STRAND1, STRAND2,
                                                NFRAGSTOTAL, NFRAGSSAMPLE");
}

sub loadIDs{
    my ($IDsFile, $IDsTable) = @_;
    loadData($IDsFile, $IDsTable, ",", "SETID, ID_");        
}

sub calculateFracBadFrags{
    my ($sampleN) = @_;
    my ($setsTable, $fIDsTable, $pIDsTable) =  getFSTables($sampleN);
    my $pairSetCountSQL = "SELECT ID_, Count(SETID) NSETSPAIR
                            FROM $pIDsTable
                            GROUP BY ID_";  #how many Sets a Pair contributes to
    my $badPairsSQL = "SELECT ID_
                        FROM ($pairSetCountSQL)
                        WHERE NSETSPAIR > 1";  #the list of Pairs which are promiscuous, i.e. in more than 1 Set
    my $badFragCountSQL = "SELECT SETID, Count(SETID) NBADFRAGS
                    FROM $pIDsTable pid, ($badPairsSQL) bp
                    WHERE pid.ID_ = bp.ID_
                    GROUP BY SETID"; #the number of Fragments in each Set which arise from promiscuous Pairs
    my $fracBadFragsSQL = "SELECT st.SETID, st.SETTYPE,
                    st.CHROMOSOME1, st.CHROMOSOME2,
                    st.SPANSTART, st.SPANEND, st.OVERLAPSTART, st.OVERLAPEND,
                    st.FRAGMENTMEAN, st.EVENTMEAN, st.STDEV, 
                    st.MEANDELTA, st.FRACTIONDELTA2, st.FRACTIONOVERLAPLEFT, st.FRACTIONOVERLAPRIGHT,
                    st.STRAND1, st.STRAND2,
                    st.NFRAGSTOTAL, st.NFRAGSSAMPLE,
                    (nvl(bf.NBADFRAGS,0)/st.NFRAGSSAMPLE) FRACBADFRAGS,
                    0 MARK, '' DESCRIPTION
                FROM $setsTable st, ($badFragCountSQL) bf
                WHERE st.SETID = bf.SETID(+)";   
    updateTable('Sets', $setsTable, $fracBadFragsSQL); 
}

sub fillSetCounts{
    my ($fragsTable, $sampleN) = @_;
    my ($setsTable, $fIDsTable, $pIDsTable) =  getFSTables($sampleN);
    my $fragSetCountSQL = "SELECT ID_ FRAGMENTID, Count(ID_) SETCOUNT
                            FROM $fIDsTable
                            GROUP BY ID_";
    my $pairSetCountSQL = "SELECT ID_ PAIRID, Count(ID_) SETCOUNT
                            FROM $pIDsTable
                            GROUP BY ID_";    
    my $sql = "SELECT f.FRAGMENTID, f.FRAGMENTTYPE,
                    f.CHROMOSOME1, f.POSITION1, f.LENGTH1, f.STRAND1, f.DISCREPANCIES1, f.NHITS1,
                    f.CHROMOSOME2, f.POSITION2, f.LENGTH2, f.STRAND2, f.DISCREPANCIES2, f.NHITS2,
                    f.FRAGMENTSIZE, f.PAIRID, f.NFRAGS,
                    f.EVENTSIZE, f.STDEVNORMAL, f.ENDTOLERANCE, 
                    f.NSETSFRAG, nvl(psc.SETCOUNT,0) NSETSPAIR
                FROM
                (SELECT f.FRAGMENTID, f.FRAGMENTTYPE,
                    f.CHROMOSOME1, f.POSITION1, f.LENGTH1, f.STRAND1, f.DISCREPANCIES1, f.NHITS1,
                    f.CHROMOSOME2, f.POSITION2, f.LENGTH2, f.STRAND2, f.DISCREPANCIES2, f.NHITS2,
                    f.FRAGMENTSIZE, f.PAIRID, f.NFRAGS,
                    f.EVENTSIZE, f.STDEVNORMAL, f.ENDTOLERANCE,
                    nvl(fsc.SETCOUNT,0) NSETSFRAG, 0 NSETSPAIR
                FROM $fragsTable f, ($fragSetCountSQL) fsc
                WHERE f.FRAGMENTID = fsc.FRAGMENTID(+) ) f, ($pairSetCountSQL) psc
                WHERE f.PAIRID = psc.PAIRID(+) ";
    updateTable('Frags', $fragsTable, $sql );
    dropTable($fIDsTable);
    dropTable($pIDsTable);
}

sub fixSetLimits{
    my ($sample) = @_;
    my $setsTable = getTableName('Sets', $sample);   
    my $fragsTable = getTableName('Frags', $sample);    
    #step 1, calculate new span and overlap limits, store in temp field (SQL can't execute all in one step...)
    my $nonPromFragSQL = "SELECT * FROM $fragsTable WHERE NSETSFRAG = 1 AND NSETSPAIR = 1"; 
    my $spanStart = "Min(f.POSITION1)";
    my $spanEnd = "Max(f.POSITION2)";
    my $overlapStart = "Max(f.POSITION1)";
    my $overlapEnd = "Min(f.POSITION2)";    
    my $newLimitsSQL = "SELECT s.SETID, s.SETTYPE, s.CHROMOSOME1, s.CHROMOSOME2,
                                s.SPANSTART, s.SPANEND, s.OVERLAPSTART, s.OVERLAPEND,    
                                s.FRAGMENTMEAN, s.EVENTMEAN, s.STDEV, s.MEANDELTA, s.FRACTIONDELTA2,   
                                0 FRACTIONOVERLAPLEFT, 0 FRACTIONOVERLAPRIGHT,
                                s.STRAND1, s.STRAND2, s.NFRAGSTOTAL, s.NFRAGSSAMPLE, s.FRACBADFRAGS, s.MARK, s.DESCRIPTION,
                                $spanStart SPANSTART_NEW, $spanEnd SPANEND_NEW, 
                                $overlapStart OVERLAPSTART_NEW, $overlapEnd OVERLAPEND_NEW       
                         FROM ($setsTable) s, ($nonPromFragSQL) f
                         WHERE f.FRAGMENTTYPE = s.SETTYPE
                           AND f.CHROMOSOME1 = s.CHROMOSOME1 AND f.CHROMOSOME2 = s.CHROMOSOME2
                           AND f.POSITION1 >= s.SPANSTART AND f.POSITION1 <= s.OVERLAPSTART
                           AND f.POSITION2 >= s.OVERLAPEND AND f.POSITION2 <= s.SPANEND 
                         GROUP BY s.SETID, s.SETTYPE, s.CHROMOSOME1, s.CHROMOSOME2,
                                  s.SPANSTART, s.SPANEND, s.OVERLAPSTART, s.OVERLAPEND,   
                                  s.FRAGMENTMEAN, s.EVENTMEAN, s.STDEV, s.MEANDELTA, s.FRACTIONDELTA2, 
                                  s.STRAND1, s.STRAND2, s.NFRAGSTOTAL, s.NFRAGSSAMPLE, s.FRACBADFRAGS, s.MARK, s.DESCRIPTION ";
    my $newLimitsTable = getTableName('Sets', "$sample\_TMP");   
    dropTable($newLimitsTable);                            
    runSQL("CREATE TABLE $newLimitsTable AS $newLimitsSQL");
    #step 2, calculate fractionOverlaps from new limits and finalize all in revised Sets table
    my $leftRef =  "(CASE s.STRAND1 WHEN 2 THEN s.SPANSTART_NEW ELSE s.OVERLAPSTART_NEW END)";
    my $rightRef = "(CASE s.STRAND2 WHEN 1 THEN s.SPANEND_NEW   ELSE s.OVERLAPEND_NEW   END)";                                  
    my $modeNormal = "(f.FRAGMENTSIZE + f.EVENTSIZE)";                                  
    my $fragExcess = "(Abs($leftRef - f.POSITION1) + Abs($rightRef - f.POSITION2))";
    my $fragOverlap = "($modeNormal - $fragExcess)";  
    my $overlap = "Greatest(Avg($fragOverlap),0)";  #overlap can be negative! avoid the confusion this creates by disallowing    
    my $excessLeft =  "(Abs(s.OVERLAPSTART_NEW - s.SPANSTART_NEW) * 2)";
    my $excessRight = "(Abs(s.SPANEND_NEW - s.OVERLAPEND_NEW) * 2)"; 
    my $fractionOverlapLeft =  "(1 - ($excessLeft /  ($overlap + $excessLeft )))"; 
    my $fractionOverlapRight = "(1 - ($excessRight / ($overlap + $excessRight)))"; 
    $newLimitsSQL =    "SELECT s.SETID, s.SETTYPE, s.CHROMOSOME1, s.CHROMOSOME2,
                                s.SPANSTART_NEW SPANSTART, s.SPANEND_NEW SPANEND, s.OVERLAPSTART_NEW OVERLAPSTART, s.OVERLAPEND_NEW OVERLAPEND,    
                                s.FRAGMENTMEAN, s.EVENTMEAN, s.STDEV, s.MEANDELTA, s.FRACTIONDELTA2, 
                                $fractionOverlapLeft FRACTIONOVERLAPLEFT, $fractionOverlapRight FRACTIONOVERLAPRIGHT,  
                                s.STRAND1, s.STRAND2, s.NFRAGSTOTAL, s.NFRAGSSAMPLE, s.FRACBADFRAGS, s.MARK, s.DESCRIPTION 
                         FROM ($newLimitsTable) s, ($nonPromFragSQL) f
                         WHERE f.FRAGMENTTYPE = s.SETTYPE
                           AND f.CHROMOSOME1 = s.CHROMOSOME1 AND f.CHROMOSOME2 = s.CHROMOSOME2
                           AND f.POSITION1 >= s.SPANSTART AND f.POSITION1 <= s.OVERLAPSTART
                           AND f.POSITION2 >= s.OVERLAPEND AND f.POSITION2 <= s.SPANEND   
                         GROUP BY s.SETID, s.SETTYPE, s.CHROMOSOME1, s.CHROMOSOME2,
                                s.FRAGMENTMEAN, s.EVENTMEAN, s.STDEV, s.MEANDELTA, s.FRACTIONDELTA2,  
                                s.STRAND1, s.STRAND2, s.NFRAGSTOTAL, s.NFRAGSSAMPLE, s.FRACBADFRAGS, s.MARK, s.DESCRIPTION,
                                s.SPANSTART_NEW, s.SPANEND_NEW, s.OVERLAPSTART_NEW, s.OVERLAPEND_NEW ";                     
    updateTable('Sets', $setsTable, $newLimitsSQL); 
    dropTable($newLimitsTable);
}

sub findZSets{
    my ($sample) = @_; 
    my $setsTable = getTableName('Sets', $sample);
    runSQL("DELETE FROM $setsTable WHERE SETTYPE = $types{Sets}{Z}");
    runSQL("SELECT Max(SETID) FROM $setsTable", \my($setID));    
    fetchRow();
    my $fieldNames = $fieldNames{Sets};
    $fieldNames =~ s/, DESCRIPTION//;      
    my $ZMapTable = createZMap($sample);
    foreach my $chrom (1..nChrom()){    
        my $setsFile = "$setsTable.csv"; 
        open my $setsFileH, ">", $setsFile; 
        runSQL("SELECT z.POSITION, z.COVERAGE, z.MEAN, z.Z
                FROM $ZMapTable z, $setsTable s    
                WHERE z.CHROMOSOME = $chrom
                  AND z.CHROMOSOME = s.CHROMOSOME1               
                  AND z.POSITION >= s.SPANSTART
                  AND z.POSITION <= s.SPANEND
                  AND (s.SETTYPE = $types{Sets}{Deletion} OR s.SETTYPE = $types{Sets}{Insertion})
                  AND Abs(s.EVENTMEAN) < 3000
                  AND z.Z > 3
                GROUP BY z.POSITION, z.COVERAGE, z.MEAN, z.Z
                ORDER BY z.POSITION", 
                \my($pos, $coverage, $mean, $Z));    
        fetchRow();
        defined $pos or next;    
        my ($minPos, $prevPos, $maxCoverage, $maxZ, $meanSum, $nBins) = ($pos, $pos, $coverage, $Z, $mean, 1);    
        while (fetchRow()){    
            if ($pos - $prevPos > $param{binSize}){
                commitZSet($setsFileH, \$setID, $chrom, $minPos, $prevPos, $maxCoverage, $maxZ, int($meanSum / $nBins));
                ($minPos, $maxCoverage, $maxZ, $meanSum, $nBins) = ($pos, $coverage, $Z, 0, 0);
            }
            $meanSum += $mean;
            $nBins++;
            $maxCoverage >= $coverage or $maxCoverage = $coverage;
            $maxZ >= $Z or $maxZ = $Z;
            $prevPos = $pos;
        }
        commitZSet($setsFileH, \$setID, $chrom, $minPos, $prevPos, $maxCoverage, $maxZ, int($meanSum / $nBins));
        close $setsFileH;
        loadData($setsFile, $setsTable, ",", $fieldNames);         
    } 
    dropTable($ZMapTable);
}

sub createZMap{
    my ($sample) = @_; 
    my $fragsTable = getTableName('Frags', $sample);
    my $ZMapTable = newTable('ZMap', $sample);
    foreach my $chrom (1..nChrom()){    
        my (%bins); 
        #The purpose of the Z map is to find regions with Sets close to the main Normal peak
        #by considering Fragment types which could be miscalled as Normal, i.e. small Deletion and Insertion.
        #A Z score is calculated for each position bin considering small Deletion, Insertion AND Normal Fragments
        #to determine if the collection of Convergent Fragments crossing a bin 
        #deviates from the expected Normal event size distribution.
        #Deletion and Insertion Fragments are included ONLY if they were identified as being in a single Set.
        runSQL("SELECT ((trunc(POSITION1 / $param{binSize}) * $param{binSize}) + $param{binSize}) AS LOWBIN,
                       (trunc((POSITION2 + LENGTH2) / $param{binSize}) * $param{binSize}) AS HIGHBIN,
                       EVENTSIZE, STDEVNORMAL
                FROM $fragsTable    
                WHERE CHROMOSOME1 = $chrom 
                  AND ( FRAGMENTTYPE = $types{Frags}{Normal}                
                        OR ((FRAGMENTTYPE = $types{Frags}{Deletion} OR FRAGMENTTYPE = $types{Frags}{Insertion})
                            AND (NSETSFRAG = 1 AND NSETSPAIR = 1) ) )                
                  AND Abs(EVENTSIZE) <= 9 * STDEVNORMAL", 
                \my($lowBin, $highBin, $eventSize, $stdDevNormal));
        while (fetchRow()){
            for (my $bin = $lowBin; $bin<= $highBin; $bin += $param{binSize}){
                $bins{$bin}{count}++;
                $bins{$bin}{eventSizeSum} += $eventSize;                
                $bins{$bin}{popStdDevSum} += $stdDevNormal;            
            }
        }
        my $ZMapFile = "$ZMapTable.csv"; 
        open my $ZMapFileH, ">", $ZMapFile;        
        foreach my $bin(keys %bins){                     
            my ($Z, $binMean) = (0, 0);                      
            if ($bins{$bin}{count}){
                my $avePopStdDev = $bins{$bin}{popStdDevSum}/$bins{$bin}{count};
                $binMean = $bins{$bin}{eventSizeSum}/$bins{$bin}{count};
                $Z = abs($binMean) * sqrt($bins{$bin}{count}) / $avePopStdDev; #pop mean event size = 0  
            }   
            print $ZMapFileH join(",", ($chrom, $bin, $bins{$bin}{count}, int($binMean), $Z))."\n";            
        }                
        close $ZMapFileH;
        loadData($ZMapFile, $ZMapTable, ",", "CHROMOSOME, POSITION, COVERAGE, MEAN, Z");        
    }
    return $ZMapTable;
}

sub commitZSet{
    my ($setsFileH, $setIDRef, $chrom, $minPos, $maxPos, $maxCoverage, $maxZ, $mean) = @_;
    $$setIDRef++;
    print $setsFileH join(",", ($$setIDRef, $types{Sets}{Z},
                                $chrom, 0,
                                $minPos, $maxPos, $minPos, $maxPos,
                                $maxZ, $mean, $maxZ, 
                                0, 0, 0, 0,
                                1, 2,
                                $maxCoverage, $maxCoverage, 0,
                                0))."\n";
}

##this works, but in general we do not have enough data density (or consistency in the first samples?) 
##for it to work effectively.  There are too many coverge differences between samples...
#sub analayzeFMapArrays{
#    my ($refSample, $testSample) = @_;;
#    my ($refFMapTable, $testFMapTable) = (getTableName('FMap', $refSample), getTableName('FMap', $testSample));    
#    my ($refArrayTable, $testArrayTable) = (newTable('Array', arrayTableName($refSample,  'MatePair', $param{refSeq})), 
#                                            newTable('Array', arrayTableName($testSample, 'MatePair', $param{refSeq})));    
#    status("translating $refSample Fragment coverage map into Array table format\n"); 
#    runSQL("INSERT INTO $refArrayTable 
#            SELECT CHROMOSOME, POSITION,
#                    NORMALIZEDCOVERAGE RATIO, NORMALIZEDCOVERAGE NORMALIZEDRATIO,
#                    -1 BFREQUENCY, -1 ZYGOSITY, -1 NORMALIZEDZYGOSITY, -1 CNVVALUE,
#                    1 GCSCORE, 0 X, 0 Y, 0 THETA,
#                    NORMALIZEDCOVERAGE R,
#                    -1 BF
#            FROM $refFMapTable");  
#    status("translating $testSample Fragment coverage map into Array table format"); 
#    runSQL("INSERT INTO $testArrayTable 
#            SELECT r.CHROMOSOME, r.POSITION,
#                    nvl(t.NORMALIZEDCOVERAGE,0) RATIO, nvl(t.NORMALIZEDCOVERAGE,0) NORMALIZEDRATIO,
#                    -1 BFREQUENCY, -1 ZYGOSITY, -1 NORMALIZEDZYGOSITY, -1 CNVVALUE,
#                    1 GCSCORE, 0 X, 0 Y, 0 THETA,
#                    nvl(t.NORMALIZEDCOVERAGE,0) R,
#                    -1 BF
#            FROM $refFMapTable r, $testFMapTable t
#            WHERE r.CHROMOSOME = t.CHROMOSOME(+)
#              AND r.POSITION = t.POSITION(+)");   
#              
#    $param{nStdDevs} = '4';
#    $param{nProbes} = '20,50';  
#     
#    $param{arrayType} = 'MatePair';   
#    $param{arrayRef} = $refSample;
#    analyzeArray($testSample);
#}




1;



