#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));

require "$param{vampPath}/bin/findIntersection.pl";

my ($majorDiv, $minorDiv) = ("=" x 50, "-" x 50);
my ($CNVsTable, $idField, $nIterations, $padding, $maxFragSize, 
    $cnvSimTable, $genomeSize, $nChrom, $maxI, $chromsRef, $cnvsRef, $nCNVs, %results);

sub findCNVHotspots{
    ($CNVsTable, $idField, $nIterations, $padding, $maxFragSize) = @_;
    initializeCNVSimulation();       
    status("analyzing the actual data...\n$majorDiv\n");  
    #my $type = "CNV Hotspots::Max Coverage";
    #mapAllCNVS($CNVsTable, \%{$results{$type}{Actual}}); 
    my $type = "CNV Hotspots::Number of CNVs in region";
    countOverlappingCNVs($CNVsTable, \%{$results{$type}{Actual}}, 1); 
    status("$majorDiv\nanalyzing $nIterations simulations...\n");
    for my $iteration(1..$nIterations){
        status("$iteration ");
        generateCNVSimulation(); 
        #mapAllCNVS($cnvSimTable, \%{$results{$type}{Simulated}}); 
        countOverlappingCNVs($cnvSimTable, \%{$results{$type}{Simulated}}, 0); 
    } 
    status("\nreporting results...\n");
    printSimResults();
}

sub simulateAllCNVs{
    ($CNVsTable, $idField, $nIterations) = @_;
    initializeCNVSimulation();     
    status("analyzing the actual data...\n");
    runSimIteration($CNVsTable, 'Actual');     
    status("analyzing $nIterations simulations...\n");
    for my $iteration(1..$nIterations){
        status("$iteration ");
        generateCNVSimulation();  
        runSimIteration($cnvSimTable, 'Simulated');    
    }    
    status("\nreporting results...\n");
    printSimResults();
}

sub runSimIteration{
    my ($CNVsTable, $actual) = @_;
    runSimIterationElement($actual, "$CNVsTable\_Hu_Ov", 'Overlap', 10,
                           "Number of events crossing human CNVs::Score Threshold", 
                           $CNVsTable, "CHROMOSOME", "START_", "END_", "CNVS_HUMAN_HG18", "CHROMOSOME", "START_", "END_");
    runSimIterationElement($actual, "$CNVsTable\_Hu_ED", 'EndDistance', 1E3,
                           "Number of event ends MinDistance from a human CNV end::MinDistance", 
                           $CNVsTable, "CHROMOSOME", "START_", "END_", "CNVS_HUMAN_HG18", "CHROMOSOME", "START_", "END_",
                           1E4); 
}

sub runSimIterationElement{
    my ($actual, $outputTable, $elementType, $scalar, $header, @dataParam) = @_;
    my ($genRef, $findRef);
    if($elementType eq 'Overlap'){
        ($genRef, $findRef) = (\&generateTargetOverlapTable, \&findTargetOverlaps);
    } elsif ($elementType eq 'EndDistance'){
        ($genRef, $findRef) = (\&generateTargetEndDistanceTable, \&findTargetEndDistances);
    }
    $actual eq 'Actual' and &$genRef(@dataParam, 0, $outputTable);
    &$findRef(@dataParam, \%{$results{$header}{$actual}}, $scalar);   
}

sub printSimResults{
    foreach my $type(keys %results){
        my ($majorHeader, $minorHeader) = split("::", $type);
        status("\n$majorDiv\n$majorHeader\n$minorDiv\n"); 
        status("\tActual\t\tSimulated\n"); 
        status("$minorHeader\tNumber of regions\tTotal CNVs\tNumber of regions\tTotal CNVs\n");  
        ##while(fetchRow()){$$resultsRef{$maxCov}{$count}++};   
        ##$$resultsRef{$cnvCount}{$hsCount}++; 
        foreach my $cnvCount(sort {$a <=> $b} keys %{$results{$type}{Actual}}){ #$majorVal = $maxCov or $cnvCount
            my @nHsAct = keys %{$results{$type}{Actual}{$cnvCount}}; #$minorVal = $count in 1 actual iteration
            my $nHsAct = $nHsAct[0]; #iterCount == 1 for actual
            my $totalCnvsAct = $nHsAct * $cnvCount;
            status("$cnvCount\t$nHsAct\t$totalCnvsAct\t");    
            if ($results{$type}{Simulated}{$cnvCount}){ #$majorVal = $maxCov    
                my @nHsSim = keys %{$results{$type}{Simulated}{$cnvCount}}; #$minorVal = $count in 1 actual iteration
                my $nHsSim;
                foreach my $hsCount(@nHsSim){ 
                    my $iterCount = $results{$type}{Simulated}{$cnvCount}{$hsCount}; #number of iterations in which count was obtained
                    $nHsSim += ($hsCount * $iterCount);
                }
                my $totalCnvsSim = $nHsSim * $cnvCount;
                $nHsSim /= $nIterations;
                $totalCnvsSim /= $nIterations;
                status("$nHsSim\t$totalCnvsSim\n"); 
                # my @simulatedValues = keys %{$results{$type}{Simulated}{$majorVal}}; #$minorVal = $count over all simulations
                # my ($sum, $sumDev2, $iterCount);    
                # foreach my $minorVal(@simulatedValues){ 
                    # my $count = $results{$type}{Simulated}{$majorVal}{$minorVal}; #number of iterations in which count was obtained
                    # $sum += ($minorVal * $count); 
                # }
                # my $mean = $sum/$nIterations;
                # foreach my $minorVal(@simulatedValues){ 
                    # my $count = $results{$type}{Simulated}{$majorVal}{$minorVal};
                    # $iterCount += $count;
                    # $sumDev2 += ((($minorVal - $mean) ** 2) * $count); 
                # }
                # my $missingIterations = $nIterations - $iterCount;
                # $sumDev2 += (((0 - $mean) ** 2) * $missingIterations) ;
                # my $stdDev = sqrt($sumDev2 / ($nIterations - 1));
                #status("$mean\t$stdDev\n");
            } else {
                status("0\t0\n"); 
            }
        }
        status("$majorDiv\n"); 
    }
}
sub initializeCNVSimulation{
    $nChrom = nChrom();
    $param{sex} and $param{sex} eq 'XX' and $nChrom--;
    my $chromInfoTable = "CHROMINFO_$param{refSeqBase}";
    runSQL("SELECT Sum(End_) FROM $chromInfoTable WHERE CHROMOSOME <= $nChrom", \($genomeSize));
    fetchRow();
    runSQL("SELECT CHROMOSOME, END_ FROM $chromInfoTable WHERE CHROMOSOME <= $nChrom ORDER BY CHROMOSOME");
    $chromsRef = fetchAllHashRef();
    $maxI = scalar(@$chromsRef) - 1;
    runSQL("SELECT $idField, CHROMOSOME, START_, END_, (END_ - START_) CNVSIZE 
            FROM $CNVsTable
            WHERE CHROMOSOME <= $nChrom
              AND END_ - START_ < $maxFragSize");    
    $cnvsRef = fetchAllHashRef();   
    $nCNVs = scalar @$cnvsRef;
    $cnvSimTable = "$CNVsTable\_SIM";  
}

sub generateCNVSimulation {
    dropTable($cnvSimTable);
    runSQL("CREATE TABLE $cnvSimTable ($idField VARCHAR2(255), CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER)"); 
    my $cnvSimFile = "$cnvSimTable.csv";  
    open my $cnvSimFileH, ">", $cnvSimFile;   
    CNVSIM: foreach my $cnvRef(@$cnvsRef){
        my $id = $$cnvRef{$idField}; 
        my $cnvSize = $$cnvRef{CNVSIZE};
        TRY_AGAIN:
        my $cnvIndex = int(rand $genomeSize);
        my ($cumGenomeSize, $cnvChrom, $cnvStart, $cnvEnd) = (0);
        foreach my $i(0..$maxI){
            if($cnvIndex <= ($cumGenomeSize + $$chromsRef[$i]{END_})){ 
                $cnvChrom = $$chromsRef[$i]{CHROMOSOME};          
                $cnvStart = $cnvIndex - $cumGenomeSize;        
                $cnvEnd =  $cnvStart + $cnvSize;
                if(inGap($cnvChrom, $cnvStart, $cnvEnd)){goto TRY_AGAIN}  #simulated CNVs cannot reside in genome gaps
                if(inGap($cnvChrom, $cnvStart, $cnvEnd, 'omniQuad')){goto TRY_AGAIN}  #simulated CNVs cannot reside in array gaps
                if($cnvEnd <= $$chromsRef[$i]{END_}){              
                    commitCNVSim($cnvSimFileH, $id, $cnvChrom, $cnvStart, $cnvEnd); 
                } else { #simulated CNVs cannot wrap chromosomes
                    goto TRY_AGAIN;
                    # commitCNVSim($cnvSimFileH, $id, $cnvChrom, $cnvStart, $$chromsRef[$i]{END_});
                    # my $j = $i + 1;
                    # $j > $maxI and $j = 0;
                    # $cnvChrom = $$chromsRef[$j]{CHROMOSOME}; 
                    # $cnvEnd = $cnvSize - ($$chromsRef[$i]{END_} - $cnvStart);                    
                    # commitCNVSim($cnvSimFileH, $id, $cnvChrom, 1, $cnvEnd);                    
                }
                next CNVSIM;
            }
            $cumGenomeSize += $$chromsRef[$i]{END_};
        }
    }
    close $cnvSimFileH;
    loadData($cnvSimFile, $cnvSimTable, ",", "$idField, CHROMOSOME, START_, END_"); 
}

#sub inGap{ #demands that at least 67% of the CNV exist outside of a gap
#    my ($chrom, $start, $end) = @_;
#    my $size = $end - $start;
#    my $typeSQL = "SELECT START_, END_, SIZE_,
#                          CASE WHEN (START_ <= $start AND END_ >= $end) THEN 1
#                               WHEN (START_ > $start AND END_ < $end) THEN 2
#                               WHEN (START_ <= $end AND END_ <= $end) THEN 3.1
#                               ELSE 3.2 END OVERLAPTYPE
#                   FROM GAP_$param{refSeqBase}
#                   WHERE CHROMOSOME = $chrom
#                     AND START_ <= $end
#                     AND END_ >= $start ";    
#    my $overlapSQL =  "SELECT CASE OVERLAPTYPE WHEN 1 THEN $size
#                                               WHEN 2 THEN SIZE_
#                                               WHEN 3.1 THEN END_ - $start
#                                               ELSE $end - START_ END OVERLAP    
#                        FROM ($typeSQL) ";    
#    runSQL("SELECT Count(*) FROM ($overlapSQL) WHERE OVERLAP/$size > 0.3", \my($count));
#    fetchRow();
#    return $count;
#}

sub inGap{ #demands that at least 50% of the CNV exist outside of a gap
    my ($chrom, $start, $end, $refSeq) = @_;
    $refSeq or $refSeq = $param{refSeqBase};
    my $size = $end - $start;
    my $typeSQL = "SELECT START_, END_, SIZE_,
                          CASE WHEN (START_ <= $start AND END_ >= $end) THEN 1
                               WHEN (START_ > $start AND END_ < $end) THEN 2
                               WHEN (START_ <= $end AND END_ <= $end) THEN 3.1
                               ELSE 3.2 END OVERLAPTYPE
                   FROM GAP_$refSeq
                   WHERE CHROMOSOME = $chrom
                     AND START_ <= $end
                     AND END_ >= $start ";    
    my $overlapSQL =  "SELECT CASE OVERLAPTYPE WHEN 1 THEN $size
                                               WHEN 2 THEN SIZE_
                                               WHEN 3.1 THEN END_ - $start
                                               ELSE $end - START_ END OVERLAP    
                        FROM ($typeSQL) ";    
    runSQL("SELECT sum(OVERLAP) FROM ($overlapSQL)", \my($inGap));
    fetchRow();
    $inGap or return undef;
    return ($inGap/$size > 0.5);
}

sub commitCNVSim{
    my ($cnvSimFileH, $id, $cnvChrom, $cnvStart, $cnvEnd) = @_;
    print $cnvSimFileH "$id,$cnvChrom,$cnvStart,$cnvEnd\n";
}

sub countOverlappingCNVs {
    my ($CNVsTable, $resultsRef, $isActual) = @_;
    my %cnvCounts;
    $isActual and status("CHROMOSOME\tSTART\tEND\tCNV_COUNT\n");
    foreach my $chrom(1..$nChrom){
        runSQL("SELECT START_, END_ 
                   FROM $CNVsTable
                   WHERE CHROMOSOME = $chrom
                     AND END_ - START_ < $maxFragSize
                   ORDER BY START_, END_");
        my $cnvsRef = fetchAllHashRef();
        my $nChromCnvs = scalar(@$cnvsRef);
        $nChromCnvs > 0 or next;  #no CNVs on chromosome 
        my ($hsStart, $hsEnd, $cnvCount) = ($$cnvsRef[0]{START_}, $$cnvsRef[0]{END_}, 1); #initalize first hotspot
        if($nChromCnvs > 1){
            foreach my $i(1..($nChromCnvs-1)){
                if($$cnvsRef[$i]{START_} > $hsEnd + $padding){ #not in same overlap group
                    $cnvCounts{$cnvCount}++; 
                    $isActual and status("$chrom\t$hsStart\t$hsEnd\t$cnvCount\n");
                    ($hsStart, $hsEnd, $cnvCount) = ($$cnvsRef[$i]{START_}, $$cnvsRef[$i]{END_}, 1); #re-initialize
                } else { #still in same overlap group
                    $hsEnd >= $$cnvsRef[$i]{END_} or $hsEnd = $$cnvsRef[$i]{END_};
                    $cnvCount++; #the number of CNVs within this hotspot
                }     
            }  
            $cnvCounts{$cnvCount}++; 
            $isActual and status("$chrom\t$hsStart\t$hsEnd\t$cnvCount\n");  
        } else {
            $cnvCounts{$cnvCount}++; 
            $isActual and status("$chrom\t$hsStart\t$hsEnd\t$cnvCount\n");
        }
    }
    foreach my $cnvCount(keys %cnvCounts){
        my $hsCount = $cnvCounts{$cnvCount}; #the number of hotspots containing that many overlapping CNVs
        $$resultsRef{$cnvCount}{$hsCount}++; #this count always 1 for actual sample
    }
}

sub mapAllCNVS{
    my ($CNVsTable, $resultsRef) = @_;
    my $binSize = 1E4;
    
    #my $windowSize = 1E6;
    my $windowSize = 0;
    my $halfWindowSize = $windowSize / 2;
    
    my $mapTable = "$CNVsTable\_MAP";
    dropTable($mapTable);
    runSQL("CREATE TABLE $mapTable (CHROMOSOME NUMBER, POSITION NUMBER, COVERAGE NUMBER)");    
    foreach my $chrom(1..$nChrom){
        my %map;   
        runSQL("SELECT ((trunc(START_/$binSize) * $binSize) + $binSize) AS LOWBIN, (trunc(END_/$binSize) * $binSize) AS HIGHBIN
                FROM $CNVsTable
                WHERE CHROMOSOME = $chrom",
                \my($lowBin, $highBin));
        while (fetchRow()){
            for (my $bin = $lowBin - $halfWindowSize; $bin <= $highBin + $halfWindowSize; $bin += $binSize){
                $map{$bin}++;
            }
        }
        my $mapFile = "$mapTable.csv";  #load the results into map table
        open my $mapFileH, ">", $mapFile;
        foreach my $bin (keys %map){
            print $mapFileH join(",", ($chrom, $bin, $map{$bin}))."\n";
        }
        close $mapFileH;
        loadData($mapFile, $mapTable, ",", "CHROMOSOME, POSITION, COVERAGE");         
    }
    my $hsTable = "$CNVsTable\_HS";
    dropTable($hsTable);
    runSQL("CREATE TABLE $hsTable (CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER, MAXCOVERAGE NUMBER)"); 
    foreach my $chrom(1..$nChrom){
        runSQL("SELECT POSITION, COVERAGE
                FROM $mapTable
                WHERE CHROMOSOME = $chrom AND COVERAGE > 1
                ORDER BY POSITION", \my($pos, $coverage));  
        fetchRow();
        $pos or next;
        my ($start, $maxCoverage, $prevPos) = ($pos, $coverage, $pos);
        my $hsFile = "$hsTable.csv";  
        open my $hsFileH, ">", $hsFile;
        while (fetchRow()){
            if(($pos - $prevPos) > $binSize){ 
                print $hsFileH join(",", ($chrom, $start, $prevPos, $maxCoverage))."\n";
                $start = $pos;
                $maxCoverage = $coverage;  
            }
            $maxCoverage >= $coverage or $maxCoverage = $coverage;
            $prevPos = $pos;
        } 
        print $hsFileH join(",", ($chrom, $start, $prevPos, $maxCoverage))."\n";
        close $hsFileH;
        loadData($hsFile, $hsTable, ",", "CHROMOSOME, START_, END_, MAXCOVERAGE"); 
    }
    runSQL("SELECT MAXCOVERAGE, Count(MAXCOVERAGE) FROM $hsTable GROUP BY MAXCOVERAGE", \my($maxCov, $count));    
    while(fetchRow()){$$resultsRef{$maxCov}{$count}++};    
}

sub findTargetOverlaps{
    my ($queryTable, $qChromField, $qStartField, $qEndField,
        $targetTable, $tChromField, $tStartField, $tEndField,
        $resultsRef, $scalar) = @_;
    generateTargetOverlapTable($queryTable, $qChromField, $qStartField, $qEndField,
                               $targetTable, $tChromField, $tStartField, $tEndField,
                               "bestScore", "SIMULATION");
    stratifySimulationResults("SIMULATION", "SCORE_T", 1, $resultsRef, $scalar);                         
}  

sub generateTargetOverlapTable{
    my ($queryTable, $qChromField, $qStartField, $qEndField,
        $targetTable, $tChromField, $tStartField, $tEndField,
        $outputType, $outputTable) = @_;
    my @intInstrs = (["findIntersection"],
                     ["fromCode"],
                     ["query", $queryTable],
                     ["chrom", $qChromField],
                     ["start", $qStartField],
                     ["end", $qEndField],
                     ["include", $idField], ##############
                     ["target", $targetTable, "T"],
                     ["chrom", $tChromField],
                     ["start", $tStartField],
                     ["end", $tEndField],
                     ["outputType", $outputType],
                     ["output", $outputTable] );
    findIntersection(\@intInstrs);                      
}
      
sub findTargetEndDistances{
    my ($queryTable, $qChromField, $qStartField, $qEndField,
        $targetTable, $tChromField, $tStartField, $tEndField, 
        $maxDistance, $resultsRef, $scalar) = @_;
    generateTargetEndDistanceTable($queryTable, $qChromField, $qStartField, $qEndField,
                                   $targetTable, $tChromField, $tStartField, $tEndField, 
                                   $maxDistance, "bestScore", "SIMULATION");
    stratifySimulationResults("SIMULATION", "DISTANCE", 0, $resultsRef, $scalar); 
}

sub generateTargetEndDistanceTable{
    my ($queryTable, $qChromField, $qStartField, $qEndField,
        $targetTable, $tChromField, $tStartField, $tEndField, 
        $maxDistance, $outputType, $outputTable) = @_;
    my $halfMaxDist = int(($maxDistance / 2) + 0.5);         
    my $queryView = createEndSpanView($queryTable, $qChromField, $qStartField, $qEndField, $halfMaxDist, $idField); #############
    my $targetView = createEndSpanView($targetTable, $tChromField, $tStartField, $tEndField, $halfMaxDist);          
    my @intInstrs = (["findIntersection"],
                     ["fromCode"],
                     ["query", $queryView],
                     ["include", "$idField, POSFIELD"],                     
                     ["target", $targetView, "T"],
                     ["include", "POSFIELD"],  
                     ["outputType", $outputType],
                     ["output", $outputTable] );
    findIntersection(\@intInstrs);  
    runSQL("ALTER TABLE $outputTable ADD DISTANCE NUMBER");
    runSQL("UPDATE $outputTable SET DISTANCE = $maxDistance * (1 - (SCORE_T / 100))");        
}

sub createEndSpanView{
    my ($table, $chromField, $startField, $endField, $halfMaxDist, $include) = @_;
    ($include and $include = ", $include") or $include = '';
    my $view = "$table\_Ends";
    dropView($view);
    runSQL("CREATE VIEW $view AS
            SELECT $chromField CHROMOSOME, $startField - $halfMaxDist START_, $startField + $halfMaxDist END_, '$startField' POSFIELD $include
            FROM $table
            UNION ALL
            SELECT $chromField CHROMOSOME, $endField - $halfMaxDist START_, $endField + $halfMaxDist END_, '$endField' POSFIELD $include
            FROM $table ");
    return $view;
}

sub stratifySimulationResults{
    my ($inputSQL, $indexField, $sortDesc, $resultsRef, $scalar) = @_;
    my $sort = '';
    $sortDesc and $sort = 'DESC';
    my $scaleSQL = "SELECT Round($indexField / $scalar) * $scalar $indexField FROM ($inputSQL)";    
    my $countSQL = "SELECT $indexField, Count(*) N FROM ($scaleSQL) GROUP BY $indexField";     
    runSQL("SELECT $indexField, Sum(N) OVER (ORDER BY $indexField $sort ROWS UNBOUNDED PRECEDING)
            FROM ($countSQL)
            ORDER BY $indexField $sort", 
            \my($majorVal, $minorVal));                       
    while(fetchRow()){$$resultsRef{$majorVal}{$minorVal}++};    
}

#sub findCrossingHumanCNVs{
#    my ($CNVsTable, $resultsRef) = @_;
#    my @intInstrs = (["findIntersection"],
#                     ["fromCode"],
#                     ["query", $CNVsTable],
#                     ["include", $idField],
#                     ["target", "HUMAN_CNVS", "H"],
#                     ["outputType", "bestScore"],
#                     ["output", "SIMULATION"] );
#    findIntersection(\@intInstrs);   
#    my $binSQL = "SELECT Round(SCORE_H/10)*10 SCORE_H FROM SIMULATION"; 
#    my $countSQL = "SELECT SCORE_H, Count(SCORE_H) N FROM ($binSQL) GROUP BY SCORE_H";
#    runSQL("SELECT SCORE_H, Sum(N) OVER (ORDER BY SCORE_H DESC ROWS UNBOUNDED PRECEDING)
#            FROM ($countSQL)
#            ORDER BY SCORE_H DESC", 
#            \my($majorVal, $minorVal));    
#    while(fetchRow()){$$resultsRef{$majorVal}{$minorVal}++};
#}





#sub findEndingHumanCNVs{
#    my ($CNVsTable, $span, $resultsRef) = @_;
#    my $halfSpan = int(($span / 2) + 0.5);    
#    my $CNVsView = createEndSpanView($CNVsTable, 'CHROMOSOME', 'START_', 'END_', $halfSpan, $idField);
#    my $humanView = createEndSpanView("HUMAN_CNVS", 'CHROMOSOME', 'START_', 'END_', $halfSpan);       
#    my @intInstrs = (["findIntersection"],
#                     ["fromCode"],    
#                     ["query", $CNVsView],    
#                     ["include", "$idField, LIMITFIELD"],
#                     ["target", $humanView, "T"],       
#                     ["outputType", "bestScore"],
#                     ["output", "SIMULATION"] );
#    findIntersection(\@intInstrs);    
#    my $distSQL = "SELECT Round($span * (1 - (SCORE_T / 100))/1E3) MINDISTANCE 
#                   FROM SIMULATION"; 
#    my $countSQL = "SELECT MINDISTANCE, Count(MINDISTANCE) N FROM ($distSQL) GROUP BY MINDISTANCE";
#    runSQL("SELECT MINDISTANCE, Sum(N) OVER (ORDER BY MINDISTANCE ROWS UNBOUNDED PRECEDING)
#            FROM ($countSQL)
#            ORDER BY MINDISTANCE", 
#            \my($majorVal, $minorVal));    
#    while(fetchRow()){$$resultsRef{$majorVal}{$minorVal}++};
#}







1;












