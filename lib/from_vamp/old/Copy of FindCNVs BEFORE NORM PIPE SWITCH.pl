#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#microarray analysis
#########################################################################

use vars(qw(%param %types %fields %fieldNames %refSeqs));  #common parameters used by most analysis scripts

my ($testSample, $nStdDev, @nStdDevs, $nProbes, $pf, $nmaFileH, $nmaID, $nmaType);
my $scalar = 0.001;
my ($refSelectSQL, $testSelectSQL, $joinSQL, $maSQL, $nmaSQL, $nmaTable, $nmaDataTable);
my ($nmaMean, $nmaStDev, $refL2Mean, $refL2StDev, $testL2Mean, $testL2StDev, %tTable);  
my $posID = " (CHROMOSOME || 'x' || POSITION) "; 

sub importArray{
    my ($run) = @_;
    if($param{arrayType} eq 'Illumina'){
        importIllumina($run)
    }
}

sub arrayTableName{
    my ($sample, $arrayType, $refSeq) = @_;
    my %shortTypes = (Nimblegen => 'NIM', Illumina => 'ILL');
    return "$sample\_$shortTypes{$arrayType}\_$refSeq";
}

sub analyzeArray{
    ($testSample) = @_;
    initializeTTable();
    my @nProbes = split(",", $param{nProbes});
    @nStdDevs = split(",", $param{nStdDevs});
    $nmaTable = newTable('NMA', "$param{arrayRef}\_$testSample");       
    my $nmaFile = "$nmaTable.csv"; 
    open $nmaFileH, ">", $nmaFile; 
    foreach my $nmaType_('RATIO', 'BFREQUENCY'){
        $nmaType = $nmaType_;
        status("nmaType = $nmaType\n"); 
        ($refL2Mean, $refL2StDev, $testL2Mean, $testL2StDev) = (undef, undef, undef, undef);
        foreach my $nProbes_(sort {$a <=> $b} @nProbes){
            $nProbes = $nProbes_;
            $nProbes % 2 or $nProbes++; 
            $pf = ($nProbes - 1) / 2;
            status("  nProbes = $nProbes\n");  
            if ($param{arrayType} eq 'Illumina'){  
                createFindCNVSetsSQL();
                getArrayStats(); 
                findCNVSets();             
            } else {
                die "array type $param{arrayType} not supported";
            }  
        }
    }
    close $nmaFileH;
    loadData($nmaFile, $nmaTable, ",", $fieldNames{NMA});
    open $nmaFileH, ">", $nmaFile;
    findCNVBlocks($nmaTable);
    close $nmaFileH;
    loadData($nmaFile, $nmaTable, ",", $fieldNames{NMA});
}

sub createFindCNVSetsSQL{
    $refSelectSQL =  getArraySelectSQL($param{arrayRef}, 1);
    $testSelectSQL = getArraySelectSQL($testSample, 0);
    
    
    ################################3
    #no, much better to track the log2(TVAL/RVAL) for piping that simply log2(TVAL)
    #bottom line is that probes tend to correlate somewhat
    #so that the raw stdev of normalized data is narrower than unnormalized
    #this also makes the piping code simpler
    
    
    $joinSQL = " SELECT r.CHROMOSOME, r.POSITION, r.VAL RVAL, t.VAL TVAL, r.L2 RL2, t.L2 TL2
                 FROM ($refSelectSQL) r, ($testSelectSQL) t
                 WHERE r.POSID = t.POSID 
                 ORDER BY CHROMOSOME, POSITION ";     
    $maSQL =   " SELECT CHROMOSOME, POSITION, RVAL, TVAL, RL2, TL2, 
                        Avg(RVAL) OVER (ORDER BY CHROMOSOME, POSITION ROWS BETWEEN $pf PRECEDING AND $pf FOLLOWING) RMA,
                        Avg(TVAL) OVER (ORDER BY CHROMOSOME, POSITION ROWS BETWEEN $pf PRECEDING AND $pf FOLLOWING) TMA
                 FROM ($joinSQL) ";    
    $nmaSQL = " SELECT CHROMOSOME, POSITION, RVAL, TVAL, RL2, TL2,
                       Round(Log(2,(TMA/nullif(RMA,0)))/$scalar)*$scalar NMA
                FROM ($maSQL) ";
    $nmaDataTable = "\U$nmaTable\_d";
#    dropTable($nmaDataTable);
#    status("    calculating normalized moving averages...\n");  
#    runSQL("CREATE TABLE $nmaDataTable AS $nmaSQL");           
    
}

sub getArraySelectSQL{
    my ($sample, $isRef) = @_;
    my $arrayTable = getTableName('Array', arrayTableName($sample, $param{arrayType}, $param{refSeq}));
    my $refFilter = " WHERE RATIO > 0.333 AND RATIO < 3.0 ";  
    if ($nmaType eq 'BFREQUENCY'){  
        my $zygosity = " (Abs(BFREQUENCY - 0.5) + 0.5) ";
        my $filter = '';
        $isRef and $filter = " $refFilter AND BFREQUENCY > 0.25 AND BFREQUENCY < 0.75 ";
        return " SELECT $posID POSID, CHROMOSOME, POSITION, $zygosity VAL, Round(Log(2,$zygosity/0.5))/$scalar)*$scalar L2 
                 FROM $arrayTable $filter "; 
    } elsif ($nmaType eq 'RATIO'){
        my $filter = '';
        $isRef and $filter = $refFilter;
        return " SELECT $posID POSID, CHROMOSOME, POSITION, RATIO VAL, Round(Log(2,(CASE RATIO WHEN 0 THEN 0.001 ELSE RATIO END))/$scalar)*$scalar L2 
                 FROM $arrayTable $filter"; 
    }
}

sub getArrayStats{
    #my ($nmaMean, $nmaStDev, $refL2Mean, $refL2StDev, $testL2Mean, $testL2StDev);
    status("    getting array stats...\n"); 
    ($nmaMean, $nmaStDev) = getStatsFromL2('NMA');
    status("      normalized moving average mean = $nmaMean, stDev = $nmaStDev\n"); 
    unless(defined $testL2Mean){
#        ($refL2Mean, $refL2StDev) = getStatsFromL2('RL2', $joinSQL);
#        status("      $param{arrayRef} mean = $refL2Mean, stDev = $refL2StDev\n"); 
        ($testL2Mean, $testL2StDev) = getStatsFromL2('TL2');
        status("      $testSample mean = $testL2Mean, stDev = $testL2StDev\n");  
    } 
}

sub getStatsFromL2{
    my ($field) = @_;
    my $histSQL = " SELECT $field X, Count($field) Y FROM $nmaDataTable WHERE $field != 0 GROUP BY $field ";
    runSQL("SELECT Max(Y) FROM ($histSQL)", \my($maxY));
    fetchRow();  
    my $minY = 0.05 * $maxY;
    $histSQL = " SELECT X, Y FROM ($histSQL) WHERE Y >= $minY ";   
    runSQL("SELECT Min(X), Max(X) FROM ($histSQL)", \my($minX, $maxX));
    fetchRow();  
    my $sGuess = ($maxX - $minX) / (2 * 2.45); #expectations of stdDev for perfect Gaussian given minY at 5% of maxY, mean = 0
    my ($mean, $stDev) = fitGaussian($histSQL, $maxY, 0, $sGuess);   
    $mean =  int(($mean  / $scalar) + 0.5) * $scalar;
    $stDev = int(($stDev / $scalar) + 0.5) * $scalar;     
    return ($mean, $stDev);
}

sub findCNVSets{
    status("    finding CNV sets...\n");
    foreach my $nStdDev_(@nStdDevs){
        $nStdDev = $nStdDev_;
        status("      nStdDev = $nStdDev: Chr "); 
        my $threshold = $nmaStDev * $nStdDev;      
        foreach my $chrom (1..$refSeqs{$param{refSeqBase}}{nChrom}){
            status("$chrom ");
            runSQL("SELECT POSITION, RVAL, TVAL, RL2, TL2, NMA
                    FROM $nmaDataTable
                    WHERE CHROMOSOME = $chrom
                    ORDER BY POSITION");    
            my $probesRef = fetchAllHashRef();       
            my ($minI, $maxI, $indexNMA, $lowI, $highI) = (0, (scalar @$probesRef) - 1, 0);                    
            foreach my $i($minI..$maxI){
                my $nma = ${$$probesRef[$i]}{NMA};
                my $overThreshold = (abs($nma - $nmaMean) >= $threshold);
                if (!$overThreshold){
                    if ($indexNMA){
                        checkCNVSet($chrom, $probesRef, $minI, $maxI, $lowI, $i - 1, $threshold);
                        $indexNMA = 0;
                    }                
                } else {
                    if ($indexNMA){
                        if (($nma/$indexNMA) < 0){
                            checkCNVSet($chrom, $probesRef, $minI, $maxI, $lowI, $i - 1, $threshold);
                            $lowI = $i; 
                            $indexNMA = $nma;
                        }
                    } else {
                        $lowI = $i; #start new set
                        $indexNMA = $nma;
                    }  
                }                 
            }
            $indexNMA and checkCNVSet($chrom, $probesRef, $minI, $maxI, $lowI, $maxI, $threshold);  
        }  
        status("\n");        
    }
    #dropTable($nmaDataTable);
}

sub checkCNVSet{
    my ($chrom, $probesRef, $minI, $maxI, $lowI, $highI, $threshold) = @_;
    ($lowI, $highI) = trimCNVSetEnds($probesRef, $minI, $maxI, $lowI, $highI, $threshold);  
    $lowI or return;
    #could take test mean and refmean from trimsetends, but doesn't return stddev, so pValue is tough
    #problem with stdev is that the span will likely be over reached
    #could take pValue from population stdev
    
    
      
    my ($n, $testMean, $refMean, $testSD, $refSD, $pValue) = getCNVSetStats($probesRef, $lowI, $highI);
    $n >= $nProbes or return; #is this necessary? desirable?
    
    #status(".");
    
    $refMean == 0 and $refMean = 0.001;  
    $testMean == 0 and $testMean = 0.001;  
    my $normMean = logBase2($testMean / $refMean);    
    #abs($normMean - $nmaMean) >= $threshold or return;     
    $nmaID++;    
    my ($normRatio, $normZyg);
    my $refTable =  getTableName('Array', arrayTableName($param{arrayRef}, $param{arrayType}, $param{refSeq}));
    my $testTable = getTableName('Array', arrayTableName($testSample,      $param{arrayType}, $param{refSeq}));
    my ($minPos, $maxPos) = (${$$probesRef[$lowI]}{POSITION}, ${$$probesRef[$highI]}{POSITION});
    if ($nmaType eq 'RATIO'){
        ($normRatio, $normZyg) = ($normMean, getCrossZygosity($refTable, $testTable, $chrom, $minPos, $maxPos));
    } else {
        ($normRatio, $normZyg) = (getCrossRatio($refTable, $testTable, $chrom, $minPos, $maxPos), $normMean);
    }
    my $start = ${$$probesRef[$lowI]}{POSITION};
    my $end = ${$$probesRef[$highI]}{POSITION};
    my $size = $end - $start + 1; 
    my $copyNumber = getCopyNumber($testMean, $normMean);
    print $nmaFileH join(",", $nmaID, $nmaType, $nProbes, $nStdDev,
                              $chrom, $start, $end, $size, $n, 
                              $testMean, $testSD, $refMean, $refSD, 
                              $normRatio, $normZyg, $pValue,
                              $copyNumber)."\n";  
}     

sub getCrossZygosity{
    my ($refTable, $testTable, $chrom, $minPos, $maxPos) = @_;
    runSQL("SELECT nvl(Log(2, Avg(Abs(t.BFREQUENCY - 0.5) + 0.5)/Avg(Abs(r.BFREQUENCY - 0.5) + 0.5)), 0)
                  FROM $refTable r, $testTable t
                  WHERE r.CHROMOSOME = t.CHROMOSOME AND r.POSITION = t.POSITION
                    AND r.CHROMOSOME = $chrom AND r.POSITION >= $minPos and r.POSITION <= $maxPos
                    AND r.RATIO > 0.333 AND r.RATIO < 3.0
                    AND r.BFREQUENCY > 0.25 AND r.BFREQUENCY < 0.75 ", \my($zygosity));
    fetchRow();
    return $zygosity;
}

sub getCrossRatio{
    my ($refTable, $testTable, $chrom, $minPos, $maxPos) = @_;
    runSQL("SELECT nvl(Log(2, nullif(Avg(t.RATIO),0)/nullif(Avg(r.RATIO),0)), 0)
                  FROM $refTable r, $testTable t
                  WHERE r.CHROMOSOME = t.CHROMOSOME AND r.POSITION = t.POSITION
                    AND r.CHROMOSOME = $chrom AND r.POSITION >= $minPos and r.POSITION <= $maxPos 
                    AND r.RATIO > 0.333 AND r.RATIO < 3.0", \my($log2R));
    fetchRow();
    return $log2R;
}

sub trimCNVSetEnds {
    my ($probesRef, $minI, $maxI, $lowI, $highI, $threshold) = @_;
    $lowI -= $pf;    
    $lowI >= $minI or $lowI = $minI;
    $highI += $pf;
    $highI <= $maxI or $highI = $maxI; 
    my $n = ($highI - $lowI + 1);
    
#    status("\nIs: $lowI, $highI\n");
    
    my ($tL2Mean, $tL2StDev) = getEventStats($probesRef, $lowI, $highI, $n, 'TL2'); 
    
#    status("test stats: $tL2Mean, $tL2StDev\n");
#    exit;
    
    my ($tL2Min, $tL2Max) = ($tL2Mean - (2 * $tL2StDev),  $tL2Mean + (2 * $tL2StDev));    
    my $tL2Step = ($tL2Max - $tL2Min)/100;    
    my %candCounts;
    for (my $candMean = $tL2Min; $candMean <= $tL2Max; $candMean += $tL2Step){
        my ($candMin, $candMax) = ($candMean - (2 * $testL2StDev),  $candMean + (2 * $testL2StDev));  
        my ($candCount, $candMinI, $candMaxI) = getCandCount($probesRef, $lowI, $highI, $candMin, $candMax, 'TL2');
        push @{$candCounts{$candCount}}, [$candMean, $candMinI, $candMaxI];
    }
    my @candCounts = sort {$b <=> $a} keys %candCounts; 
    my $bestCandCount = $candCounts[0];
    my %deltas;
    foreach my $candRef (@{$candCounts{$bestCandCount}}){
        my ($delta, $rL2Mean) = getCandDelta($probesRef, @$candRef);
        $delta >= $threshold or next;  
        push @{$deltas{$delta}}, [@$candRef, $rL2Mean];
    }
    my @deltas = sort {$b <=> $a} keys %deltas; 
    scalar @deltas or return undef;
    my $bestDelta = $deltas[0];
    ($tL2Mean, my($j, $k, $rL2Mean)) = @{${$deltas{$bestDelta}}[0]};
    while ($j > $lowI and abs(${$$probesRef[$j - 1]}{TL2} - $tL2Mean) < abs(${$$probesRef[$j - 1]}{TL2} - $rL2Mean)){ $j-- }
    while ($k < $highI and abs(${$$probesRef[$k + 1]}{TL2} - $tL2Mean) < abs(${$$probesRef[$k + 1]}{TL2} - $rL2Mean)){ $k++ }
    return ($j, $k);
}

sub getEventStats{
    my ($probesRef, $lowI, $highI, $n, $field) = @_;
    my ($sum, $sumDev2);
    foreach my $i($lowI..$highI){ $sum += ${$$probesRef[$i]}{$field} }
    my $mean = $sum/$n;
    foreach my $i($lowI..$highI){ $sumDev2 += ((${$$probesRef[$i]}{$field} - $mean) ** 2) }
    my $stdDev = sqrt($sumDev2 / ($n - 1));
    return ($mean, $stdDev); 
}

sub getCandCount{
    my ($probesRef, $lowI, $highI, $min, $max, $field) = @_;
    my ($count, $minI, $maxI) = (0);
    
    #this will overinclude adjacent baseline regions from lowI to highI
    #however, the mean will be more accurate
    #no! test mean is getting recalculated!
    
    foreach my $i($lowI..$highI){ 
        my $value = ${$$probesRef[$i]}{$field};
        if ($value >= $min and $value <= $max){
            $minI or $minI = $i;
            $maxI = $i;
            $count++;
        }  
    }
    return ($count, $minI, $maxI);
}

sub getCandDelta{
    my ($probesRef, $candMean, $lowI, $highI) = @_;
    my $n = ($highI - $lowI + 1);
    $n > 1 or return 0;
    my ($rL2Mean, $rL2StDev) = getEventStats($probesRef, $lowI, $highI, $n, 'RL2'); 
    $candMean == 0 and $candMean = 0.001;
    $rL2Mean == 0 and $rL2Mean = 0.001;
    return (abs(logBase2((2 ** $candMean) / (2 ** $rL2Mean)) - $nmaMean), $rL2Mean);
}





sub trimCNVSetEndsXXX {
    my ($probesRef, $minI, $maxI, $lowI, $highI) = @_;
    my $minJ = $lowI - $pf;
    $minJ >= $minI or $minJ = $minI;
    my $maxJ = $lowI + $pf;
    $maxJ <= $maxI or $maxJ = $maxI;
    my $minK = $highI - $pf;
    $minK >= $minI or $minK = $minI;
    my $maxK = $highI + $pf;
    $maxK <= $maxI or $maxK = $maxI;
    my %normMeans;
     foreach my $j($minJ..$maxJ){
        foreach my $k($minK..$maxK){   
            my $n = ($k - $j + 1);
            $n > 1 or next;        
            my ($testMean, $refMean) = getSampleMeans($probesRef, $j, $k, $n); 
            my $normMean = abs(logBase2($testMean / $refMean));             
            $normMeans{$normMean} = [$j, $k, $n, $testMean, $refMean];    
        }
    }    
    my @normMeans = sort {$b <=> $a} keys %normMeans; 
    my $bestNormMean = $normMeans[0];      
    my ($j, $k, $n, $testMean, $refMean) = @{$normMeans{$bestNormMean}}; 
    #################################################################################
    #decide for sure which expansion elaboration is better, by test SD or proximity to means
    #################################################################################
    #my ($discard1, $testSD) = getSampleStdDev($probesRef, $j, $k, $n, 'TVAL', $testMean);
    #my ($discard2, $refSD) =  getSampleStdDev($probesRef, $j, $k, $n, 'RVAL', $refMean);
    #my $test2SD = $testSD * 2;
    #while ($j > $minJ and abs(${$$probesRef[$j - 1]}{TVAL} - $testMean) <= $test2SD){ $j-- }
    #while ($k < $maxK and abs(${$$probesRef[$k + 1]}{TVAL} - $testMean) <= $test2SD){ $k++ }    
    while ($j > $minJ and abs(${$$probesRef[$j - 1]}{TVAL} - $testMean) < abs(${$$probesRef[$j - 1]}{TVAL} - $refMean)){ $j-- }
    while ($k < $maxK and abs(${$$probesRef[$k + 1]}{TVAL} - $testMean) < abs(${$$probesRef[$k + 1]}{TVAL} - $refMean)){ $k++ }
    return ($j, $k);   
}

sub getCNVSetStats{
    my ($probesRef, $lowI, $highI) = @_;
    my $n = ($highI - $lowI + 1);
    $n > 1 or return 0;
    my ($testMean, $refMean) = getSampleMeans($probesRef, $lowI, $highI, $n); 
    my ($testSD2, $testSD) = getSampleStdDev($probesRef, $lowI, $highI, $n, 'TVAL', $testMean);
    my ($refSD2, $refSD) =   getSampleStdDev($probesRef, $lowI, $highI, $n, 'RVAL', $refMean);    
    my $pValue = runTTest($n, $testMean, $refMean, $testSD2, $refSD2);
    return ($n, $testMean, $refMean, $testSD, $refSD, $pValue);
}

sub getSampleMeans{
    my ($probesRef, $lowI, $highI, $n) = @_;
    my ($testSum, $refSum);
    foreach my $i($lowI..$highI){
        $testSum += ${$$probesRef[$i]}{TVAL};
        $refSum += ${$$probesRef[$i]}{RVAL}; 
    }
    return ($testSum/$n, $refSum/$n);
}

sub getSampleStdDev{
    my ($probesRef, $lowI, $highI, $n, $field, $mean) = @_;
    my $sumDev2;    
    foreach my $i($lowI..$highI){ $sumDev2 += ((${$$probesRef[$i]}{$field} - $mean) ** 2) }
    my $stdDev2 = $sumDev2 / ($n - 1);
    return ($stdDev2, sqrt $stdDev2);    
}

sub runTTest{
    my ($n, $testMean, $refMean, $testSD2, $refSD2) = @_;
    my $tNum = abs($testMean - $refMean);
    my $tDenom = sqrt(($testSD2/$n) + ($refSD2/$n));
    $tDenom or return 0;
    my $tValue = $tNum / $tDenom;
    #adjust this if tTable is made 2 sig digits
    my $tScalar = 10 ** (length(int($tValue)) - 1);
    $tValue = int($tValue / $tScalar) * $tScalar;    
    my $dof = $n - 1 ; 
    #change this is tTable is expanded to larger N    
    $dof > 200 and $dof = 200;   
    my $pValue;
    if(defined $tTable{$dof}{$tValue}){
        $pValue = $tTable{$dof}{$tValue} 
    } else {
        $pValue = 0 #$tValue was very large...
    }
    return $pValue;
}

sub initializeTTable{
    my $tTableFile = "$param{vampPath}/bin/tTable.csv";
    open my $tTableFileH, "<", $tTableFile;
    my $tValues = <$tTableFileH>;
    my ($emptyCell, @tValues) = split(",", $tValues);
    my $maxI = (scalar @tValues) - 1;
    while (<$tTableFileH>){
        my ($dof, @pValues) = split(",", $_);
        foreach my $i(0..$maxI){ $tTable{$dof}{$tValues[$i]} = $pValues[$i] }
    }
    close $tTableFileH;
}

sub logBase2 {
    my ($value) = @_;
    $value or return log(0.0001)/log(2); 
    return log($value)/log(2);
}

sub getCopyNumber{
    my ($testMean, $normMean) = @_;
    if($nmaType eq 'RATIO'){
        $normMean < -1 and return 0;
        $normMean < -0.2 and return 1;
        $normMean < 0.2 and return 2;
        return 3;
    } else { #BFREQUENCY, $testMean = zygosity score
        $testMean > 0.8 and return 1; #very high zygosity = copy number loss, i.e deletion to 1
        my $copyNumber = 1 / (1 - $testMean);
        return int($copyNumber + 0.5);
    }
}

sub findCNVBlocks{
    my ($nmaTable) = @_;
    status("finding overlapping blocks of CNV sets...\n");
    foreach my $nStdDev(@nStdDevs){
        status("  nStdDev = $nStdDev: Chr "); 
        foreach my $chrom (1..$refSeqs{$param{refSeqBase}}{nChrom}){
            status("$chrom ");
            runSQL("SELECT NMATYPE, START_, END_, NPROBESINSET, COPYNUMBER, NORMALIZEDRATIO, NORMALIZEDZYGOSITY, PVALUE
                    FROM $nmaTable 
                    WHERE CHROMOSOME = $chrom
                      AND NSTDDEV = $nStdDev
                    ORDER BY START_", \my($type, $start, $end, $nProbesInSet, $copyNumber, $normRatio, $normZyg, $pValue));  
            fetchRow() or next;                    
            my ($blockStart, $blockEnd, $maxNPIS,      $copyNumberSum, $bestNormRatio, $bestNormZyg, $bestPValue, $nEvents) = 
               ($start,      $end,      $nProbesInSet, $copyNumber,    $normRatio,     $normZyg,     $pValue,         1);
            while (fetchRow()){
                if($start > $blockEnd){
                    $nmaID++;
                    my $blockSize = $blockEnd - $blockStart + 1;
                    my $blockCopyNumber = int(($copyNumberSum / $nEvents) + 0.5);
                    print $nmaFileH join(",", $nmaID, 'BLOCK', 0, $nStdDev,
                                              $chrom, $blockStart, $blockEnd, $blockSize, $maxNPIS, 
                                              0, 0, 0, 0,
                                              $bestNormRatio, $bestNormZyg, $bestPValue,
                                              $blockCopyNumber)."\n"; 
                    ($blockStart, $blockEnd, $maxNPIS,      $copyNumberSum, $bestNormRatio, $bestNormZyg, $bestPValue, $nEvents) = 
                    ($start,      $end,      $nProbesInSet, $copyNumber,    $normRatio,     $normZyg,     $pValue,         1);
                } else {
                    $blockEnd > $end or $blockEnd = $end; 
                    $maxNPIS > $nProbesInSet or $maxNPIS = $nProbesInSet; 
                    $bestPValue < $pValue or $bestPValue = $pValue;
                    $copyNumberSum += $copyNumber;
                    abs($bestNormRatio) > abs($normRatio) or $bestNormRatio = $normRatio;
                    abs($bestNormZyg) > abs($normZyg) or $bestNormZyg = $normZyg;
                    $nEvents++;
                }    
            }
            $nmaID++;
            my $blockSize = $blockEnd - $blockStart + 1;
            my $blockCopyNumber = int(($copyNumberSum / $nEvents) + 0.5);
            print $nmaFileH join(",", $nmaID, 'BLOCK', 0, $nStdDev,
                                      $chrom, $blockStart, $blockEnd, $blockSize, $maxNPIS, 
                                      0, 0, 0, 0,
                                      $bestNormRatio, $bestNormZyg, $bestPValue,
                                      $blockCopyNumber)."\n"; 
        }
        status("\n");
    } 
}

#    #my %normMeans;
#    my %pValues;
#    foreach my $j($minJ..$maxJ){
#        foreach my $k($minK..$maxK){
##            my $n = ($k - $j + 1);
##            $n > 1 or next;
#            my ($n, $testMean, $refMean, $testSD, $refSD, $pValue) = getCNVSetStats($probesRef, $j, $k);
#            $n or next;
#            #my ($testMean, $refMean) = getSampleMeans($probesRef, $j, $k, $n); 
#            #my $normMean = abs(logBase2($testMean / $refMean));             
#            #$normMeans{$normMean} = [$j, $k, $n, $testMean];
#            push @{$pValues{$pValue}}, [$j, $k, $testMean, $refMean];
#        }
#    }
#    my @pValues = sort {$a <=> $b} keys %pValues;
#    my $bestPValue = $pValues[0];
#    my %normMeans;
#    foreach my $arrayRef(@{$pValues{$bestPValue}}){ 
#        my ($j, $k, $testMean, $refMean) = @$arrayRef;
#        my $normMean = abs(logBase2($testMean / $refMean));  
#        push @{$normMeans{$normMean}}, [$j, $k];
#    }    
#    my @normMeans = sort {$b <=> $a} keys %normMeans; 
#    my $bestNormMean = $normMeans[0];   
#    my ($bestJ, $bestK) = (1E9, 0);
#    foreach my $arrayRef(@{$normMeans{$bestNormMean}}){ 
#        my ($j, $k) = @$arrayRef;
#        $bestJ <= $j or $bestJ = $j;
#        $bestK >= $k or $bestK = $k;
#    }
#    return ($bestJ, $bestK);

#    $lowI = getBestLeftEnd($probesRef, $minI, $maxI, $lowI, $highI);
#    ($highI, $n, $testMean, $refMean, $testSD, $refSD, $pValue) = getBestRightEnd($probesRef, $minI, $maxI, $lowI, $highI);
        
##    $lowI -= $pf;
##    $lowI >= $minI or $lowI = $minI;
##    $highI += $pf;
##    $highI <= $maxI or $highI = $maxI;
##    my ($testMean, $refMean);    
##    ($lowI, $highI, $testMean, $refMean) = trimCNVSetEnds($probesRef, $lowI, $highI);    
#    $lowI or return;
#    #could purge internal outliers here... 
#    my ($testSD, $refSD, $pValue) = getCNVStats($probesRef, $lowI, $highI, $nProbesInSet, $testMean, $refMean);

#sub getBestLeftEnd{
#    my ($probesRef, $minI, $maxI, $lowI, $highI) = @_;
#    my $minJ = $lowI - (2 * $pf);
#    $minJ >= $minI or $minJ = $minI;
#    my $maxJ = $lowI + $pf;
#    $maxJ <= $maxI or $maxJ = $maxI;
#    my %pValues;
#    foreach my $j($minJ..$maxJ){
#        my ($n, $testMean, $refMean, $testSD, $refSD, $pValue) = getCNVSetStats($probesRef, $j, $highI);
#        push @{$pValues{$pValue}}, [$j, [$n, $testMean, $refMean, $testSD, $refSD, $pValue]];
#    }
#    my @pValues = sort {$a <=> $b} keys %pValues;
#    my $bestPValue = $pValues[0];
#    my %Js;
#    foreach my $jRef(@{$pValues{$bestPValue}}){ $Js{$$jRef[0]} = $$jRef[1] }
#    my @Js = sort {$a <=> $b} keys %Js;
#    my $bestJ = $Js[0];
#    return $bestJ, @{$Js{$bestJ}}
#}

#sub getBestLeftEnd{
#    my ($probesRef, $minI, $maxI, $lowI, $highI) = @_;
#    my $minJ = $lowI - (2 * $pf);
#    $minJ >= $minI or $minJ = $minI;
#    my $maxJ = $lowI + $pf;
#    $maxJ <= $maxI or $maxJ = $maxI;
#    my %pValues;
#    foreach my $j($minJ..$maxJ){
#        my ($n, $testMean, $refMean, $testSD, $refSD, $pValue) = getCNVSetStats($probesRef, $j, $highI);
#        $n or next;
#        push @{$pValues{$pValue}}, $j;
#    }
#    my @pValues = sort {$a <=> $b} keys %pValues;
#    my $bestPValue = $pValues[0];
#    my @Js = sort {$a <=> $b} @{$pValues{$bestPValue}};
#    my $bestJ = $Js[0];
#    return $bestJ;   
#}

#sub getBestRightEnd{
#    my ($probesRef, $minI, $maxI, $lowI, $highI) = @_;
#    my $minK = $highI - $pf;    
#    $minK >= $minI or $minK = $minI;    
#    my $maxK = $highI + (2 * $pf);    
#    $maxK <= $maxI or $maxK = $maxI; 
#    my %pValues;
#    foreach my $k($minK..$maxK){   
#        my ($n, $testMean, $refMean, $testSD, $refSD, $pValue) = getCNVSetStats($probesRef, $lowI, $k);    
#        $n or next;
#        push @{$pValues{$pValue}}, [$k, [$n, $testMean, $refMean, $testSD, $refSD, $pValue]];
#    } 
#    my @pValues = sort {$a <=> $b} keys %pValues;
#    my $bestPValue = $pValues[0];
#    my %Ks;
#    foreach my $kRef(@{$pValues{$bestPValue}}){ $Ks{$$kRef[0]} = $$kRef[1] }   
#    my @Ks = sort {$b <=> $a} keys %Ks;     
#    my $bestK = $Ks[0];
#    return $bestK, @{$Ks{$bestK}}
#}


#sub getCNVStats{ #by T Test
#    my ($probesRef, $lowI, $highI, $nProbesInSet, $testMean, $refMean) = @_;
##    my ($testSD2, $testSD) = getSampleStdDev($probesRef, $lowI, $highI, $nProbesInSet, 'TVAL', $testMean);
##    my ($refSD2, $refSD) =   getSampleStdDev($probesRef, $lowI, $highI, $nProbesInSet, 'RVAL', $refMean);
#    my $tNum = abs($testMean - $refMean);
#    my $tDenom = sqrt(($testSD2/$nProbesInSet) + ($refSD2/$nProbesInSet));
#    my $tValue = $tNum / $tDenom;
#    my $tScalar = 10 ** (length(int($tValue)) - 1);
#    $tValue = int($tValue / $tScalar) * $tScalar;    
#    my $dof = $nProbesInSet - 1 ;     
#    $dof > 200 and $dof = 200;   
#    my $pValue;
#    if(defined $tTable{$dof}{$tValue}){
#        $pValue = $tTable{$dof}{$tValue} 
#    } else {
#        $pValue = 999 #an error flag
#    }
#    return ($testSD, $refSD, $pValue);
#}

#sub trimCNVSetEnds{
#    my ($probesRef, $lowI, $highI) = @_;
#    my ($testMean, $refMean) = getSampleMeans($probesRef, $lowI, $highI);   
#    #TO do:  if the terminal probe passes, should try adding to the event
#    #i.e. should expand as well as contract??
#    #when expanding might possibly want to consider whether the new probe drops the SD
#    while(isBaselineProbe(${$$probesRef[$lowI]}{TVAL}, $testMean, $refMean)){
#        $lowI++;
#        $lowI > $highI - $nProbes + 1 and return undef;        
#        ($testMean, $refMean) = getSampleMeans($probesRef, $lowI, $highI);    
#    }  
#    while(isBaselineProbe(${$$probesRef[$highI]}{TVAL}, $testMean, $refMean)){
#        $highI--;
#        $highI < $lowI + $nProbes - 1 and return undef;               
#        ($testMean, $refMean) = getSampleMeans($probesRef, $lowI, $highI);    
#    } 
#    return ($lowI, $highI, $testMean, $refMean);
#}

#sub isBaselineProbe{
#    my ($testValue, $testMean, $refMean) = @_;
#    return (abs($testValue - $testMean) > abs($testValue - $refMean));
#}

#sub checkWithOutliersRemoved{ #this sub purges false CNV probe sets caused by extreme outlier values
#    my ($probesRef, $lowI, $highI, $indexNMA, $mean, $stDev, $nStdDev)= @_;
#    my ($nProbes, $isPos, %normVals, @normVals, $TValSum, $RValSum, $nKeptProbes) = (($highI - $lowI + 1), ($indexNMA > 0));
#    my $normValsToKeep = int(($nProbes * 0.8) + 0.5); #discard top 20% of values  
#    foreach my $i($lowI..$highI){
#        my $normVal = logBase2(${$$probesRef[$i]}{TVAL} / ${$$probesRef[$i]}{RVAL});        
#        push @{$normVals{$normVal}}, $$probesRef[$i];  
#    }
#    if ($isPos){
#        @normVals = sort {$a <=> $b} keys %normVals;
#    } else {
#        @normVals = sort {$b <=> $a} keys %normVals;
#    }   
#    @normVals = splice(@normVals, 0, $normValsToKeep);
#    foreach my $normVal(@normVals){
#        my $probesRef = $normVals{$normVal};
#        foreach my $probeRef(@$probesRef){
#            $TValSum += $$probeRef{TVAL};
#            $RValSum += $$probeRef{RVAL};
#            $nKeptProbes++;
#        }
#    }
#    my ($TMean, $RMean) = ($TValSum/$nKeptProbes, $RValSum/$nKeptProbes);
#    $RMean == 0 and $RMean = 0.001;
#    my $normMean = logBase2($TMean / $RMean); 
#    my $threshold = $stDev * ($nStdDev - 1);  #check against a more permissive criterion  
#    return (abs($normMean - $mean) >= $threshold);
#}


1;





















