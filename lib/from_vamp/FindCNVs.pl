#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#microarray analysis
#########################################################################

use vars(qw(%param %types %fields %fieldNames %refSeqs %reverseRefSeqs));  #common parameters used by most analysis scripts

my ($testSample, $nStdDev, @nStdDevs, $nProbes, $pf, $nmaFileH, $nmaID, $nmaType);
my $scalar = 0.001;
my ($refSelectSQL, $testSelectSQL, $nValSQL, $maSQL, $nmaSQL, $nmaTable, $nmaDataTable, $nvlTable);
my ($nmaMean, $nmaStDev, $nvalMean, $nvalStDev, %zTable);  
my ($refTable, $testTable);
my $posID = " (CHROMOSOME || 'x' || POSITION) "; 
my $noZeroRatio = " (CASE RATIO WHEN 0 THEN $scalar ELSE RATIO END) "; #prevent divide by and log of zero errors
my $noZeroR = " (CASE R WHEN 0 THEN $scalar ELSE R END) "; #prevent divide by and log of zero errors
my %snpTypes = (Illumina => 1, MatePair => 1, Affymetrix => 1);
my %imputeLOH = (Illumina => 1, MatePair => 1); #could allow imputation on Affy if can define a usable LOH rule
my $isSNPs = 0;
my ($minChrom, $maxChrom, $isHemizygous);
my @nProbes = ();
my %tCritical = ();

sub importArray{
    my ($run) = @_;
    if($param{arrayType} eq 'Illumina'){
        importIllumina($run)
    } elsif ($param{arrayType} eq 'Nimblegen') {
        importNimblegen($run)
    } elsif ($param{arrayType} eq 'Affymetrix') {
        importAffy($run)
    }
}

sub importExternalCNVs{
    my ($run) = @_;
    if($param{arrayType} eq 'Illumina'){
        importIlluminaCNVs($run)
#    } elsif ($param{arrayType} eq 'Nimblegen') {
#        importNimblegenCNVs($run)
    }
}

sub arrayTableName{
    my ($sample, $arrayType, $refSeq) = @_;
    my %shortTypes = (Nimblegen => 'NIM', Illumina => 'ILL', MatePair => 'MP', Affymetrix => 'AFFY');
    return "$sample\_$shortTypes{$arrayType}\_$refSeq";
}

sub nmaTableName{
    my ($testSample) = @_;
    if ($param{arrayRef}) {
        return "$param{arrayRef}\_$testSample";
    } else {
        return $testSample;
    }
}

sub analyzeArray{
    ($testSample) = @_;
    $testSample or die "test array sample must be specified for command analyzeArray\n";
    initializeTCritical();
    $snpTypes{$param{arrayType}} and $isSNPs = 1;
    $isSNPs or $param{arrayRef} = undef;
    @nProbes = split(",", $param{nProbes});
    @nProbes = sort {$a <=> $b} @nProbes;
    @nStdDevs = split(",", $param{nStdDevs});
    $nmaTable = newTable('NMA', nmaTableName($testSample)); 
    $nvlTable = "\U$nmaTable\_j";
    $nmaDataTable = "\U$nmaTable\_d";   
    my $nmaFile = "$nmaTable.csv"; 
    open $nmaFileH, ">", $nmaFile; 
    $fieldNames{NMA} =~ s/, DESCRIPTION//;
    getSexChromBlocks(\my@sexChromBlocks);    
    foreach my $sexChromBlockRef(@sexChromBlocks){   
        ($minChrom, $maxChrom, $isHemizygous) = @$sexChromBlockRef;    
        status("\nanalyzing chromosomes $reverseRefSeqs{$param{refSeqBase}}{$minChrom} to $reverseRefSeqs{$param{refSeqBase}}{$maxChrom}\n");
        (!$param{arrayRef} and $imputeLOH{$param{arrayType}}) and imputeLOH();
        NMATYPE: foreach my $nmaType_('RATIO', 'BFREQUENCY'){  
            $nmaType = $nmaType_;   
            if($nmaType eq 'BFREQUENCY'){($isSNPs and !$isHemizygous) or next NMATYPE} 
            if($nmaType eq 'BFREQUENCY'){$param{arrayType} eq 'MatePair' and next NMATYPE}      
            status("nmaType = $nmaType\n");
            calculateNVals() or next NMATYPE; 
            #$nmaType eq 'RATIO' and parseExternalCNVs();
            NPROBES: foreach my $nProbes_(@nProbes){        
                $nProbes = $nProbes_;
                $nProbes % 2 or $nProbes++; 
                $pf = ($nProbes - 1) / 2;
                status("  nProbes = $nProbes\n");       
                calculateNMAs() or next NPROBES;
                findCNVSets();
                dropTable($nmaDataTable);    
            }
            dropTable($nvlTable);    
        }
    }
    close $nmaFileH;
    loadData($nmaFile, $nmaTable, ",", $fieldNames{NMA});
    open $nmaFileH, ">", $nmaFile;
    findCNVBlocks($nmaTable);
    close $nmaFileH;
    loadData($nmaFile, $nmaTable, ",", $fieldNames{NMA});
}     

sub getSexChromBlocks{
    my ($sexChromBlocksRef) = @_;
    
    # this need to become aware of reference sex if different than test sex!
    
    if($param{sex}){
        if($param{sex} eq 'XX'){
            push @$sexChromBlocksRef, [1, $refSeqs{$param{refSeqBase}}{nChrom} - 1, 0]; #autosomes and disomic X
        } elsif ($param{sex} eq 'XY'){
            push @$sexChromBlocksRef, [1, $refSeqs{$param{refSeqBase}}{nChrom} - 2, 0]; #autosomes
            push @$sexChromBlocksRef, [$refSeqs{$param{refSeqBase}}{nChrom} - 1, $refSeqs{$param{refSeqBase}}{nChrom}, 1]; #monosomic X and Y
        } else {
            die "sex $param{sex} not supported"
        }
    } else {
        push @$sexChromBlocksRef, [1, $refSeqs{$param{refSeqBase}}{nChrom}, 0]; #no sex specified, use all chromosomes
    }
}

sub imputeLOH{
    $isHemizygous and return;
    status("imputing LOH...\n"); 
    my $afTable = getTableName("AF", "CEU");   
    $refTable = undef;
    $testTable = getTableName('Array', arrayTableName($testSample, $param{arrayType}, $param{refSeq}));   
    status("  Chr ");   
    foreach my $chrom ($minChrom..$maxChrom){
        status("$chrom ");
        runSQL("SELECT t.POSITION, t.ZYGOSITY, nvl(af.UNINFORMATIVEFREQ,1) UNINFORMATIVEFREQ
                    FROM $testTable t, $afTable af
                    WHERE t.CHROMOSOME = $chrom
                      AND t.CHROMOSOME = af.CHROMOSOME(+)
                      AND t.POSITION = af.POSITION(+)
                      AND t.GCSCORE >= 0.15 AND t.R > 0.1 AND t.R < 3.0
                    ORDER BY t.POSITION");
        my $probesRef = fetchAllHashRef();    
        my ($cumUnInfFreq, $nProbesInSet, $minPos, $prevPos, $eHet) = (1, 0, 0, 0, 0);           
        foreach my $__(@$probesRef){
            if($$__{ZYGOSITY} > 0.9){  #this rule only works for Illumina where BAF is easily binnable
                $minPos or $minPos = $$__{POSITION};
                $cumUnInfFreq *= $$__{UNINFORMATIVEFREQ};
                $nProbesInSet++; 
                $eHet += (1 - $$__{UNINFORMATIVEFREQ});
            } else {
                commitLOH($chrom, $cumUnInfFreq, $nProbesInSet, $minPos, $prevPos, $eHet);
                ($cumUnInfFreq, $nProbesInSet, $minPos, $eHet) = (1, 0, 0, 0);  
            }
            $prevPos = $$__{POSITION};
        }
        commitLOH($chrom, $cumUnInfFreq, $nProbesInSet, $minPos, $prevPos, $eHet);
    } 
    status("\n");
}

sub commitLOH{
    my ($chrom, $cumUnInfFreq, $nProbesInSet, $minPos, $maxPos, $eHet) = @_;
    $cumUnInfFreq < 0.001 or return;
    my ($log2R, $chi2, $size);    
    my @nProbes = split(",", $param{nProbes});      
    foreach my $nProbes_(sort {$a <=> $b} @nProbes){
        my $nProbes = $nProbes_;
        $nProbes % 2 or $nProbes++;  
        if($nProbesInSet >= $nProbes){  
            $nmaID++;   
            defined $log2R or $log2R = getCrossRatio($chrom, $minPos, $maxPos); 
            defined $chi2 or $chi2 = ($nProbesInSet * $eHet)/($nProbesInSet - $eHet); 
            defined $size or $size = $maxPos - $minPos + 1;  
            foreach my $nStdDev(@nStdDevs){      
                print $nmaFileH join(",", $nmaID, 'LOH', $nProbes, $nStdDev,
                                          $chrom, $minPos, $maxPos, $size, $nProbesInSet, 
                                          $cumUnInfFreq, 0, 0, 0, 
                                          $log2R, logBase2(1 / 0.5), $chi2, 1)."\n";                                                                                            
            }                
        }
    }
}

sub calculateNVals{
    status("  calculating normalized values...\n"); 
    ($refTable, $testTable) = (undef, undef);
    $param{arrayRef} and $refTable =  getTableName('Array', arrayTableName($param{arrayRef}, $param{arrayType}, $param{refSeq}));
    $testTable = getTableName('Array', arrayTableName($testSample, $param{arrayType}, $param{refSeq}));
    $refSelectSQL = getArraySelectSQL($param{arrayRef}, 1);
    $testSelectSQL = getArraySelectSQL($testSample, 0);  
    ($refSelectSQL and $testSelectSQL) or return undef;
    if($nmaType eq 'BFREQUENCY'){
        my $bfr = "(t.BF/r.BF)";
        my $bfrSQL = " SELECT $bfr BFR
                       FROM ($refSelectSQL) r, ($testSelectSQL) t
                       WHERE r.POSID = t.POSID "; 
        
        #how well do these limits work for WGA?               
        runSQL("SELECT Avg(BFR) FROM ($bfrSQL) WHERE BFR > 0.66 AND BFR < 1.34", \my($bfrMean));
        
        fetchRow();
        $bfr = "($bfr + 1 - $bfrMean)";    
        $nValSQL = " SELECT r.CHROMOSOME, r.POSITION, 
                            (Abs(t.BFREQUENCY - 0.5) + 0.5) TVAL, (Abs(r.BFREQUENCY - 0.5) + 0.5) RVAL, 
                            Log(2, (Abs($bfr - 1) + 1)) NVAL  
                     FROM ($refSelectSQL) r, ($testSelectSQL) t
                     WHERE r.POSID = t.POSID ";                             
    } else {
        $nValSQL = " SELECT r.CHROMOSOME, r.POSITION, 
                            t.RATIO TVAL, r.RATIO RVAL, 
                            Log(2, (t.R/r.R)) NVAL
                     FROM ($refSelectSQL) r, ($testSelectSQL) t
                     WHERE r.POSID = t.POSID ";  
    }  
    dropTable($nvlTable);    
    runSQL("CREATE TABLE $nvlTable AS $nValSQL");
    ($nvalMean, $nvalStDev) = getStatsFromL2("SELECT Round(NVAL/$scalar)*$scalar NVAL FROM $nvlTable",'NVAL');
    defined $nvalMean or return undef;
    status("    normalized value mean = $nvalMean, stDev = $nvalStDev\n");  
    $nmaID++;
    foreach my $chrom ($minChrom..$maxChrom){
        print $nmaFileH join(",", $nmaID, "$nmaType\_STATS_$chrom", 0, 0,
                                  0, 0, 0, 0, 0, 
                                  $nvalMean, $nvalStDev, 0, 0,
                                  0, 0, 0,
                                  -1)."\n";     
    }
    return 1;
}

sub calculateNMAs{
    status("  calculating normalized moving averages...\n"); 
    $maSQL =   " SELECT CHROMOSOME, POSITION, TVAL, RVAL, NVAL,  
                        Avg(NVAL) OVER (ORDER BY CHROMOSOME, POSITION ROWS BETWEEN $pf PRECEDING AND $pf FOLLOWING) NMA
                 FROM ($nvlTable) "; 
    $nmaSQL = " SELECT CHROMOSOME, POSITION, TVAL, RVAL, Round(NVAL/$scalar)*$scalar NVAL, Round(NMA/$scalar)*$scalar NMA  
                FROM ($maSQL) ";            
    dropTable($nmaDataTable);
    runSQL("CREATE TABLE $nmaDataTable AS $nmaSQL"); 
    ($nmaMean, $nmaStDev) = getStatsFromL2($nmaDataTable, 'NMA');
    defined $nmaMean or return undef;
    status("    normalized moving average mean = $nmaMean, stDev = $nmaStDev\n");  
    return 1;
}

sub getArraySelectSQL{
    my ($sample, $isRef) = @_;
    my $chromFilter = "  WHERE CHROMOSOME >= $minChrom and CHROMOSOME <= $maxChrom ";    
    my $refFilter = " AND GCSCORE >= 0.15 AND R > 0.1 AND R < 3.0 ";
    if ($sample){ #sample is testSample or sample is a defined reference (true only for SNP array pairs)
        my $arrayTable = getTableName('Array', arrayTableName($sample, $param{arrayType}, $param{refSeq}));
        if ($nmaType eq 'BFREQUENCY'){ 
            my $bf = ' BF ';
            $param{arrayRef} or $bf = ' BFREQUENCY ';
            my $filter = ''; 
            my $informativeFilter = getInformativeFilter($bf, "");
            $isRef and $filter = " $refFilter $informativeFilter ";  #restrict to inferred or provided informative probes
            return " SELECT $posID POSID, CHROMOSOME, POSITION, BFREQUENCY, $bf BF
                     FROM $arrayTable $chromFilter $filter"; 
        } elsif ($nmaType eq 'RATIO'){
            my $r = $noZeroR;
            $param{arrayRef} or $r = $noZeroRatio;
            my $filter = '';
            $isRef and $filter = $refFilter; 
            return " SELECT $posID POSID, CHROMOSOME, POSITION, $noZeroRatio RATIO, $r R
                     FROM $arrayTable $chromFilter $filter";                    
        } 
    } else {  #sample is a non-existent reference, return testSample probes positions with idealized reference values
        my $arrayTable = getTableName('Array', arrayTableName($testSample, $param{arrayType}, $param{refSeq}));
        if ($nmaType eq 'BFREQUENCY'){   
            my $informativeFilter = getInformativeFilter("BFREQUENCY", "");                   
            my $filter = " $refFilter $informativeFilter ";                     
            return " SELECT $posID POSID, CHROMOSOME, POSITION, 0.5 BFREQUENCY, 0.5 BF
                     FROM $arrayTable $chromFilter $filter ";             
        } elsif ($nmaType eq 'RATIO'){
            return " SELECT $posID POSID, CHROMOSOME, POSITION, 1 RATIO, 1 R
                     FROM $arrayTable $chromFilter"; 
        } 
    }
}

sub getInformativeFilter{
    my ($bf, $prefix) =@_;
    my $prefixBF = "$prefix$bf";
    my $informativeFilter = " AND $prefixBF > 0.25 AND $prefixBF < 0.75 "; 
    my $prefixPosID = " ($prefix"."CHROMOSOME || 'x' || $prefix"."POSITION) "; 
    $param{infPosTable} and $informativeFilter = " AND $prefixPosID IN (SELECT $posID FROM $param{infPosTable}) "; 
    return $informativeFilter;
}

sub getStatsFromL2{
    my ($table, $field) = @_;
    my $histSQL = " SELECT $field X, Count($field) Y FROM ($table) WHERE $field != 0 GROUP BY $field ";
    runSQL("SELECT Max(Y) FROM ($histSQL)", \my($maxY));
    fetchRow();  
    defined $maxY or return undef;
    my $minY = 0.05 * $maxY;
    $histSQL = " SELECT X, Y FROM ($histSQL) WHERE Y >= $minY ";   
    runSQL("SELECT Min(X), Max(X) FROM ($histSQL)", \my($minX, $maxX));
    fetchRow();
    (defined $minX and defined $maxX) or return undef;
    my ($mean, $stDev);
    if ($nmaType eq 'BFREQUENCY'){
        ($mean, $stDev) = (0, $maxX / 2.45); 
    } else {
        my $sGuess = ($maxX - $minX) / (2 * 2.45); #expectations of stdDev for perfect Gaussian given minY at 5% of maxY, mean = 0    
        my $sGuess3 = $sGuess * 3;
        runSQL("SELECT Avg($field), StdDev($field) 
                FROM ($table)
                WHERE $field <= $sGuess3 AND $field >= -$sGuess3", \($mean, $stDev));
        fetchRow();     
    }        
    $mean =  int(($mean  / $scalar) + 0.5) * $scalar;
    $stDev = int(($stDev / $scalar) + 0.5) * $scalar;  
    return ($mean, $stDev);
}

sub findCNVSets{
    status("    finding CNV sets...\n");
    foreach my $nStdDev_(@nStdDevs){
        $nStdDev = $nStdDev_;
        status("      nStdDev = $nStdDev: Chr "); 
        my $nmaThreshold = $nmaStDev * $nStdDev;      
        foreach my $chrom ($minChrom..$maxChrom){
            status("$chrom ");
            runSQL("SELECT POSITION, TVAL, RVAL, NVAL, NMA
                    FROM $nmaDataTable
                    WHERE CHROMOSOME = $chrom
                    ORDER BY POSITION");    
            my $probesRef = fetchAllHashRef();       
            my ($minI, $maxI, $indexNMA, $lowI, $highI) = (0, (scalar @$probesRef) - 1, 0);                    
            foreach my $i($minI..$maxI){
                my $nma = ${$$probesRef[$i]}{NMA};
                my $overThreshold = (abs($nma - $nmaMean) >= $nmaThreshold);
                if (!$overThreshold){
                    if ($indexNMA){
                        checkCNVSet($chrom, $probesRef, $minI, $maxI, $lowI, $i - 1);
                        $indexNMA = 0;
                    }                
                } else {
                    if ($indexNMA){
                        if (($nma/$indexNMA) < 0){
                            checkCNVSet($chrom, $probesRef, $minI, $maxI, $lowI, $i - 1);
                            $lowI = $i; 
                            $indexNMA = $nma;
                        }
                    } else {
                        $lowI = $i; #start new set
                        $indexNMA = $nma;
                    }  
                }                 
            }
            $indexNMA and checkCNVSet($chrom, $probesRef, $minI, $maxI, $lowI, $maxI);  
        }  
        status("\n");        
    }
}

sub checkCNVSet{
    my ($chrom, $probesRef, $minI, $maxI, $lowI, $highI) = @_;
    my ($j, $k) = trimCNVSetEnds($probesRef, $minI, $maxI, $lowI, $highI);   
    my ($minPos, $maxPos) = (${$$probesRef[$j]}{POSITION}, ${$$probesRef[$k]}{POSITION});
    (defined $minPos and defined $maxPos) or return;
    my ($n, $testMean, $testStdDev, $refMean, $refStdDev, $normMean, $normStdDev, $normMax, $normMin, $Z, %outliers, $maxDevNVal);
    while(checkForOutliers($n, $normMean, $normStdDev, $normMax, $normMin, \$maxDevNVal)){
        defined $maxDevNVal and $outliers{$maxDevNVal} = 1;        
        ($n, $testMean, $testStdDev, $refMean, $refStdDev, 
         $normMean, $normStdDev, $normMax, $normMin, $Z) = getCNVSetStats($probesRef, $j, $k, \%outliers);    
         $n or return;
    }
    while(defined $outliers{${$$probesRef[$j]}{NVAL}}){$j++}
    while(defined $outliers{${$$probesRef[$k]}{NVAL}}){$k--}
    $nmaID++;         
    my $size = $maxPos - $minPos + 1;    
    my $copyNumber = getCopyNumber($testMean, $normMean);
    my ($normRatio, $normZyg);
    if ($nmaType eq 'RATIO'){
        ($normRatio, $normZyg) = ($normMean, getCrossZygosity($chrom, $minPos, $maxPos));
    } else {
        ($normRatio, $normZyg) = (getCrossRatio($chrom, $minPos, $maxPos), $normMean);
    }
    print $nmaFileH join(",", $nmaID, $nmaType, $nProbes, $nStdDev,
                              $chrom, $minPos, $maxPos, $size, $n, 
                              $testMean, $testStdDev, $refMean, $refStdDev, 
                              $normRatio, $normZyg, $Z,
                              $copyNumber)."\n";                        
}     

sub trimCNVSetEnds {
    my ($probesRef, $minI, $maxI, $lowI, $highI) = @_;
    my $j = $lowI + $pf; 
    my $k = $highI - $pf;    
    $j <= $k or ($j, $k) = ($k, $j);    
    $j < $minI and $j = $minI;    
    $k > $maxI and $k = $maxI;    
    $lowI -= $pf;  
    $lowI < $minI and $lowI = $minI;
    $highI += $pf;
    $highI > $maxI and $highI = $maxI; 
    my $nValThreshold = 2 * $nvalStDev;  
    while ($j > $lowI and $j >= ($minI + 2) and
           (abs(${$$probesRef[$j - 1]}{NVAL} - $nvalMean) >= $nValThreshold or
            abs(${$$probesRef[$j - 2]}{NVAL} - $nvalMean) >= $nValThreshold)){ $j-- }
    while ($k < $highI and $k <= ($maxI - 2) and 
           (abs(${$$probesRef[$k + 1]}{NVAL} - $nvalMean) >= $nValThreshold or
            abs(${$$probesRef[$k + 2]}{NVAL} - $nvalMean) >= $nValThreshold)){ $k++ }
    return ($j, $k);
}

sub getCNVSetStats{
    my ($probesRef, $lowI, $highI, $outliersRef) = @_;
    my $nOutliers = scalar(keys %$outliersRef);    
    my $n = ($highI - $lowI + 1) - $nOutliers;
    $n >= $nProbes or return 0;  
    my ($testMean, $testStdDev) = getEventStats($probesRef, $lowI, $highI, $n, 'TVAL', $outliersRef); 
    my ($refMean, $refStdDev) =   getEventStats($probesRef, $lowI, $highI, $n, 'RVAL', $outliersRef);
    my ($normMean, $normStdDev, $normMax, $normMin) = getEventStats($probesRef, $lowI, $highI, $n, 'NVAL', $outliersRef);
    my $Z = abs($normMean - $nvalMean) * sqrt($n) / $nvalStDev;  
    return ($n, $testMean, $testStdDev, $refMean, $refStdDev, $normMean, $normStdDev, $normMax, $normMin, $Z);
}

sub getEventStats{
    my ($probesRef, $lowI, $highI, $n, $field, $outliersRef) = @_;
    my ($max, $min, $sum, $sumDev2) = (-1E9, 1E9);
    foreach my $i($lowI..$highI){     
        defined $$outliersRef{${$$probesRef[$i]}{NVAL}} and next;
        my $value = ${$$probesRef[$i]}{$field};
        $sum += $value;
        $max >= $value or $max = $value;
        $min <= $value or $min = $value; 
    }
    my $mean = $sum/$n;
    foreach my $i($lowI..$highI){ $sumDev2 += ((${$$probesRef[$i]}{$field} - $mean) ** 2) }
    my $stdDev = sqrt($sumDev2 / ($n - 1));
    return ($mean, $stdDev, $max, $min); 
}

sub checkForOutliers{
    #using Grubb's test method
    my ($n, $normMean, $normStdDev, $normMax, $normMin, $maxDevNValRef) = @_;
    defined $n or return 1;  
    $normStdDev or return 0; 
    $$maxDevNValRef = $normMax;    
    $normMean < 0 and $$maxDevNValRef = $normMin;
    my $maxDev = abs($$maxDevNValRef - $normMean);
    my $G = $maxDev / $normStdDev;
    my $alpha = $param{outlierAlpha} / $n;  
    if ($alpha > 0.05){$alpha = 0.1}
    elsif ($alpha > 0.025){$alpha = 0.05}
    elsif ($alpha > 0.01){$alpha = 0.025}
    elsif ($alpha > 0.005){$alpha = 0.01}    
    elsif ($alpha > 0.001){$alpha = 0.005}    
    else {$alpha = 0.001}   
    my $dof = $n - 2;
    $dof > 100 and $dof = 999;
    my $t = $tCritical{$alpha}{$dof};
    my $t2 = $t ** 2;
    my $GLimit = ( ($n - 1) / sqrt($n) ) * sqrt( $t2 / ($n - 2 + $t2) );
    return ($G > $GLimit);
}

sub initializeTCritical{
    my $tCriticalFile = "$param{vampPath}/bin/tCritical.csv";
    open my $tCriticalFileH, "<", $tCriticalFile;
    my $line = <$tCriticalFileH>;
    chomp $line;
    $line =~ s/\r//g;
    my ($emptyCell, @alphas) = split(",", $line);
    my $maxI = (scalar @alphas) - 1;
    while (<$tCriticalFileH>){
        chomp $_;
        $_ =~ s/\r//g;
        my ($dof, @ts) = split(",", $_);
        foreach my $i (0..$maxI){ $tCritical{$alphas[$i]}{$dof} = $ts[$i] }
    }
    close $tCriticalFileH;
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

sub getCrossZygosity{
    $isHemizygous and return 0;
    my ($chrom, $minPos, $maxPos) = @_;
    my $zygosity = 0;
    if ($refTable){
        my $informativeFilter = getInformativeFilter("BFREQUENCY", "r.");  
        runSQL("SELECT nvl(Log(2, Avg(Abs(t.BFREQUENCY - 0.5) + 0.5)/Avg(Abs(r.BFREQUENCY - 0.5) + 0.5)), 0)
                FROM $refTable r, $testTable t
                WHERE r.CHROMOSOME = t.CHROMOSOME AND r.POSITION = t.POSITION
                  AND r.CHROMOSOME = $chrom AND r.POSITION >= $minPos and r.POSITION <= $maxPos
                  AND r.RATIO > 0.333 AND r.RATIO < 3.0
                  $informativeFilter ", \($zygosity));
                  
    } else {
        my $informativeFilter = getInformativeFilter("BFREQUENCY", "");  
        runSQL("SELECT nvl(Log(2, Avg(Abs(BFREQUENCY - 0.5) + 0.5)/0.5), 0)
                FROM $testTable
                WHERE CHROMOSOME = $chrom AND POSITION >= $minPos and POSITION <= $maxPos
                  AND RATIO > 0.333 AND RATIO < 3.0
                  $informativeFilter ", \($zygosity));     
                                   
    }
    fetchRow();
    return $zygosity;       
}

sub getCrossRatio{
    my ($chrom, $minPos, $maxPos) = @_;
    my $log2R = 0;
    if ($refTable){
        runSQL("SELECT nvl(Log(2, nullif(Avg(t.RATIO),0)/nullif(Avg(r.RATIO),0)), 0)
                FROM $refTable r, $testTable t
                WHERE r.CHROMOSOME = t.CHROMOSOME AND r.POSITION = t.POSITION
                  AND r.CHROMOSOME = $chrom AND r.POSITION >= $minPos and r.POSITION <= $maxPos 
                  AND r.RATIO > 0.333 AND r.RATIO < 3.0", \($log2R));   
    } else {
        runSQL("SELECT nvl(Log(2, Avg $noZeroRatio), 0)
                FROM $testTable
                WHERE CHROMOSOME = $chrom AND POSITION >= $minPos and POSITION <= $maxPos
                  AND GCSCORE >= 0.15 AND RATIO > 0.333 AND RATIO < 3.0 ", \($log2R));         
    }
    fetchRow();
    return $log2R; 
}

sub findCNVBlocks{
    my ($nmaTable) = @_;
    status("finding overlapping blocks of CNV sets...\n");
    foreach my $nStdDev(@nStdDevs){
        status("  nStdDev = $nStdDev: Chr "); 
        foreach my $chrom (1..$refSeqs{$param{refSeqBase}}{nChrom}){
            status("$chrom ");
            runSQL("SELECT NMATYPE, START_, END_, NPROBESINSET, COPYNUMBER, NORMALIZEDRATIO, NORMALIZEDZYGOSITY, ZSCORE
                    FROM $nmaTable 
                    WHERE CHROMOSOME = $chrom
                      AND NSTDDEV = $nStdDev
                      AND NMATYPE <> 'LOH'
                    ORDER BY START_", \my($type, $start, $end, $nProbesInSet, $copyNumber, $normRatio, $normZyg, $zScore));  
            fetchRow() or next;                    
            my ($blockStart, $blockEnd, $maxNPIS,      $copyNumberSum, $bestNormRatio, $bestNormZyg, $bestZScore, $nEvents) = 
               ($start,      $end,      $nProbesInSet, $copyNumber,    $normRatio,     $normZyg,     $zScore,         1);
            while (fetchRow()){
                if($start > $blockEnd){
                    $nmaID++;
                    my $blockSize = $blockEnd - $blockStart + 1;
                    my $blockCopyNumber = int(($copyNumberSum / $nEvents) + 0.5);
                    print $nmaFileH join(",", $nmaID, 'BLOCK', 0, $nStdDev,
                                              $chrom, $blockStart, $blockEnd, $blockSize, $maxNPIS, 
                                              0, 0, 0, 0,
                                              $bestNormRatio, $bestNormZyg, $bestZScore,
                                              $blockCopyNumber)."\n"; 
                    ($blockStart, $blockEnd, $maxNPIS,      $copyNumberSum, $bestNormRatio, $bestNormZyg, $bestZScore, $nEvents) = 
                    ($start,      $end,      $nProbesInSet, $copyNumber,    $normRatio,     $normZyg,     $zScore,         1);
                } else {
                    $blockEnd > $end or $blockEnd = $end; 
                    $maxNPIS > $nProbesInSet or $maxNPIS = $nProbesInSet; 
                    $bestZScore > $zScore or $bestZScore = $zScore;
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
                                      $bestNormRatio, $bestNormZyg, $bestZScore,
                                      $blockCopyNumber)."\n"; 
        }
        status("\n");
    } 
}

sub parseExternalCNVs{
    $param{arrayRef} and return;
    status("parsing external CNVs in $testSample array table...\n");
    # need to know expected CNV value??
    $nProbes = $nProbes[0];
    status("  Chr "); 
    foreach my $chrom ($minChrom..$maxChrom){
        status("$chrom ");
        runSQL("SELECT POSITION, CNVVALUE, RATIO TVAL, 1 RVAL, RATIO NVAL
                FROM $testTable 
                WHERE CHROMOSOME = $chrom
                  AND CNVVALUE <> -1
                ORDER BY POSITION");  
        my $probesRef = fetchAllHashRef();                  
        my ($minI, $maxI, $lowI) = (1, (scalar @$probesRef) - 1, 0);                
        foreach my $i($minI..$maxI){    
             if (${$$probesRef[$i]}{CNVVALUE} != ${$$probesRef[$i - 1]}{CNVVALUE}){                
                commitExternalCNV($lowI, $i - 1, $probesRef, $chrom);
                $lowI = $i;
            }               
        }
        commitExternalCNV($lowI, $maxI, $probesRef, $chrom);
    }
    status("\n");
}

sub commitExternalCNV{
    my ($j, $k, $probesRef, $chrom) = @_;
    my $copyNumber = ${$$probesRef[$j]}{CNVVALUE};
    $copyNumber == 2 and return;
    my ($minPos, $maxPos) = (${$$probesRef[$j]}{POSITION}, ${$$probesRef[$k]}{POSITION});
    ($minPos and $maxPos) or return;    
    my ($n, $testMean, $testStdDev, $refMean, $refStdDev, $normMean, $Z) = getCNVSetStats($probesRef, $j, $k);     
    if ($n){    
        $nmaID++;      
        my $size = $maxPos - $minPos + 1;  
        my ($normRatio, $normZyg) = ($normMean, 0);
        ($isSNPs and !$isHemizygous) and $normZyg = getCrossZygosity($chrom, $minPos, $maxPos); 
        print $nmaFileH join(",", $nmaID, 'External', $nProbes, -1,
                                  $chrom, $minPos, $maxPos, $size, $n, 
                                  $testMean, $testStdDev, $refMean, $refStdDev, 
                                  $normRatio, $normZyg, $Z,
                                  $copyNumber)."\n";    
    }   
}

#sub runTTest{
#    my ($n, $testMean, $refMean, $testSD2, $refSD2) = @_;
#    my $tNum = abs($testMean - $refMean);
#    my $tDenom = sqrt(($testSD2/$n) + ($refSD2/$n));
#    $tDenom or return 0;
#    my $tValue = $tNum / $tDenom;
#    #adjust this if tTable is made 2 sig digits
#    my $tScalar = 10 ** (length(int($tValue)) - 1);
#    $tValue = int($tValue / $tScalar) * $tScalar;    
#    my $dof = $n - 1 ; 
#    #change this is tTable is expanded to larger N    
#    $dof > 200 and $dof = 200;   
#    my $pValue;
#    if(defined $tTable{$dof}{$tValue}){
#        $pValue = $tTable{$dof}{$tValue} 
#    } else {
#        $pValue = 0 #$tValue was very large...
#    }
#    return $pValue;
#}

#sub initializeZTable{
#    my $zTableFile = "$param{vampPath}/bin/zTable.csv";
#    open my $zTableFileH, "<", $zTableFile;
#    my $line = <$zTableFileH>;
#    chomp $line;
#    $line =~ s/\r//g;
#    my ($emptyCell, @lsd) = split(",", $line);
#    my $maxI = (scalar @lsd) - 1;
#    while (<$zTableFileH>){
#        chomp $_;
#        $_ =~ s/\r//g;
#        my ($msd, @halfAreas) = split(",", $_);
#        foreach my $i (0..$maxI){ 
#            my $Z = $msd + $lsd[$i];  
#            my $p = 1 - (2 * $halfAreas[$i]);         
#            $zTable{$Z} = $p; 
#        }
#    }
#    close $zTableFileH;
#}

1;



