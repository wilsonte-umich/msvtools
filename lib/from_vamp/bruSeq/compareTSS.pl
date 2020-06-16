#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %fieldNames));

#callable parameters and commands in this script
defined $command{compareTSS} or $command{compareTSS} = ['singleThread', '12:00:00', 2000, 0];

my ($ratio1Name, $ratio2Name, $frac1Name);

sub compareTSS {
    my ($nonUVsample1, $uvSample1, $nonUVsample2, $uvSample2) = @_;
    setTssNames($uvSample1, $uvSample2);
    checkTssDirs();
    my $tssTmpTable = getMergedTssList($nonUVsample1, $uvSample1, $nonUVsample2, $uvSample2);
    my $tssTable = runTssCompare($tssTmpTable, $uvSample1, $uvSample2);
    plotCompareTss($tssTable, $uvSample1, $uvSample2);
    printCompareTss($tssTable, $uvSample1, $uvSample2);
}
sub setTssNames {
    my ($uvSample1, $uvSample2) = @_;
    $ratio1Name = "$uvSample1\_$uvSample2";
    $ratio2Name = "$uvSample2\_$uvSample1";
    $frac1Name = "fraction_$uvSample1";
}
sub checkTssDirs {
    my @dirs = qw(  plots
                    plots/correlation
                    plots/histogram
                    plots/correlation/correlateRatios
                    plots/correlation/maxSampleND
                    plots/correlation/twoSample
                    plots/correlation/maxSampleND/fraction
                    plots/correlation/maxSampleND/ratio
                    plots/histogram/fraction
                    plots/histogram/normalizedDensity
                    plots/histogram/ratio
                    tables
                    tables/all);
    foreach my $dir(@dirs) { 
        my $path = "$param{inputPath}/$dir";
        -d $path or mkdir $path;
    }
}
sub getMergedTssList { #generate the list of all tss in either sample, overlapping tss are merged into one larger tss
    my ($nonUVsample1, $uvSample1, $nonUVsample2, $uvSample2) = @_;
    status("creating merged tss list\n"); 
    my $tssTmpTable = "tsstmp$nonUVsample1$nonUVsample2";    
    dropTable($tssTmpTable);
    runSQL("CREATE TABLE $tssTmpTable (chromosome NUMBER, start_ NUMBER, end_ NUMBER, inGene NUMBER)"); 
    my $tssTmpFile = "$tssTmpTable.csv";
    open my $tssTmpFileH, ">", $tssTmpFile or die "could not open $tssTmpFile for writing: $!\n";      
    my $tssTable1 = getTableName("HB", "$uvSample1\_tss");
    my $tssTable2 = getTableName("HB", "$uvSample2\_tss");    
    foreach my $chrom(1..nChrom()){    
        my $tssSql1 = "SELECT start_, end_, inGene FROM $tssTable1 WHERE chromosome = $chrom";
        my $tssSql2 = "SELECT start_, end_, inGene FROM $tssTable2 WHERE chromosome = $chrom";
        my $unionSql =   "SELECT start_, end_, inGene FROM ($tssSql1) UNION ALL ($tssSql2)";
        my $groupBySql = "SELECT start_, end_, inGene FROM ($unionSql) 
                          GROUP BY start_, end_, inGene
                          ORDER BY start_, end_";    
        runSQL($groupBySql, \my($start,$end,$inGene));
        my ($tssStart, $tssEnd, $tssInGene) = (0, 0, 0); 
        while(fetchRow()){
            if ($tssEnd and $start > $tssEnd){ 
                print $tssTmpFileH "$chrom,$tssStart,$tssEnd,$tssInGene\n";
                ($tssStart, $tssEnd, $tssInGene) = (0, 0, 0); 
            }
            $tssStart or $tssStart = $start;
            $tssEnd >= $end or $tssEnd = $end;
            $tssInGene = ($tssInGene or $inGene);
        }
        $tssEnd and print $tssTmpFileH "$chrom,$tssStart,$tssEnd,$tssInGene\n";
    }
    close $tssTmpFileH;
    loadData($tssTmpFile, $tssTmpTable, ",", "CHROMOSOME, START_, END_, INGENE");
    return $tssTmpTable;
}
sub runTssCompare {
    my ($tssTmpTable, $uvSample1, $uvSample2) = @_;
    status("parsing sample hit counts\n"); 
    my ($hitsTable1, $expectedDensity1) = getTssHitInfo($uvSample1);
    my ($hitsTable2, $expectedDensity2) = getTssHitInfo($uvSample2);   
    my $tssTable = getTableName('TSS', "$uvSample1\_$uvSample2");
    dropTable($tssTable);
    runSQL("CREATE TABLE $tssTable (chromosome NUMBER, start_ NUMBER, end_ NUMBER, inGene NUMBER,
                                    $uvSample1 NUMBER(*,5), $uvSample2 NUMBER(*,5),
                                    $ratio1Name NUMBER(*,5), $ratio2Name NUMBER(*,5), 
                                    $frac1Name NUMBER(*,5))");   
    my $tssFile = "$tssTable.csv";                                                       
    open my $tssFileH, ">", $tssFile or die "could not open $tssFile for writing: $!\n";                                                            
    status("Chr: ");
    foreach my $chrom(1..nChrom()){
    #foreach my $chrom(1..1){
        status(" $chrom");
        my $ndSql1 = getTssNDSql($hitsTable1, $tssTmpTable, $expectedDensity1, $chrom);
        my $ndSql2 = getTssNDSql($hitsTable2, $tssTmpTable, $expectedDensity2, $chrom);
        my $joinSql = "SELECT nd1.start_, nd1.end_, nd1.inGene,
                              nd1.ND $uvSample1, nd2.ND $uvSample2,
                              nd1.ND/nullif(nd2.ND,0) $ratio1Name, 
                              nd2.ND/nullif(nd1.ND,0) $ratio2Name, 
                              nd1.ND/(nd1.ND + nd2.ND) $frac1Name
                       FROM ($ndSql1) nd1, ($ndSql2) nd2
                       WHERE nd1.start_ = nd2.start_
                         AND nd1.end_ = nd2.end_"; 
        runSQL($joinSql, \my($start,$end,$inGene,$nd1,$nd2,$r1,$r2,$f1));
        while (fetchRow()){ 
            $r1 or $r1 = "";
            $r2 or $r2 = "";
            print $tssFileH "$chrom,$start,$end,$inGene,$nd1,$nd2,$r1,$r2,$f1\n" 
        }
    }
    status("\n");    
    close $tssFileH;
    my $fieldNames = "CHROMOSOME, START_, END_, INGENE, $uvSample1, $uvSample2, $ratio1Name, $ratio2Name, $frac1Name";
    loadData($tssFile, $tssTable, ",", $fieldNames);   
    dropTable($tssTmpTable);
    return $tssTable;
}
sub getTssHitInfo { #get the information needed to parse the hit densities for uv samples 1 and 2
    my ($uvSample) = @_;
    my $hitsTable = getTableName("Hits", $uvSample);
    my $genomeSize = getGenomeSize();
    my $nHits = getGenomeHitCount($hitsTable);
    my $expectedDensity = $nHits / $genomeSize;   
    return ($hitsTable, $expectedDensity);
}
sub getTssNDSql {
    my ($hitsTable, $tssTmpTable, $expectedDensity, $chrom) = @_;
    my $agg = "count(*)";
    $param{keepDups} and $agg = "sum(h.count_)";
    my $hitsSql = "SELECT * FROM $hitsTable   WHERE chromosome = $chrom";    
    my $tssSql  = "SELECT * FROM $tssTmpTable WHERE chromosome = $chrom";
    my $ndSql =  "SELECT tss.start_, tss.end_, tss.inGene,
                         (nvl($agg,0)/(tss.end_ - tss.start_)) / $expectedDensity ND
                     FROM ($hitsSql) h, ($tssSql) tss
                     WHERE h.position(+) BETWEEN tss.start_ AND tss.end_
                     GROUP BY tss.start_, tss.end_, tss.inGene";   
    return $ndSql;
}
sub plotCompareTss {
    my ($tssTable, $uvSample1, $uvSample2) = @_;                
    my $include = [];
    my $ndMax = 1000; 
    my @ndAxis = (1/$ndMax, $ndMax, 1);
    my $maxNDSql = "greatest($uvSample1, $uvSample2)";
    my $maxNDName = "maxSampleND";
    my $maxNDLabel = "Max Normalized Density";
    my $ndBinSize = 1/3;
    
    my $ratioMax = 100;
    my @ratioAxis = (1/$ratioMax, $ratioMax, 1);

    my $ratio1Label = "$uvSample1 / $uvSample2";
    my $ratio2Label = "$uvSample2 / $uvSample1";
    my $ratioBinSize = 1/10;
    
    my @fractionAxis = (0, 1, undef);    
    my $fractionLabel = "Fraction $uvSample1";
    my $fractionBinSize = 0.05;

    my $foldInduction = 2;
    my $ratioLines = [1/$foldInduction,1,$foldInduction];
    my $fractionLines = [1/(1+$foldInduction), 1/2, $foldInduction/(1+$foldInduction)];      

    status("creating correlation plots\n");
    plotCorrelation($tssTable, undef, $include, "plots/correlation/twoSample/",
                    $uvSample1, undef, "$uvSample1 Normalized Density", @ndAxis,
                    $uvSample2, undef, "$uvSample2 Normalized Density", @ndAxis,
                    1, [], []);
    plotCorrelation($tssTable, undef, $include, "plots/correlation/maxSampleND/ratio/",
                    $maxNDSql, $maxNDName, $maxNDLabel, @ndAxis,
                    $ratio1Name, undef, $ratio1Label, @ratioAxis,
                    undef, $ratioLines, []);
    plotCorrelation($tssTable, undef, $include, "plots/correlation/maxSampleND/ratio/",
                    $maxNDSql, $maxNDName, $maxNDLabel, @ndAxis,
                    $ratio2Name, undef, $ratio2Label, @ratioAxis,
                    undef, $ratioLines, []);        
    plotCorrelation($tssTable, undef, $include, "plots/correlation/maxSampleND/fraction/",
                    $maxNDSql, $maxNDName, $maxNDLabel, @ndAxis,
                    $frac1Name, undef, $fractionLabel, @fractionAxis,
                    undef, $fractionLines, []);

    status("creating histogram plots\n");
    plotHistogram($tssTable, undef, "plots/histogram/ratio/",
                  $ratio1Name, undef, $ratio1Label, @ratioAxis, $ratioBinSize,
                  $ratioLines);
    plotHistogram($tssTable, undef, "plots/histogram/ratio/",
                  $ratio2Name, undef, $ratio2Label, @ratioAxis, $ratioBinSize,
                  $ratioLines);
    plotHistogram($tssTable, undef, "plots/histogram/fraction/",
                  $frac1Name, undef, $fractionLabel, @fractionAxis, $fractionBinSize,
                  $fractionLines);    
    plotPairedHistograms($tssTable, undef, $uvSample1, undef, 
                         $tssTable, undef, $uvSample2, undef, 
                         'normDens', 'Normalized Density', @ndAxis, $ndBinSize,
                         "plots/histogram/normalizedDensity/", [1]);   
    plotPairedHistograms($tssTable, "inGene > 0", $uvSample1, 'inGene', 
                         $tssTable, "inGene = 0", $uvSample1, 'notInGene', 
                         'normDens', 'Normalized Density', @ndAxis, 1/10,
                         "plots/histogram/normalizedDensity/$uvSample1", []);  
    plotPairedHistograms($tssTable, "inGene > 0", $uvSample2, 'inGene', 
                         $tssTable, "inGene = 0", $uvSample2, 'notInGene', 
                         'normDens', 'Normalized Density', @ndAxis, 1/10,
                         "plots/histogram/normalizedDensity/$uvSample2", []);                         
}
sub printCompareTss {  
    my ($tssTable, $uvSample1, $uvSample2) = @_;
    my $sql = "SELECT tss.*, end_ - start_ size_
               FROM $tssTable tss
               ORDER BY $ratio1Name";
    return printPlotData($sql, "tables/all/$tssTable", 
                         "CHROMOSOME", "START_", "END_", "SIZE_". "INGENE", $uvSample1, $uvSample2, $ratio1Name, $ratio2Name, $frac1Name);
}

1;

