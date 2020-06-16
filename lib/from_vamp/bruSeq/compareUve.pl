#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %fieldNames));

#callable parameters and commands in this script
defined $command{compareUVE} or $command{compareUVE} = ['singleThread', '12:00:00', 2000, 0];

my ($ratio1Name, $ratio2Name, $frac1Name);

sub compareUVE {
    my ($nonUVsample1, $uvSample1, $nonUVsample2, $uvSample2) = @_;
    setUveNames($uvSample1, $uvSample2);
    checkUveDirs();
    my $uveTmpTable = getMergedUveList($nonUVsample1, $uvSample1, $nonUVsample2, $uvSample2);
    my $uveTable = runUveCompare($uveTmpTable, $uvSample1, $uvSample2);
    plotCompareUve($uveTable, $uvSample1, $uvSample2);
    printCompareUve($uveTable, $uvSample1, $uvSample2);
}
sub setUveNames {
    my ($uvSample1, $uvSample2) = @_;
    $ratio1Name = "$uvSample1\_$uvSample2";
    $ratio2Name = "$uvSample2\_$uvSample1";
    $frac1Name = "fraction_$uvSample1";
}
sub checkUveDirs {
    my @dirs = qw(  plots
                    plots/correlation
                    plots/histogram
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
sub getMergedUveList { #generate the list of all uve in either sample, overlapping uve are merged into one larger uve
    my ($nonUVsample1, $uvSample1, $nonUVsample2, $uvSample2) = @_;
    status("creating merged UVE list\n"); 
    my $uveTmpTable = "uvetmp$nonUVsample1$nonUVsample2";    
    dropTable($uveTmpTable);
    runSQL("CREATE TABLE $uveTmpTable (chromosome NUMBER, start_ NUMBER, end_ NUMBER, inGene NUMBER)"); 
    my $uveTmpFile = "$uveTmpTable.csv";
    open my $uveTmpFileH, ">", $uveTmpFile or die "could not open $uveTmpFile for writing: $!\n";      
    my $uveTable1 = getTableName("HB", "$uvSample1\_uve");
    my $uveTable2 = getTableName("HB", "$uvSample2\_uve");    
    foreach my $chrom(1..nChrom()){    
        my $uveSql1 = "SELECT start_, end_, inGene FROM $uveTable1 WHERE chromosome = $chrom";
        my $uveSql2 = "SELECT start_, end_, inGene FROM $uveTable2 WHERE chromosome = $chrom";
        my $unionSql =   "SELECT start_, end_, inGene FROM ($uveSql1) UNION ALL ($uveSql2)";
        my $groupBySql = "SELECT start_, end_, inGene FROM ($unionSql) 
                          GROUP BY start_, end_, inGene
                          ORDER BY start_, end_";    
        runSQL($groupBySql, \my($start,$end,$inGene));
        my ($uveStart, $uveEnd, $uveInGene) = (0, 0, 0); 
        while(fetchRow()){
            if ($uveEnd and $start > $uveEnd){ 
                print $uveTmpFileH "$chrom,$uveStart,$uveEnd,$uveInGene\n";
                ($uveStart, $uveEnd, $uveInGene) = (0, 0, 0); 
            }
            $uveStart or $uveStart = $start;
            $uveEnd >= $end or $uveEnd = $end;
            $uveInGene = ($uveInGene or $inGene);
        }
        $uveEnd and print $uveTmpFileH "$chrom,$uveStart,$uveEnd,$uveInGene\n";
    }
    close $uveTmpFileH;
    loadData($uveTmpFile, $uveTmpTable, ",", "CHROMOSOME, START_, END_, INGENE");
    return $uveTmpTable;
}
sub runUveCompare {
    my ($uveTmpTable, $uvSample1, $uvSample2) = @_;
    status("parsing sample hit counts\n"); 
    my ($hitsTable1, $expectedDensity1) = getUveHitInfo($uvSample1);
    my ($hitsTable2, $expectedDensity2) = getUveHitInfo($uvSample2);   
    my $uveTable = getTableName('UVE', "$uvSample1\_$uvSample2");
    dropTable($uveTable);
    runSQL("CREATE TABLE $uveTable (chromosome NUMBER, start_ NUMBER, end_ NUMBER, inGene NUMBER,
                                    $uvSample1 NUMBER(*,5), $uvSample2 NUMBER(*,5),
                                    $ratio1Name NUMBER(*,5), $ratio2Name NUMBER(*,5), 
                                    $frac1Name NUMBER(*,5))");   
    my $uveFile = "$uveTable.csv";                                                       
    open my $uveFileH, ">", $uveFile or die "could not open $uveFile for writing: $!\n";                                                            
    status("Chr: ");
    foreach my $chrom(1..nChrom()){
    #foreach my $chrom(1..1){
        status(" $chrom");
        my $ndSql1 = getUveNDSql($hitsTable1, $uveTmpTable, $expectedDensity1, $chrom);
        my $ndSql2 = getUveNDSql($hitsTable2, $uveTmpTable, $expectedDensity2, $chrom);
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
            print $uveFileH "$chrom,$start,$end,$inGene,$nd1,$nd2,$r1,$r2,$f1\n" 
        }
    }
    status("\n");    
    close $uveFileH;
    my $fieldNames = "CHROMOSOME, START_, END_, INGENE, $uvSample1, $uvSample2, $ratio1Name, $ratio2Name, $frac1Name";
    loadData($uveFile, $uveTable, ",", $fieldNames);   
    dropTable($uveTmpTable);
    return $uveTable;
}
sub getUveHitInfo { #get the information needed to parse the hit densities for uv samples 1 and 2
    my ($uvSample) = @_;
    my $hitsTable = getTableName("Hits", $uvSample);
    my $genomeSize = getGenomeSize();
    my $nHits = getGenomeHitCount($hitsTable);
    my $expectedDensity = $nHits / $genomeSize;   
    return ($hitsTable, $expectedDensity);
}
sub getUveNDSql {
    my ($hitsTable, $uveTmpTable, $expectedDensity, $chrom) = @_;
    my $agg = "count(*)";
    $param{keepDups} and $agg = "sum(h.count_)";
    my $hitsSql = "SELECT * FROM $hitsTable   WHERE chromosome = $chrom";    
    my $uveSql  = "SELECT * FROM $uveTmpTable WHERE chromosome = $chrom";
    my $ndSql =  "SELECT uve.start_, uve.end_, uve.inGene,
                         (nvl($agg,0)/(uve.end_ - uve.start_)) / $expectedDensity ND
                     FROM ($hitsSql) h, ($uveSql) uve
                     WHERE h.position(+) BETWEEN uve.start_ AND uve.end_
                     GROUP BY uve.start_, uve.end_, uve.inGene";   
    return $ndSql;
}
sub plotCompareUve {
    my ($uveTable, $uvSample1, $uvSample2) = @_;                
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
    plotCorrelation($uveTable, undef, $include, "plots/correlation/twoSample/",
                    $uvSample1, undef, "$uvSample1 Normalized Density", @ndAxis,
                    $uvSample2, undef, "$uvSample2 Normalized Density", @ndAxis,
                    1, [], []);
    plotCorrelation($uveTable, undef, $include, "plots/correlation/maxSampleND/ratio/",
                    $maxNDSql, $maxNDName, $maxNDLabel, @ndAxis,
                    $ratio1Name, undef, $ratio1Label, @ratioAxis,
                    undef, $ratioLines, []);
    plotCorrelation($uveTable, undef, $include, "plots/correlation/maxSampleND/ratio/",
                    $maxNDSql, $maxNDName, $maxNDLabel, @ndAxis,
                    $ratio2Name, undef, $ratio2Label, @ratioAxis,
                    undef, $ratioLines, []);        
    plotCorrelation($uveTable, undef, $include, "plots/correlation/maxSampleND/fraction/",
                    $maxNDSql, $maxNDName, $maxNDLabel, @ndAxis,
                    $frac1Name, undef, $fractionLabel, @fractionAxis,
                    undef, $fractionLines, []);

    status("creating histogram plots\n");
    plotHistogram($uveTable, undef, "plots/histogram/ratio/",
                  $ratio1Name, undef, $ratio1Label, @ratioAxis, $ratioBinSize,
                  $ratioLines);
    plotHistogram($uveTable, undef, "plots/histogram/ratio/",
                  $ratio2Name, undef, $ratio2Label, @ratioAxis, $ratioBinSize,
                  $ratioLines);
    plotHistogram($uveTable, undef, "plots/histogram/fraction/",
                  $frac1Name, undef, $fractionLabel, @fractionAxis, $fractionBinSize,
                  $fractionLines);    
    plotPairedHistograms($uveTable, undef, $uvSample1, undef, 
                         $uveTable, undef, $uvSample2, undef, 
                         'normDens', 'Normalized Density', @ndAxis, $ndBinSize,
                         "plots/histogram/normalizedDensity/", [1]);   
    plotPairedHistograms($uveTable, "inGene > 0", $uvSample1, 'inGene', 
                         $uveTable, "inGene = 0", $uvSample1, 'notInGene', 
                         'normDens', 'Normalized Density', @ndAxis, 1/10,
                         "plots/histogram/normalizedDensity/$uvSample1", []);  
    plotPairedHistograms($uveTable, "inGene > 0", $uvSample2, 'inGene', 
                         $uveTable, "inGene = 0", $uvSample2, 'notInGene', 
                         'normDens', 'Normalized Density', @ndAxis, 1/10,
                         "plots/histogram/normalizedDensity/$uvSample2", []);                         
}
sub printCompareUve {  
    my ($uveTable, $uvSample1, $uvSample2) = @_;
    my $sql = "SELECT uve.*, end_ - start_ size_
               FROM $uveTable uve
               ORDER BY $ratio1Name";
    return printPlotData($sql, "tables/all/$uveTable", 
                         "CHROMOSOME", "START_", "END_", "INGENE", $uvSample1, $uvSample2, $ratio1Name, $ratio2Name, $frac1Name, "GENESIZE");
}

1;

