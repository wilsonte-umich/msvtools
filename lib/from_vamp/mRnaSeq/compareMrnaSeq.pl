#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));
use Math::Trig;

requireFolders_("$param{vampPath}/bin/crosstab");
requireFolders_("$param{vampPath}/bin/plot");
    
#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{xx} or $command{xx} = ['multiThread', '24:00:00', 5000, 0];
defined $command{compareMrnaSeq} or $command{compareMrnaSeq} = ['singleThread', '1:00:00', 2000, 0];
defined $command{correlateSampleRatios} or $command{correlateSampleRatios} = ['singleThread', '1:00:00', 2000, 0];

my $gmapTrendlineTable = "GMAP_TRENDLINES";

#TODO:
#could and probably should add propagation of coverage counts 
#and calculation of eTest p-values for gene comparisons relative to genome at large
#in a fashion similar to IRF and ERF
#not made a priority since p-values will probably tend to be very favorable
#given the high counts involved - this does NOT mean they would be reproducible
#between same-sample comparisons - error could be in what was counted (the library), not the counting

sub compareMrnaSeq {
    my ($sample1, $sample2, $sample3, $sample4) = @_;
    #$sample1 and $sample2 are the samples to be correlated
    #$sample3 and $sample4 are optional analogous references samples for scaling the results
    #e.g. $sample1=treated 6h, $sample2=untreated 6h, $sample3=treated 0h, $sample4=untreated 0h
    #results are reported as ($sample1/$sample2) / ($sample3/$sample4) if 3 and 4 are provided, 
    #else as $sample1/$sample2 if only 1 and 2 are provided 
    
    checkGMapDirs();
    my ($crosstabTable, $ratio1Name, $ratio2Name, $fractionName) = crosstabGMaps($sample1, $sample2, $sample3, $sample4);
    plotGMapCrosstab($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName);
    printGMapCrosstab($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName);
    printGMapCrosstabFiltered($sample1, $sample2, $crosstabTable);
    printGMapCrosstabOutliers($sample1, $sample2, $crosstabTable);
}

sub checkGMapDirs {
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
                    tables/all
                    tables/correlateRatios
                    tables/filtered
                    tables/outliers);
    foreach my $dir(@dirs) { 
        my $path = "$param{inputPath}/$dir";
        -d $path or mkdir $path;
    }
}

sub crosstabGMaps {
    my ($sample1, $sample2, $sample3, $sample4) = @_;
    my ($crosstabTable, $gMapTable1, $gMapTable2) = getGMapCrosstabTable($sample1, $sample2);
    my ($qCorr, $b) = getGMapTrendline($gMapTable1, $gMapTable2);

    status("getting $sample1 $sample2 crosstab\n");  
    my @gmaps = getGMapArray($sample1, $sample2);
    my ($crosstabSql, $groupByName) = getSampleCrosstabSql("name2", "normalizedDensity", undef, undef, @gmaps);
    my $ctTmpTable = "ctgmaptmp$sample1$sample2";
    updateCrosstabTable($ctTmpTable, $crosstabSql);

    status("transforming to log space\n");  
    my $logSql = "SELECT $groupByName, log(10,nullif($sample1,0)) $sample1, 
                                       log(10,nullif($sample2,0)) $sample2
                  FROM ($ctTmpTable)
                  WHERE $sample1 > 0 
                     OR $sample2 > 0";
    
    status("normalizing so that trendline is unity\n");
    my $xySql = "SELECT $groupByName, $sample1 x, $sample2 - $b y FROM ($logSql)";
    my $polarSql =  "SELECT $groupByName, sqrt(power(x,2) + power(y,2)) r, atan2(y, x) - $qCorr q FROM ($xySql)";
    $xySql = "SELECT $groupByName, r*cos(q) $sample1, r*sin(q) $sample2 FROM ($polarSql)";
    
    status("transforming back to linear space\n"); 
    my $linearSql = "SELECT $groupByName, power(10, $sample1) $sample1, 
                                          power(10, $sample2) $sample2
                     FROM ($xySql)";  
    
    #an initial value of zero in _either_ sample will give null transformations for _both_ samples
    #to preserve that gene's information, recover original values for all null transformations
    status("recovering original values for null transformations\n"); 
    my $recoverSql = "SELECT tr.$groupByName, nvl(tr.$sample1, ct.$sample1) $sample1, 
                                              nvl(tr.$sample2, ct.$sample2) $sample2
                      FROM ($linearSql) tr, ($ctTmpTable) ct
                      WHERE tr.$groupByName = ct.$groupByName";  
                    
    status("calculating sample ratios\n"); 
    my ($ratioSql, $ratio1Name, $ratio2Name, $fractionName) = getTwoSampleRatioSql($recoverSql, $groupByName, $sample1, $sample2);  

    if($sample3) {  #assumes compare previously performed on $sample3, $sample4 alone
        status("applying $sample3 $sample4 ratio normalization\n"); 
        my ($crosstabTableRef) = getGMapCrosstabTable($sample3, $sample4);
        my $ratio3Name = "$sample3\_$sample4";
        my $ratio4Name = "$sample4\_$sample3";
        $ratioSql = "SELECT t.$groupByName, t.$sample1, t.$sample2, 
                              t.$ratio1Name / nullif(r.$ratio3Name,0) $ratio1Name,
                              t.$ratio2Name / nullif(r.$ratio4Name,0) $ratio2Name,
                              t.$fractionName
                       FROM ($ratioSql) t, $crosstabTableRef r
                       WHERE t.$groupByName = r.$groupByName";
    }

    status("ranking by sample ratio\n"); 
    my $rankSql = getCrosstabRankSql($ratioSql, $groupByName, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName);    

    status("associating gene information\n"); 
    my $refGeneGenesTable = getRefGeneTableName();
    my $geneSQL = "SELECT g.chromosome, g.corrstart_, g.end_, g.strand, r.*
                   FROM ($rankSql) r, $refGeneGenesTable g
                   WHERE r.name2 = g.name2
                   ORDER BY r.rank";
    updateCrosstabTable($crosstabTable, $geneSQL);
    dropTable($ctTmpTable);
    return ($crosstabTable, $ratio1Name, $ratio2Name, $fractionName);
}
sub getGMapCrosstabTable {
    my ($sample1, $sample2) = @_;
    my $gMapTable1 = getTableName('GMap', "$sample1\_$param{geneType}");  
    my $gMapTable2 = getTableName('GMap', "$sample2\_$param{geneType}");  
    my $crosstabTable = "gmap_$sample1\_$sample2\_$param{geneType}_ct";
    return ($crosstabTable, $gMapTable1, $gMapTable2);
}
sub getGMapTrendline {
    my ($gMapTable1, $gMapTable2) = @_;
    status("retrieving correlation trendline\n"); 
    runSQL("SELECT x1, y1, x2, y2
            FROM $gmapTrendlineTable 
            WHERE gmapTable1 = '$gMapTable1'
              AND gmapTable2 = '$gMapTable2'", \my($x1, $y1, $x2, $y2));
    fetchRow();
    defined $x1 or die "  could not find entry for $gMapTable1 $gMapTable2 in table $gmapTrendlineTable\n".
                       "  use Graphs.pl to establish correlation trend line and Save to the database\n";
    my $m = ($y2-$y1)/($x2-$x1);
    my $b = $y2 - $m*$x2;
    my $q = atan($m);
    my $qCorr = $q - 3.14159265/4;
    status("  m\t$m\n  b\t$b\n  q\t$q\n"); 
    return ($qCorr, $b);                    
}
sub getGMapArray {
    my (@samples) = @_;
    my @gmaps;
    foreach my $sample (@samples) {
        $sample or next;
        my $gMapTable = getTableName('GMap', "$sample\_$param{geneType}");
        push @gmaps, [$sample, $gMapTable];
    }
    return @gmaps;
}

sub plotGMapCrosstab {
    my ($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName) = @_;
    
    my $include = ['name2','rank'];
    my $ndMax = 1000; 
    my @ndAxis = (1/$ndMax, $ndMax, 1);
    my $maxNDSql = "greatest($sample1, $sample2)";
    my $maxNDName = "maxSampleND";
    my $maxNDLabel = "Max Normalized Density";
    my $ndBinSize = 1/3;
    
    my $ratioMax = 100;
    my @ratioAxis = (1/$ratioMax, $ratioMax, 1);
    my $ratio1Label = "$sample1 / $sample2";
    my $ratio2Label = "$sample2 / $sample1";
    my $ratioBinSize = 1/10;
    
    my @fractionAxis = (0, 1, undef);    
    my $fractionLabel = "Fraction $sample1";
    my $fractionBinSize = 0.05;
    
    my $foldInduction = 2;
    my $ratioLines = [1/$foldInduction,1,$foldInduction];
    my $fractionLines = [1/(1+$foldInduction), 1/2, $foldInduction/(1+$foldInduction)];

    status("creating correlation plots\n");
    plotCorrelation($crosstabTable, undef, $include, "plots/correlation/twoSample/",
                    $sample1, undef, "$sample1 Normalized Density", @ndAxis,
                    $sample2, undef, "$sample2 Normalized Density", @ndAxis,
                    1, [], []);
    plotCorrelation($crosstabTable, undef, $include, "plots/correlation/maxSampleND/ratio/",
                    $maxNDSql, $maxNDName, $maxNDLabel, @ndAxis,
                    $ratio1Name, undef, $ratio1Label, @ratioAxis,
                    undef, $ratioLines, []);
    plotCorrelation($crosstabTable, undef, $include, "plots/correlation/maxSampleND/ratio/",
                    $maxNDSql, $maxNDName, $maxNDLabel, @ndAxis,
                    $ratio2Name, undef, $ratio2Label, @ratioAxis,
                    undef, $ratioLines, []);        
    plotCorrelation($crosstabTable, undef, $include, "plots/correlation/maxSampleND/fraction/",
                    $maxNDSql, $maxNDName, $maxNDLabel, @ndAxis,
                    $fractionName, undef, $fractionLabel, @fractionAxis,
                    undef, $fractionLines, []);

    status("creating histogram plots\n");
    plotHistogram($crosstabTable, undef, "plots/histogram/ratio/",
                  $ratio1Name, undef, $ratio1Label, @ratioAxis, $ratioBinSize,
                  $ratioLines);
    plotHistogram($crosstabTable, undef, "plots/histogram/ratio/",
                  $ratio2Name, undef, $ratio2Label, @ratioAxis, $ratioBinSize,
                  $ratioLines);
    plotHistogram($crosstabTable, undef, "plots/histogram/fraction/",
                  $fractionName, undef, $fractionLabel, @fractionAxis, $fractionBinSize,
                  $fractionLines);    
    plotPairedHistograms($crosstabTable, undef, $sample1, undef, 
                         $crosstabTable, undef, $sample2, undef, 
                         'normDens', 'Normalized Density', @ndAxis, $ndBinSize,
                         "plots/histogram/normalizedDensity/", [1]);               
}
         
sub printGMapCrosstab {  #print comparison crosstab to csv file mainly for non-db-friendly users
    my ($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName) = @_;
    my $sql = "SELECT CHROMOSOME,CORRSTART_,END_,STRAND,RANK,NAME2,
                      $sample1,$sample2,$ratio1Name,$ratio2Name,$fractionName
               FROM $crosstabTable
               ORDER BY RANK";
    return printPlotData($sql, "tables/all/$crosstabTable", 
                         'CHROMOSOME','CORRSTART_','END_','STRAND','RANK','NAME2',$sample1,$sample2,$ratio1Name,$ratio2Name,$fractionName);
}
sub printGMapCrosstabFiltered {
    my ($sample1, $sample2, $crosstabTable) = @_;
    my $filteredCrosstabSql = getFilteredCrosstabSql($sample1, $sample2, $crosstabTable);
    return printPlotData($filteredCrosstabSql, "tables/filtered/$crosstabTable.filtered", 
                         'RANK','GENE','CHROMOSOME','START','END','STRAND',$sample1,$sample2,"$sample1 \/ $sample2","$sample2 \/ $sample1");
}
sub printGMapCrosstabOutliers {
    my($sample1, $sample2, $crosstabTable) = @_;
    my $filteredCrosstabSql = getFilteredCrosstabSql($sample1, $sample2, $crosstabTable);
    my $minRatio = 2;
    my $outlierSql = "SELECT *
                      FROM ($filteredCrosstabSql) 
                      WHERE ($sample1\_$sample2 <= 1/$minRatio OR $sample2\_$sample1 <= 1/$minRatio)
                      ORDER BY $sample1\_$sample2";
    return printPlotData($outlierSql, "tables/outliers/$crosstabTable.outlier", 'RANK','GENE','CHROMOSOME','START','END','STRAND',
                                                                 $sample1,$sample2,"$sample1 \/ $sample2","$sample2 \/ $sample1");             
}
sub getFilteredCrosstabSql {
    my($sample1, $sample2, $crosstabTable, $minND, $minSize) = @_;
    defined $minND or $minND = 0.3;
    defined $minSize or $minSize = 500;
    return "SELECT sum(1) over (order by $sample1\_$sample2 rows unbounded preceding) RANK, 
                   name2,chromosome,corrstart_,end_,strand,$sample1,$sample2,$sample1\_$sample2,$sample2\_$sample1
            FROM $crosstabTable 
            WHERE ($sample1 >= $minND OR $sample2 >= $minND)
              AND end_ - corrstart_ > $minSize
            ORDER BY $sample1\_$sample2";
}

sub correlateSampleRatios {
    my ($sample1, $sample2, $sample3, $sample4, $minND, $minSize) = @_;
    
    checkCorrSampleRatiosDirs();
    
    #join two prior crosstab results
    my ($crosstabTable12) = getGMapCrosstabTable($sample1, $sample2);
    my ($crosstabTable34) = getGMapCrosstabTable($sample3, $sample4);
    my $filteredCrosstabSql12 = getFilteredCrosstabSql($sample1, $sample2, $crosstabTable12, $minND, $minSize);
    my $filteredCrosstabSql34 = getFilteredCrosstabSql($sample3, $sample4, $crosstabTable34, $minND, $minSize); 
    my ($sample1Name, $sample2Name, $sample3Name, $sample4Name) = ($sample1, $sample2, $sample3, $sample4);
    $sample1Name eq $sample3Name and $sample3Name = "$sample3\_";
    $sample2Name eq $sample4Name and $sample4Name = "$sample4\_"; 
    my $joinSql = "SELECT nvl(ct12.name2,ct34.name2) name2,
                          ct12.$sample1 $sample1Name, ct12.$sample2 $sample2Name, ct34.$sample3 $sample3Name, ct34.$sample4 $sample4Name,
                          ct12.$sample1\_$sample2, ct12.$sample2\_$sample1, 
                          ct34.$sample3\_$sample4, ct34.$sample4\_$sample3
                   FROM ($filteredCrosstabSql12) ct12 FULL OUTER JOIN ($filteredCrosstabSql34) ct34
                     ON (ct12.name2 = ct34.name2)";
                     
    #add the gene information   
    my $refGeneGenesTable = getRefGeneTableName();
    my $geneSQL = "SELECT g.chromosome, g.corrstart_, g.end_, g.strand, j.*
                   FROM ($joinSql) j, $refGeneGenesTable g
                   WHERE j.name2 = g.name2";         
                 
    #export the merged table
    printPlotData($geneSQL, "tables/$sample1\_$sample2.$sample3\_$sample4.merged", 
                            "Chromosome", "Start", "End", "Strand",
                            'Gene',$sample1Name, $sample2Name, $sample3Name, $sample4Name,
                            "$sample1 \/ $sample2","$sample2 \/ $sample1","$sample3 \/ $sample4","$sample4 \/ $sample3");   

    #plot the resulting ratio correlation
    #assumes user is looking for induced stability and or synthesis is sample1/3
    #so put it in the denominator as the larger value                
    my $include = ['name2'];
    my $ratioMax = 1000;
    my @ratioAxis = (1/$ratioMax, $ratioMax, 1);
    my $foldInduction = 2;
    my $ratioLines = [1/$foldInduction,1,$foldInduction];
    plotCorrelation($joinSql, undef, $include, "plots/",
                    "$sample2\_$sample1", undef, undef, @ratioAxis,
                    "$sample4\_$sample3", undef, undef, @ratioAxis,
                    undef, $ratioLines, $ratioLines);
}
sub checkCorrSampleRatiosDirs {
    my @dirs = qw(  plots
                    tables );
    foreach my $dir(@dirs) { 
        my $path = "$param{inputPath}/$dir";
        -d $path or mkdir $path;
    }
}

1;
    
