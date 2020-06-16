#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
defined $param{geneType} or $param{geneType} = "Unq"; 
defined $param{minND} or $param{minND} = 0.1;  #minimum normalized density value for allowed genes
defined $param{minIntronSize} or $param{minIntronSize} = 500; 
defined $param{minExonSize} or $param{minExonSize} = 50; 
defined $command{compareIRF} or $command{compareIRF} = ['singleThread', '5:00:00', 2000, 0];
defined $command{compareERF} or $command{compareERF} = ['singleThread', '5:00:00', 2000, 0];

my $featureCov;
    
sub compareIRF {
    my ($sample1, $sample2) = @_;
    compareFRF($sample1, $sample2, 'intron', 'IRF', $param{minIntronSize});
}
sub compareERF {
    my ($sample1, $sample2) = @_;
    compareFRF($sample1, $sample2, 'exon', 'ERF', $param{minExonSize});
}
sub compareFRF {
    my ($sample1, $sample2, $featureType, $tableRoot, $minSize) = @_;
    checkFRFDirs($featureType, $tableRoot);
    my ($crosstabTable, $fieldNames) = crosstabFRFs($sample1, $sample2, $featureType, $tableRoot, $minSize);
    printFRFCrosstab($sample1, $sample2, $featureType, $tableRoot, $crosstabTable, $fieldNames);   
    plotFRFCrosstab($crosstabTable, $sample1, $sample2, $tableRoot, $featureType);
}
sub crosstabFRFs {
    my ($sample1, $sample2, $featureType, $tableRoot, $minSize) = @_;
    
    status("performing $tableRoot comparison\n"); 
    
    status("  getting $sample1 $sample2 crosstab\n");  
    my $featureCoverage = $featureType."Coverage";
    $featureCov = $featureType."Cov"; #prevent too-long names
    my $groupBys = [ ['name2'], ['chromosome'], ['start_'], ['end_'] ];
    my $sums = [ ['geneNormalizedDensity', 'geneND'], ['geneCoverage', 'geneCov'], [$featureCoverage, $featureCov], [$tableRoot] ];
    my @frfs = getFRFArray($tableRoot, $minSize, $sample1, $sample2);
    my $where = "  $sample1\_geneND > $param{minND} AND $sample2\_geneND > $param{minND}";
    my ($frfCrosstabSql, $groupByNames, $sampleSumNames) = getSampleMultiCrosstabSql($groupBys, $sums, \@frfs, $where);
    my $ctTmpTable = "ct$tableRoot"."tmp$sample1$sample2";
    updateCrosstabTable($ctTmpTable, $frfCrosstabSql);

    status("  correcting for different sample backgrounds\n"); 
    my $sampleFRF1 = "$sample1\_$tableRoot";
    my $sampleFRF2 = "$sample2\_$tableRoot";   
    runSQL("SELECT avg(log(10,nullif($sampleFRF1,0))) - 
                   avg(log(10,nullif($sampleFRF2,0))) corr
            FROM $ctTmpTable",\my($sample1Corr));
    fetchRow();
    runSQL("UPDATE $ctTmpTable SET $sampleFRF1 = power(10,log(10,nullif($sampleFRF1,0)) - $sample1Corr)");
    
    status("  calculating sample comparisons\n"); 
    my $ratio1Name = "$sample1\_$sample2";
    my $ratio2Name = "$sample2\_$sample1";
    my $fractionName = "fraction_$sample1"; 
    my $statsSql = "SELECT ct.*,
                           $sampleFRF1/nullif($sampleFRF2,0) $ratio1Name,
                           $sampleFRF2/nullif($sampleFRF1,0) $ratio2Name,
                           $sampleFRF2/nullif($sampleFRF1 + $sampleFRF2,0) $fractionName
                    FROM ($ctTmpTable) ct";  

    my $fieldNames = join(",", $groupByNames, $sampleSumNames, $ratio1Name, $ratio2Name, $fractionName);
    my $crosstabTable = "$tableRoot\_$sample1\_$sample2\_$param{geneType}_ct"; 
    updateCrosstabTable($crosstabTable, $statsSql); 
    dropTable($ctTmpTable);
    return ($crosstabTable, $fieldNames);   
}
sub getFRFArray {
    my ($tableRoot, $minSize, @samples) = @_;
    my @frfs;
    foreach my $sample (@samples) {
        $sample or next;
        my $frfTable = getTableName($tableRoot, "$sample\_$param{geneType}");
        my $frfSql = "SELECT * FROM $frfTable WHERE end_ - start_ > $minSize";
        push @frfs, [$sample, $frfSql];
    }
    return @frfs;
}
sub printFRFCrosstab {  #print comparison crosstab to csv file mainly for non-db-friendly users
    my ($sample1, $sample2, $featureType, $tableRoot, $crosstabTable, $fieldNames) = @_;
    status("calculating p-values and printing comparison results\n");   
    my %fieldNames;
    my @fieldNames = split(",", $fieldNames);
    foreach my $i(0..(scalar(@fieldNames)-1)){ $fieldNames{$fieldNames[$i]} = $i }
    runSQL("SELECT $fieldNames FROM $crosstabTable ORDER BY $fieldNames");  
    my $featureFile = "$param{inputPath}/$tableRoot/tables/by_$featureType/comparison/$crosstabTable.csv";         
    open my $featureFileH, ">", $featureFile or die "could not open $featureFile for writing: $!\n"; 
    print $featureFileH "$fieldNames,p-value\n";
    my $geneFile = $featureFile;
    $geneFile =~ s/by_$featureType/by_gene/;
    open my $geneFileH, ">", $geneFile or die "could not open $geneFile for writing: $!\n"; 
    my @groupByFields = ("name2","$sample1\_geneND", "$sample2\_geneND");
    my $groupByFields = join(",", @groupByFields);    
    print $geneFileH "$groupByFields,smallest p-value,$sample1\_$sample2 ...\n";    
    my ($prevGene, $minP, @ratios) = ("");
    while (my @vals = fetchRowArray()){
        my $pValue = getFRFpValue($sample1, $sample2, \@vals, \%fieldNames, $tableRoot);
        fixOracleNulls(\@vals);        
        print $featureFileH join(",", @vals, $pValue)."\n";
        my @groupByVals;
        foreach my $groupByField(@groupByFields){ push @groupByVals, $vals[$fieldNames{$groupByField}] } 
        my $groupByVals = join(",", @groupByVals); 
        if($prevGene and $prevGene ne $groupByVals){
            print $geneFileH join(",", $groupByVals, $minP, @ratios)."\n";
            $minP = undef;
            @ratios = ();
        }
        (defined $minP and $minP <= $pValue) or $minP = $pValue;        
        my $ratio = $vals[$fieldNames{"$sample1\_$sample2"}];
        push @ratios, $ratio;
        $prevGene = $groupByVals;
    }      
    close $featureFileH;  
    close $geneFileH;        
    return $featureFile;                                       
}
sub getFRFpValue {
    my ($sample1, $sample2, $vals, $fieldNames, $tableRoot) = @_;
    my $k1 = $$vals[$$fieldNames{"$sample1\_$featureCov"}];
    my $n1 = $$vals[$$fieldNames{"$sample1\_geneCov"}];
    my $k2 = $$vals[$$fieldNames{"$sample2\_$featureCov"}];
    my $n2 = $$vals[$$fieldNames{"$sample2\_geneCov"}];
    $k1 or $k1 = 0;
    $k2 or $k2 = 0;      
    $n1 or $n1 = 0;
    $n2 or $n2 = 0;       
    if($tableRoot eq 'IRF'){ #since intron hits not included in isMrna geneCov
        $n1 += $k1;
        $n2 += $k2;
    }
    my $pValue = -1;
    ($n1 and $n2) and $pValue = qx/$param{etestPath}etest 1 $k1 $n1 $k2 $n2/;
    chomp $pValue;
    $pValue == -1 and $pValue = "na";
    return $pValue;
}
sub plotFRFCrosstab {
    my ($crosstabTable, $sample1, $sample2, $tableRoot, $featureType) = @_;

    my $include = [];
    my $ratioMax = 100;
    my @ratioAxis = (1/$ratioMax, $ratioMax, 1);
    my $sample1FRF = "$sample1\_$tableRoot";
    my $sample2FRF = "$sample2\_$tableRoot";
    my $ratio1Name = "$sample1\_$sample2";
    my $ratio2Name = "$sample2\_$sample1";    
    my $ratio1Label = "$sample1 / $sample2";
    my $ratio2Label = "$sample2 / $sample1";
    my $ratioBinSize = 1/10;
    my @frfAxis = (1/10000, 10, 1);    
    my $frfLabel = $tableRoot;
    my $frfBinSize = 1/10;
    my $foldInduction = 2;
    my $ratioLines = [1/$foldInduction,1,$foldInduction];
    my $fractionLines = [1/(1+$foldInduction), 1/2, $foldInduction/(1+$foldInduction)];
                 
    status("creating correlation plots\n");
    plotCorrelation($crosstabTable, undef, $include, "$tableRoot/plots/correlation/",
                    $sample1FRF, undef, "$sample1 $frfLabel", @frfAxis,
                    $sample2FRF, undef, "$sample2 $frfLabel", @frfAxis,
                    1, [], []);

    status("creating histogram plots\n");
    plotHistogram($crosstabTable, undef, "$tableRoot/plots/histogram/comparison/",
                  $ratio1Name, undef, $ratio1Label, @ratioAxis, $ratioBinSize,
                  $ratioLines);
    plotHistogram($crosstabTable, undef, "$tableRoot/plots/histogram/comparison/",
                  $ratio2Name, undef, $ratio2Label, @ratioAxis, $ratioBinSize,
                  $ratioLines);    
    plotPairedHistograms($crosstabTable, undef, $sample1FRF, undef, 
                         $crosstabTable, undef, $sample2FRF, undef, 
                         $tableRoot, $frfLabel, @frfAxis, $frfBinSize,
                         "$tableRoot/plots/histogram/by_sample/", []);        
}

#sub compareFRF {
#    my ($sample1, $sample2, $featureType, $tableRoot, $minSize) = @_;
#    checkFRFDirs($featureType, $tableRoot);
#    my ($crosstabTable, $ratio1Name, $ratio2Name, $fractionName, $groupByName) = crosstabFRFs($sample1, $sample2, $featureType, $tableRoot, $minSize);
#    plotFRFCrosstab($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName, $tableRoot, $featureType);
#    printFRFCrosstab($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName, $tableRoot, $featureType);
#    #printFRFCrosstabFiltered($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName, $tableRoot, $featureType);
#}
#sub crosstabFRFs {
#    my ($sample1, $sample2, $featureType, $tableRoot, $minSize) = @_;
#    my @frfs = getFRFArray($tableRoot, $minSize, $sample1, $sample2);
#    my $groupBySql = "name2 || ',' || chromosome || ',' || start_ || ',' || end_";
#    my $groupByName = $featureType."ID";
#    
#    status("getting $sample1 $sample2 $tableRoot crosstab\n");  
#    my ($frfCrosstabSql) = getSampleCrosstabSql($groupBySql, $tableRoot, $groupByName, undef, @frfs);    
#
#    status("getting $sample1 $sample2 geneNormalizedDensity crosstab\n");  
#    my ($ndCrosstabSql) = getSampleCrosstabSql($groupBySql, 'geneNormalizedDensity', $groupByName, undef, @frfs);   
#
#    status("filtering $tableRoot for geneNormalizedDensity > $param{minND}\n");  
#    my $minNDSql = "SELECT * FROM ($ndCrosstabSql) WHERE $sample1 > $param{minND} AND $sample2 > $param{minND}";
#    my $filterSql = "SELECT frf.*
#                     FROM ($frfCrosstabSql) frf, ($minNDSql) nd
#                     WHERE frf.$groupByName = nd.$groupByName";

#    status("correcting for different sample backgrounds\n"); 
#    my $ctTmpTable = "ct$tableRoot"."tmp$sample1$sample2";
#    updateCrosstabTable($ctTmpTable, $filterSql);
#    runSQL("SELECT avg(log(10,nullif($sample1,0))) - 
#                   avg(log(10,nullif($sample2,0))) corr
#            FROM $ctTmpTable",\my($sample1Corr));
#    fetchRow();
#    runSQL("UPDATE $ctTmpTable SET $sample1 = power(10,log(10,nullif($sample1,0)) - $sample1Corr)");
#
#    status("calculating sample ratios\n"); 
#    my ($ratioSql, $ratio1Name, $ratio2Name, $fractionName) = getTwoSampleRatioSql($ctTmpTable, $groupByName, $sample1, $sample2);  
#
#    status("ranking by sample ratio\n"); 
#    my $rankSql = getCrosstabRankSql($ratioSql, $groupByName, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName);    
#    my $crosstabTable = "$tableRoot\_$sample1\_$sample2\_$param{geneType}_ct";           
#    updateCrosstabTable($crosstabTable, $rankSql); 
#    dropTable($ctTmpTable);
#    return ($crosstabTable, $ratio1Name, $ratio2Name, $fractionName, $groupByName);
#}
#sub printFRFCrosstab {  #print comparison crosstab to csv file mainly for non-db-friendly users
#    my ($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName, $tableRoot, $featureType) = @_;
#    my $sql = "SELECT RANK,$groupByName,$sample1,$sample2,$ratio1Name,$ratio2Name,$fractionName
#               FROM $crosstabTable
#               ORDER BY RANK";
#    my $crosstabFile = printPlotData($sql, "$tableRoot/tables/by_$featureType/comparison/$crosstabTable", 
#                         'RANK','GENE','CHROMOSOME','START','END',$sample1,$sample2,$ratio1Name,$ratio2Name,$fractionName);
#
#    open my $crosstabFileH, "<", $crosstabFile or die "could not open $crosstabFile: $!\n";      
#    my (%mins, %maxs);
#    while (my $line = <$crosstabFileH>){
#        chomp $line;
#        my ($rank,$name2,$chrom,$start,$end,$sample1,$sample2,$ratio1,$ratio2,$fraction) = split(",", $line);
#        (defined $mins{$name2} and $mins{$name2} <= $fraction) or $mins{$name2} = $fraction;
#        (defined $maxs{$name2} and $maxs{$name2} >= $fraction) or $maxs{$name2} = $fraction;
#    }            
#    close $crosstabFileH;   
#
#    $crosstabFile =~ s/by_$featureType/by_gene/;
#    open my $outH, ">", $crosstabFile or die "could not open $crosstabFile: $!\n"; 
#    print $outH "GENE,min_$fractionName,max_$fractionName\n";
#    foreach my $name2(sort {$a cmp $b} keys %mins){
#        print $outH "$name2,$mins{$name2},$maxs{$name2}\n";
#    }
#    close $outH;
#                                           
#}
#sub printFRFCrosstabFiltered {
#    my ($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName, $tableRoot, $featureType) = @_;
#    my $filteredCrosstabSql = getFRFFilteredCrosstabSql($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName);
#    printPlotData($filteredCrosstabSql, "$tableRoot/tables/by_$featureType/comparison/$crosstabTable.filtered", 
#                         'RANK','GENE','CHROMOSOME','START','END',,$sample1,$sample2,$ratio1Name,$ratio2Name,$fractionName);
#}
#sub getFRFFilteredCrosstabSql {
#    my($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName) = @_;
#    return "SELECT sum(1) over (order by $ratio1Name rows unbounded preceding) RANK, 
#                   $groupByName,$sample1,$sample2,$ratio1Name,$ratio2Name,$fractionName
#            FROM $crosstabTable 
#            WHERE ($sample1 >= $param{minIRF} OR $sample2 >= $param{minIRF})
#            ORDER BY $ratio1Name";
#}

1;

