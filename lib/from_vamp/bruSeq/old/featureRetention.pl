#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
defined $param{geneType} or $param{geneType} = "Unq"; 
defined $param{minND} or $param{minND} = 0.1;  #minimum normalized density value for allowed genes

defined $param{minIntronSize} or $param{minIntronSize} = 500; 

defined $param{minIRF} or $param{minIRF} = 0.05;  #minimum intron retention fraction for filtered list

defined $command{calculateIRF} or $command{calculateIRF} = ['multiThread', '1:00:00', 1000, 0];
defined $command{calculateERF} or $command{calculateERF} = ['multiThread', '1:00:00', 1000, 0];
defined $command{compareIRF} or $command{compareIRF} = ['singleThread', '1:00:00', 1000, 0];
defined $command{compareERF} or $command{compareERF} = ['singleThread', '1:00:00', 1000, 0];
    
#samples expected to be 6h sample, i.e. isMrna = TRUE during parseGeneHits
sub calculateIRF { #IRF = intron retention fraction, for finding unspliced introns
    my ($sample) = @_; 
    calculateFRF($sample, 'intron', 'IMap', 'IRF');
}
sub calculateERF { #ERF = exon retention fraction, for finding alternatively spliced transcripts
    my ($sample) = @_;
    calculateFRF($sample, 'exon', 'EMap', 'ERF');  
}
sub calculateFRF { #FRF - generic Feature Retention Fraction
    my ($sample, $featureType, $mapType, $tableRoot) = @_; 
    
    print "$sample, $featureType, $mapType, $tableRoot";
    
    my ($frfTable, $featureDensity) = createFRFTable($sample, $featureType, $mapType, $tableRoot);
    checkIRFDirs($featureType, $tableRoot);
    printFRFTables($featureType, $tableRoot, $frfTable, $featureDensity);
}
sub createFRFTable {
    my ($sample, $featureType, $mapType, $tableRoot) = @_; 
    status("calculating $sample $featureType retention fraction\n");  
    my $fmapTable = getTableName($mapType, $sample);
    my $gmapTable = getTableName('GMap', "$sample\_$param{geneType}");  
    my $frfTable = getTableName($tableRoot, "$sample\_$param{geneType}");
    my $featureDensity = $featureType."Density";
    my $joinSql = "SELECT f.name2, f.chromosome, f.start_, f.end_, 
                          f.density $featureDensity, g.density geneDensity, g.normalizedDensity geneNormalizedDensity,
                          f.density/nullif(g.density,0) $tableRoot
                   FROM $fmapTable f, $gmapTable g
                   WHERE f.name2 = g.name2";            
    dropTable($frfTable);      
    runSQL("CREATE TABLE $frfTable AS $joinSql");
    return ($frfTable, $featureDensity);
}
sub checkIRFDirs {
    my ($featureType, $tableRoot) = @_; 
    my $frfPath = "$param{inputPath}/$tableRoot";
    -d $frfPath or mkdir $frfPath;
    my @dirs = ("plots",
                "plots/correlation",
                "plots/histogram",
                "plots/histogram/by_sample",
                "plots/histogram/comparison",
                "tables",   
                "tables/by_$featureType",
                "tables/by_$featureType/by_sample",
                "tables/by_$featureType/comparison",
                "tables/by_gene",
                "tables/by_gene/by_sample", 
                "tables/by_gene/comparison");
    foreach my $dir(@dirs) { 
        my $path = "$frfPath/$dir";
        -d $path or mkdir $path;
    }
}
sub printFRFTables {
    my ($featureType, $tableRoot, $frfTable, $featureDensity) = @_;
    status("printing output tables\n");  
    my $printSql = "SELECT sum(1) over (order by $tableRoot DESC rows unbounded preceding) RANK, 
                           name2, chromosome, start_, end_, end_ - start_ + 1 size_,
                          $featureDensity, geneDensity, geneNormalizedDensity, $tableRoot
                    FROM $frfTable
                    WHERE geneNormalizedDensity > $param{minND}
                    ORDER BY $tableRoot DESC";
    printPlotData($printSql, "$tableRoot/tables/by_$featureType/by_sample/$frfTable", 
                         'RANK','GENE','CHROMOSOME','START','END','SIZE',
                          $featureDensity,'geneDensity','geneNormalizedDensity',$tableRoot);
 
    my $geneTable = "refgene_$param{geneType}_$param{refSeq}";
    $printSql = "SELECT g.name2, g.chromosome, g.corrstart_ start_, g.end_, g.genesize, g.mrnasize,
                        f.geneDensity, f.geneNormalizedDensity,
                        min($featureDensity) min_$featureDensity, max($featureDensity) max_$featureDensity,
                        min($tableRoot) min_$tableRoot, max($tableRoot) max_$tableRoot
                 FROM $frfTable f, $geneTable g
                 WHERE f.name2 = g.name2
                   AND f.geneNormalizedDensity > $param{minND}
                 GROUP BY g.name2, g.chromosome, g.corrstart_, g.end_, g.genesize, g.mrnasize,
                          f.geneDensity, f.geneNormalizedDensity";   
    $printSql = "SELECT sum(1) over (order by p.max_$tableRoot rows unbounded preceding) RANK, p.* 
                 FROM ($printSql) p
                 ORDER BY p.max_$tableRoot";      
    printPlotData($printSql, "$tableRoot/tables/by_gene/by_sample/$frfTable", 
                         'RANK','GENE','CHROMOSOME','START','END','GENESIZE','MRNASIZE',
                         'geneDensity','geneNormalizedDensity',
                         "min_$featureDensity", "max_$featureDensity",
                         "min_$tableRoot", "max_$tableRoot");                                
}






sub compareIRF {
    my ($sample1, $sample2) = @_;
    checkIRFDirs();
    my ($crosstabTable, $ratio1Name, $ratio2Name, $fractionName, $groupByName) = crosstabIRFs($sample1, $sample2);
    plotIRFCrosstab($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName);
    printIRFCrosstab($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName);
    printIRFCrosstabFiltered($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName);
}






sub crosstabIRFs {
    my ($sample1, $sample2) = @_;
    my @irfs = getIRFArray($sample1, $sample2);
    my $groupBySql = "name2 || ',' || chromosome || ',' || start_ || ',' || end_";
    my $groupByName = 'intronID';
    
    status("getting $sample1 $sample2 retentionFraction crosstab\n");  
    my ($irfCrosstabSql) = getSampleCrosstabSql($groupBySql, 'retentionFraction', $groupByName, undef, @irfs);    

    status("getting $sample1 $sample2 normalizedDensity crosstab\n");  
    my ($ndCrosstabSql) = getSampleCrosstabSql($groupBySql, 'exonNormalizedDensity', $groupByName, undef, @irfs);   

    status("filtering retentionFraction for normalizedDensity > $param{minND}\n");  
    my $minNDSql = "SELECT * FROM ($ndCrosstabSql) WHERE $sample1 > $param{minND} AND $sample2 > $param{minND}";
    my $filterSql = "SELECT irf.*
                     FROM ($irfCrosstabSql) irf, ($minNDSql) nd
                     WHERE irf.$groupByName = nd.$groupByName";

    status("correcting for different sample backgrounds\n"); 
    my $ctTmpTable = "ctirftmp$sample1$sample2";
    updateCrosstabTable($ctTmpTable, $filterSql);
    runSQL("SELECT avg(log(10,nullif($sample1,0))) - 
                   avg(log(10,nullif($sample2,0))) corr
            FROM $ctTmpTable",\my($sample1Corr));
    fetchRow();
    runSQL("UPDATE $ctTmpTable SET $sample1 = power(10,log(10,nullif($sample1,0)) - $sample1Corr)");

    status("calculating sample ratios\n"); 
    my ($ratioSql, $ratio1Name, $ratio2Name, $fractionName) = getTwoSampleRatioSql($ctTmpTable, $groupByName, $sample1, $sample2);  

    status("ranking by sample ratio\n"); 
    my $rankSql = getCrosstabRankSql($ratioSql, $groupByName, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName);    

    my $crosstabTable = "IRF_$sample1\_$sample2\_$param{geneType}_ct";               
    updateCrosstabTable($crosstabTable, $rankSql); 
    dropTable($ctTmpTable);
    return ($crosstabTable, $ratio1Name, $ratio2Name, $fractionName, $groupByName);
}
sub getIRFArray {
    my (@samples) = @_;
    my @irfs;
    foreach my $sample (@samples) {
        $sample or next;
        my $irfTable = getTableName('IRF', "$sample\_$param{geneType}");
        my $irfSql = "SELECT * FROM $irfTable WHERE end_ - start_ > $param{minIntronSize}";
        push @irfs, [$sample, $irfSql];
    }
    return @irfs;
}

sub printIRFCrosstab {  #print comparison crosstab to csv file mainly for non-db-friendly users
    my ($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName) = @_;
    my $sql = "SELECT RANK,$groupByName,$sample1,$sample2,$ratio1Name,$ratio2Name,$fractionName
               FROM $crosstabTable
               ORDER BY RANK";
    return printPlotData($sql, "IRF/tables/all/$crosstabTable", 
                         'RANK','GENE','CHROMOSOME','START','END',$sample1,$sample2,$ratio1Name,$ratio2Name,$fractionName);
}
sub printIRFCrosstabFiltered {
    my ($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName) = @_;
    my $filteredCrosstabSql = getIRFFilteredCrosstabSql($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName);
    return printPlotData($filteredCrosstabSql, "IRF/tables/filtered/$crosstabTable.filtered", 
                         'RANK','GENE','CHROMOSOME','START','END',,$sample1,$sample2,$ratio1Name,$ratio2Name,$fractionName);
}
sub getIRFFilteredCrosstabSql {
    my($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName) = @_;
    return "SELECT sum(1) over (order by $ratio1Name rows unbounded preceding) RANK, 
                   $groupByName,$sample1,$sample2,$ratio1Name,$ratio2Name,$fractionName
            FROM $crosstabTable 
            WHERE ($sample1 >= $param{minIRF} OR $sample2 >= $param{minIRF})
            ORDER BY $ratio1Name";
}

sub plotIRFCrosstab {
    my ($crosstabTable, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName, $groupByName) = @_;
    
    my $include = ['RANK'];

    my $ratioMax = 100;
    my @ratioAxis = (1/$ratioMax, $ratioMax, 1);
    my $ratio1Label = "$sample1 / $sample2";
    my $ratio2Label = "$sample2 / $sample1";
    my $ratioBinSize = 1/10;

    my @irfAxis = (1/10000, 10, 1);    
    my $irfLabel = "Intron Retention Fraction";
    my $irfBinSize = 1/10;

    my $foldInduction = 2;
    my $ratioLines = [1/$foldInduction,1,$foldInduction];
    my $fractionLines = [1/(1+$foldInduction), 1/2, $foldInduction/(1+$foldInduction)];

    status("creating correlation plots\n");
    plotCorrelation($crosstabTable, undef, $include, 'IRF/plots/correlation/',
                    $sample1, undef, "$sample1 $irfLabel", @irfAxis,
                    $sample2, undef, "$sample2 $irfLabel", @irfAxis,
                    1, [], []);

    status("creating histogram plots\n");
    plotHistogram($crosstabTable, undef, 'IRF/plots/histogram/ratio/',
                  $ratio1Name, undef, $ratio1Label, @ratioAxis, $ratioBinSize,
                  $ratioLines);
    plotHistogram($crosstabTable, undef, 'IRF/plots/histogram/ratio/',
                  $ratio2Name, undef, $ratio2Label, @ratioAxis, $ratioBinSize,
                  $ratioLines);    
    plotPairedHistograms($crosstabTable, undef, $sample1, undef, 
                         $crosstabTable, undef, $sample2, undef, 
                         'retentionFraction', $irfLabel, @irfAxis, $irfBinSize,
                         'IRF/plots/histogram/retentionFraction/', []);        
}

1;

