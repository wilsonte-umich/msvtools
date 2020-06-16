#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
defined $param{geneType} or $param{geneType} = "Unq"; 
defined $param{minND} or $param{minND} = 0.1;  #minimum normalized density value for allowed genes
defined $command{calculateIRF} or $command{calculateIRF} = ['multiThread', '1:00:00', 1000, 0];
defined $command{calculateERF} or $command{calculateERF} = ['multiThread', '1:00:00', 1000, 0];
    
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
    my ($frfTable, $featureCoverage, $featureDensity) = createFRFTable($sample, $featureType, $mapType, $tableRoot);
    checkFRFDirs($featureType, $tableRoot);
    printFRFTables($featureType, $tableRoot, $frfTable, $featureCoverage, $featureDensity);
}
sub createFRFTable {
    my ($sample, $featureType, $mapType, $tableRoot) = @_; 
    status("calculating $sample $featureType retention fraction\n");  
    my $fmapTable = getTableName($mapType, $sample);
    my $gmapTable = getTableName('GMap', "$sample\_$param{geneType}");  
    my $frfTable = getTableName($tableRoot, "$sample\_$param{geneType}");
    my $featureDensity = $featureType."Density";
    my $featureCoverage = $featureType."Coverage";
    my $joinSql = "SELECT f.name2, f.chromosome, f.start_, f.end_, 
                          f.coverage $featureCoverage, f.density $featureDensity, 
                          g.coverage geneCoverage, g.density geneDensity, g.normalizedDensity geneNormalizedDensity,
                          f.density/nullif(g.density,0) $tableRoot
                   FROM $fmapTable f, $gmapTable g
                   WHERE f.name2 = g.name2";            
    dropTable($frfTable);      
    runSQL("CREATE TABLE $frfTable AS $joinSql");
    return ($frfTable, $featureCoverage, $featureDensity);
}
sub checkFRFDirs {
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
    my ($featureType, $tableRoot, $frfTable, $featureCoverage, $featureDensity) = @_;
    status("printing output tables\n");  
    my $printSql = "SELECT sum(1) over (order by $tableRoot DESC rows unbounded preceding) RANK, 
                           name2, chromosome, start_, end_, end_ - start_ + 1 size_,
                           $featureCoverage, $featureDensity, geneCoverage, geneDensity, geneNormalizedDensity, $tableRoot
                    FROM $frfTable
                    WHERE geneNormalizedDensity > $param{minND}
                    ORDER BY $tableRoot DESC";
    printPlotData($printSql, "$tableRoot/tables/by_$featureType/by_sample/$frfTable", 
                         'RANK','GENE','CHROMOSOME','START','END','SIZE',
                          $featureCoverage, $featureDensity,'geneCoverage','geneDensity','geneNormalizedDensity',$tableRoot);
 
    my $geneTable = "refgene_$param{geneType}_$param{refSeq}";
    $printSql = "SELECT g.name2, g.chromosome, g.corrstart_ start_, g.end_, g.genesize, g.mrnasize,
                        f.geneCoverage, f.geneDensity, f.geneNormalizedDensity,
                        min($featureDensity) min_$featureDensity, max($featureDensity) max_$featureDensity,
                        min($tableRoot) min_$tableRoot, max($tableRoot) max_$tableRoot
                 FROM $frfTable f, $geneTable g
                 WHERE f.name2 = g.name2
                   AND f.geneNormalizedDensity > $param{minND}
                 GROUP BY g.name2, g.chromosome, g.corrstart_, g.end_, g.genesize, g.mrnasize,
                          f.geneCoverage, f.geneDensity, f.geneNormalizedDensity";   
    $printSql = "SELECT sum(1) over (order by p.max_$tableRoot rows unbounded preceding) RANK, p.* 
                 FROM ($printSql) p
                 ORDER BY p.max_$tableRoot";      
    printPlotData($printSql, "$tableRoot/tables/by_gene/by_sample/$frfTable", 
                         'RANK','GENE','CHROMOSOME','START','END','GENESIZE','MRNASIZE',
                         'geneCoverage','geneDensity','geneNormalizedDensity',
                         "min_$featureDensity", "max_$featureDensity",
                         "min_$tableRoot", "max_$tableRoot");                                
}

1;

