#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#binSize is required but already defined by VAMP
defined $param{minND} or $param{minND} = 0.3;  #minimum normalized density value for allowed genes
defined $param{uvSpan} or $param{uvSpan} = 5000;  #region after TSS used in fracUV calculation
defined $param{padding} or $param{padding} = 2500;  #required non-overlap region before and after uvSpan
defined $command{uvRepression} or $command{uvRepression} = ['singleThread', '48:00:00', 2000, 0];

my ($mergedSample);

sub uvRepression {
    my ($nonUVSample, $uvSample) = @_;  
    $mergedSample = "$nonUVSample$uvSample";    
    checkUvrDirs();
    my $uvrTable = calculateUvr($nonUVSample, $uvSample);
    plotUvrHistogram($uvrTable);
    printUvrTable($uvrTable);
}
sub checkUvrDirs {
    my @dirs = qw(  plots
                    tables   );
    foreach my $dir(@dirs) { 
        my $path = "$param{inputPath}/$dir";
        -d $path or mkdir $path;
    }
}
sub calculateUvr {
    my ($nonUVSample, $uvSample) = @_;  
    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $uvHitCount) = getTSSScalars($nonUVSample, $uvSample); 
    my $uvrTable = "UVR_$nonUVSample\_$uvSample";
    dropTable($uvrTable);
    runSQL("CREATE TABLE $uvrTable (name2 varchar2(255), nDensity number, geneSize number,
                                    uCount number, nCount number, fracUV number)");    
    my $uvrFile = "$uvrTable.csv";
    open my $uvrFileH, ">", $uvrFile or die "could not open $uvrFile: $!"; 
    foreach my $chrom(1..nChrom()){
    #foreach my $chrom(19..19){     
        status("  chrom $chrom\n");
        my $hitsTable = createUvrHitsTable($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $chrom);
        my $geneTable = createUvrGeneTable($nonUVSample, $chrom);
        my $joinSQL = "SELECT g.name2, g.nDensity, g.geneSize, sum(h.uCount) uCount, sum(h.nCount) nCount
                          FROM $hitsTable h, $geneTable g
                          WHERE h.position BETWEEN g.dataStart AND g.dataEnd
                          GROUP BY g.name2, g.nDensity, g.geneSize";        
        my $fracUVSql = "SELECT name2, nDensity, geneSize, uCount, nCount, uCount/nvl(uCount + nCount,0) fracUV
                         FROM ($joinSQL)";
        status("    extracting fracUV data\n");                         
        runSQL($fracUVSql,\my($name2,$nDens,$geneSize,$uCount,$nCount,$fracUV));             
        while (fetchRow()){ 
            $fracUV or $fracUV = "";
            print $uvrFileH "$name2,$nDens,$geneSize,$uCount,$nCount,$fracUV\n";
        }     
        dropTable($hitsTable, $geneTable);    
    }
    status("  loading fracUV data\n");
    close $uvrFileH;
    loadData($uvrFile, $uvrTable, ",", "NAME2, NDENSITY, GENESIZE, UCOUNT, NCOUNT, FRACUV"); 
    return $uvrTable;
}
sub createUvrHitsTable {
    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $chrom) = @_;
    my ($nAgg, $uAgg) = ("decode(n.count_,null,0,1)", "decode(u.count_,null,0,1)");
    $param{keepDups} and ($nAgg, $uAgg) = ("nvl(n.count_,0)", "nvl(u.count_,0)");
    my $nonUVSQL = "SELECT * FROM $nonUVHitsTable WHERE chromosome = $chrom";
    my $uvSQL = "SELECT * FROM $uvHitsTable WHERE chromosome = $chrom";
    my $hitJoinSQL =  "SELECT nvl(n.position,u.position) position, 
                              $nAgg * $nonUVScalar nCount, $uAgg * $uvScalar uCount
                       FROM ($nonUVSQL) n FULL OUTER JOIN ($uvSQL) u
                         ON (n.position = u.position AND n.strand = u.strand)";                    
    return createUvrTable('hits', $hitJoinSQL);
}
sub createUvrGeneTable {
    my ($nonUVSample, $chrom) = @_;  
    my $tssSQL = "SELECT name2, strand, corrstart_, end_, end_ - corrstart_ + 1 geneSize, decode(strand,1,corrstart_,end_) tss
                   FROM refGene_unq_$param{refSeq}
                   WHERE chromosome = $chrom "; #ensure that all return genes sample the full uvSpan
    my $dataSQL = "SELECT name2, tss, geneSize,
                   decode(strand,1,tss,greatest(tss-$param{uvSpan},corrstart_)) dataStart,
                   decode(strand,1,least(tss+$param{uvSpan},end_),tss) dataEnd
                   FROM ($tssSQL)
                   WHERE geneSize >= $param{uvSpan}";                
    my $geneSQL = "SELECT name2, tss, geneSize, dataStart, dataEnd, 
                   dataStart-$param{padding} spanStart,
                   dataEnd+$param{padding} spanEnd
                   FROM ($dataSQL)";  
    my $ovrSQL = "SELECT g1.name2
                  FROM ($geneSQL) g1, ($geneSQL) g2
                  WHERE g1.spanStart <= g2.spanEnd
                    AND g2.spanStart <= g1.spanEnd
                    AND g1.name2 != g2.name2";        
    my $joinSQL = "SELECT g.name2, m.normalizedDensity nDensity, g.tss, g.geneSize, g.dataStart, g.dataEnd
                   FROM ($geneSQL) g, gmap_$nonUVSample\_unq m
                   WHERE g.name2 = m.name2
                     AND m.normalizedDensity >= $param{minND}
                     AND g.name2 NOT IN ($ovrSQL)";  
    return createUvrTable('genes', $joinSQL);
}
sub createUvrTable {
    my ($tableSuffix, $sql) = @_;  
    status("    creating $tableSuffix table\n");
    my $outTable = "uvr$tableSuffix$mergedSample";                                                                  
    dropTable($outTable);                 
    runSQL("CREATE TABLE $outTable AS $sql"); 
    return $outTable;      
}
sub plotUvrHistogram {
    my ($uvrTable) = @_; 
    my $include = ['name2'];
    my $ndMax = 1000; 
    my @ndAxis = (1/$ndMax, $ndMax, 1);
    my $ndBinSize = 1/3;
    my $ndLabel = "Control Normalized Density";
    my @fractionAxis = (0, 1, undef);  
    my $fractionLabel = "Fraction UV";
    my $fractionBinSize = 0.05;
    my $foldInduction = 2;
    my $fractionLines = [1/(1+$foldInduction), 1/2, $foldInduction/(1+$foldInduction)];
    status("creating plots\n");
    plotCorrelation($uvrTable, undef, $include, "plots/$uvrTable",
                    'nDensity', undef, $ndLabel, @ndAxis,
                    'fracUV', undef, $fractionLabel, @fractionAxis,
                    undef, $fractionLines, []);                    
    plotHistogram($uvrTable, undef, "plots/$uvrTable",
                  'fracUV', undef, $fractionLabel, @fractionAxis, $fractionBinSize,
                  $fractionLines);    
}
sub printUvrTable { 
    my ($uvrTable) = @_;
    my $sql = "SELECT name2, geneSize, nDensity, uCount, nCount, fracUV
               FROM $uvrTable
               ORDER BY fracUV";
    return printPlotData($sql, "tables/$uvrTable", 'Gene', 'Size', 'Control ND', 'UV Count','Control Count','Fraction UV');
}

1;


