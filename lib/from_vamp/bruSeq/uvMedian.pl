#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#binSize is required but already defined by VAMP
defined $param{minND} or $param{minND} = 1;  #minimum normalized density value for allowed genes
defined $param{padding} or $param{padding} = 10000;  #bp to pad onto the indicated end of allowed genes
defined $param{uvPlotDistance} or $param{uvPlotDistance} = 50000;  #bp within the gene to plot
defined $param{uvSpan} or $param{uvSpan} = 15000; #the anticipated span in bp of the tss uv enhancement region
defined $param{combineUVStrands} or $param{combineUVStrands} = 0; #whether or not to maintain strand information
defined $param{geneEnd} or $param{geneEnd} = 'Start'; #whether or not to maintain strand information
defined $param{minBinHits} or $param{minBinHits} = 10; #the minimum number of combined bin hits required for fracUV calculation (NOT counts)
defined $command{uvMedian} or $command{uvMedian} = ['singleThread', '48:00:00', 2000, 0];

my ($nonUVSample, $uvSample, $mergedSample);

#P=padding, RP=refPos, UVP=uvPlotDistance
#FORWARD GENE:
#              Start                                End
#   ------------->===================================X--------------
#       |--------|-------------|       |-------------|--------|
#         P      RP    UVP                   UVP     RP    P
#REVERSE GENE:
#             End                                  Start
#   -----------X=====================================<--------------
#     |--------|-------------|         |-------------|--------|
#        P     RP    UVP                   UVP       RP    P

sub uvMedian {
    ($nonUVSample, $uvSample) = @_;
    $mergedSample = "$nonUVSample$uvSample";
    checkUvMedianDirs();
    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $uvHitCount) = getTSSScalars($nonUVSample, $uvSample); 
    status("  nonUV scalar = $nonUVScalar\n  UV scalar = $uvScalar\n");
    my $densityTable = createUVDensityTable($nonUVSample, $uvSample, $nonUVScalar, $uvScalar); 
    my $fracUVTable = "uvfracUV$param{geneEnd}$mergedSample$param{combineUVStrands}";
    dropTable($fracUVTable);
    runSQL("CREATE TABLE $fracUVTable (name2 varchar2(255), strand number, bin number, nRelDens number, uRelDens number, fracUV number)");
    my $fracUVFile = "$fracUVTable.csv";
    open my $fracUVFileH, ">", $fracUVFile or die "could not open $fracUVFile: $!"; 
    foreach my $chrom(1..nChrom()){
    #foreach my $chrom(1..1){
        status("  chrom $chrom\n");
        my $hitsTable = createUVHitsTable($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $chrom);
        my $geneTable = createUVGeneTable($nonUVSample, $chrom);
        my $joinTable = createUVJoinTable($hitsTable, $geneTable);   
        my $binSumTable = createUVBinSumTable($joinTable);
        status("    extracting fracUV data\n");
        my $fracUVSQL = "SELECT b.name2, b.bin, b.binsStart, b.binsEnd, strand,
                                (b.nCount/$param{binSize})/d.nDensity nRelDens, 
                                (b.uCount/$param{binSize})/d.uDensity uRelDens,
                                case when (abs(b.uCount) + abs(b.nCount)) > $param{minBinHits}
                                     then abs(b.uCount)/(abs(b.uCount) + abs(b.nCount))
                                     else null end fracUV
                         FROM $binSumTable b, $densityTable d
                         WHERE b.name2 = d.name2 ";  
        runSQL($fracUVSQL, \my($name2,$bin,$binsStart,$binsEnd,$strand,$nRelDens,$uRelDens,$fracUV)); 
        my (%genes, %fracUV);
        while (fetchRow()){ 
            $genes{$name2} = [$binsStart,$binsEnd];
            $fracUV{$name2}{$bin}{$strand} = [$nRelDens,$uRelDens,$fracUV];
        }   
        foreach my $name2(keys %genes){
            my ($binsStart,$binsEnd) = @{$genes{$name2}};
            for (my $bin = $binsStart; $bin <= $binsEnd; $bin += $param{binSize}){
                foreach my $strand(-1,1){
                    my ($nRelDens,$uRelDens,$fracUV);
                    if($fracUV{$name2}{$bin} and $fracUV{$name2}{$bin}{$strand}){
                        ($nRelDens,$uRelDens,$fracUV) = @{$fracUV{$name2}{$bin}{$strand}};                         
                    } else {
                        ($nRelDens,$uRelDens) = (0, 0);    
                    }
                    $fracUV or $fracUV = -99;
                    print $fracUVFileH "$name2,$bin,$strand,$nRelDens,$uRelDens,$fracUV\n";
                }
            }              
        }    
        dropTable($hitsTable, $geneTable, $joinTable, $binSumTable);    
    }
    status("  loading fracUV data\n");
    close $fracUVFileH;
    loadData($fracUVFile, $fracUVTable, ",", "NAME2, BIN, STRAND, NRELDENS, URELDENS, FRACUV"); 
    status("  calculating means and medians\n");
    my $aggSql = "SELECT strand, bin, count(*) nGenes,
                 percentile_cont(0.05) WITHIN GROUP (ORDER BY nullif(fracUV,-99)) fracUV5, 
                 percentile_cont(0.5) WITHIN GROUP (ORDER BY nullif(fracUV,-99)) fracUV50,
                 percentile_cont(0.95) WITHIN GROUP (ORDER BY nullif(fracUV,-99)) fracUV95,
                 percentile_cont(0.05) WITHIN GROUP (ORDER BY nRelDens) nRelDens5, 
                 percentile_cont(0.5) WITHIN GROUP (ORDER BY nRelDens) nRelDens50,
                 percentile_cont(0.95) WITHIN GROUP (ORDER BY nRelDens) nRelDens95,
                 percentile_cont(0.05) WITHIN GROUP (ORDER BY uRelDens) uRelDens5, 
                 percentile_cont(0.5) WITHIN GROUP (ORDER BY uRelDens) uRelDens50,
                 percentile_cont(0.95) WITHIN GROUP (ORDER BY uRelDens) uRelDens95,
                 avg(nRelDens) nRelDensMean,
                 avg(uRelDens) uRelDensMean
                 FROM $fracUVTable
                 GROUP BY strand, bin
                 ORDER BY strand, bin";
    runSQL($aggSql, \my($strand,$bin,$nGenes,$fracUV5,$fracUV50,$fracUV95,
                                             $nRelDens5,$nRelDens50,$nRelDens95,
                                             $uRelDens5,$uRelDens50,$uRelDens95,
                                             $nRelDensMean,$uRelDensMean));
    my $strandType = "splitStrands";
    $param{combineUVStrands} and $strandType = "combinedStrands";
    my $file = "$param{inputPath}/$param{geneEnd}/$strandType/uvMedian.$nonUVSample\_$uvSample.csv";
    open my $fileH, ">", $file or die "could not open $file: $!\n";
    print $fileH "strand,bin,nGenes,fracUV5,fracUV50,fracUV95,".
                                   "nRelDens5,nRelDens50,nRelDens95,".
                                   "uRelDens5,uRelDens50,uRelDens95,".
                                   "nRelDensMean,uRelDensMean\n";
    while(fetchRow()){ 
       print $fileH "$strand,$bin";
       foreach my $val($nGenes,$fracUV5,$fracUV50,$fracUV95,
                             $nRelDens5,$nRelDens50,$nRelDens95,
                             $uRelDens5,$uRelDens50,$uRelDens95,
                             $nRelDensMean,$uRelDensMean){
           defined $val or $val = "";
           print $fileH ",$val";
       }
       print $fileH "\n";
    }
    close $fileH;            
    my ($xMin, $xMax) = (-$param{padding}, $param{uvPlotDistance});
    $param{geneEnd} eq 'End' and ($xMin, $xMax) = (-$param{uvPlotDistance}, $param{padding});
    status("  making plots\n"); 
    system("Rscript $param{vampPath}/bin/bruSeq/uvMedianPlot.R $file $xMin $xMax $param{geneEnd} $param{combineUVStrands}");
    #dropTable($densityTable, $fracUVTable);
    dropTable($densityTable);
}
sub checkUvMedianDirs {
    my @dirs = qw(  Start
                    Start/splitStrands
                    Start/combinedStrands
                    End
                    End/splitStrands
                    End/combinedStrands
    );
    foreach my $dir(@dirs) {
        my $path = "$param{inputPath}/$dir";
        -d $path or mkdir $path;
    }
}
sub createUVDensityTable {
    my ($nonUVSample, $uvSample, $nonUVScalar, $uvScalar) = @_;  
    my $densSQL = "SELECT n.name2, n.density*$nonUVScalar nDensity, u.density*$uvScalar uDensity
                   FROM gmap_$nonUVSample\_unq n, gmap_$uvSample\_unq u
                   WHERE n.name2 = u.name2";
    return createUVTable('density', $densSQL);
}
sub createUVHitsTable {
    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $chrom) = @_;
    my ($nAgg, $uAgg) = ("decode(n.count_,null,0,1)", "decode(u.count_,null,0,1)");
    $param{keepDups} and ($nAgg, $uAgg) = ("nvl(n.count_,0)", "nvl(u.count_,0)");
    my $nonUVSQL = "SELECT * FROM $nonUVHitsTable WHERE chromosome = $chrom";
    my $uvSQL = "SELECT * FROM $uvHitsTable WHERE chromosome = $chrom";
    my $hitJoinSQL =  "SELECT nvl(n.position,u.position) position, nvl(n.strand,u.strand) strand,
                              $nAgg * $nonUVScalar nCount, $uAgg * $uvScalar uCount
                       FROM ($nonUVSQL) n FULL OUTER JOIN ($uvSQL) u
                         ON (n.position = u.position AND n.strand = u.strand)";                    
    return createUVTable('hits', $hitJoinSQL);
}
sub createUVGeneTable {
    my ($nonUVSample, $chrom) = @_;  
    my $geneSQL = "SELECT name2, strand,
                          decode(strand,1,corrstart_,end_) refPos,    
                          decode(strand,1,corrstart_-$param{padding},corrstart_) geneStart, 
                          decode(strand,1,end_,end_+$param{padding}) geneEnd,
                          decode(strand,1,corrstart_-$param{padding},greatest(corrstart_,end_-$param{uvPlotDistance})) dataStart, 
                          decode(strand,1,least(end_,corrstart_+$param{uvPlotDistance}),end_+$param{padding}) dataEnd,
                          decode(strand,1,1,-1) multiplier
                    FROM refGene_unq_$param{refSeq}
                   WHERE chromosome = $chrom ";
    $param{geneEnd} eq 'End' and
    $geneSQL =    "SELECT name2, strand,
                          decode(strand,1,end_,corrstart_) refPos,
                          decode(strand,1,corrstart_,corrstart_-$param{padding}) geneStart,
                          decode(strand,1,end_+$param{padding},end_) geneEnd,
                          decode(strand,1,greatest(corrstart_+$param{uvSpan},end_-$param{uvPlotDistance}),corrstart_-$param{padding}) dataStart,
                          decode(strand,1,end_+$param{padding},least(end_-$param{uvSpan},corrstart_+$param{uvPlotDistance})) dataEnd,
                          decode(strand,1,1,-1) multiplier
                   FROM refGene_unq_$param{refSeq}
                   WHERE chromosome = $chrom
                     AND ((strand = 1 AND corrstart_ + $param{uvSpan} < end_)
                      OR  (strand = 2 AND end_ - $param{uvSpan} > corrstart_))"; #enforce gene size limit on end to prevent start from interfering
    my $ovrSQL = "SELECT g1.name2
                  FROM ($geneSQL) g1, ($geneSQL) g2
                  WHERE g1.geneStart <= g2.geneEnd
                    AND g2.geneStart <= g1.geneEnd
                    AND g1.name2 != g2.name2";
    my $joinSQL = "SELECT g.name2, g.strand, g.refPos, g.dataStart, g.dataEnd, g.multiplier,
                          round((decode(strand,1,g.dataStart,g.dataEnd)-g.refPos) / $param{binSize}) * $param{binSize} * g.multiplier binsStart,
                          round((decode(strand,1,g.dataEnd,g.dataStart)-g.refPos) / $param{binSize}) * $param{binSize} * g.multiplier binsEnd
                   FROM ($geneSQL) g, gmap_$nonUVSample\_unq m
                   WHERE g.name2 = m.name2
                     AND m.normalizedDensity >= $param{minND}
                     AND g.name2 NOT IN ($ovrSQL)";  
    return createUVTable('genes', $joinSQL);
}
sub createUVJoinTable {
    my ($hitsTable, $geneTable) = @_;  
    my $strand = 1;             
    $param{stranded} and $strand = "case when g.strand = h.strand then 1 else -1 end";
    my $geneJoinSQL = "SELECT round((h.position-g.refPos) / $param{binSize}) * $param{binSize} * g.multiplier bin,
                              $strand strand, h.nCount, h.uCount, g.name2, g.binsStart, g.binsEnd   
                      FROM $hitsTable h, $geneTable g
                      WHERE h.position BETWEEN g.dataStart - $param{binSize} 
                                           AND g.dataEnd   + $param{binSize} ";                                                                                  
    return createUVTable('join', $geneJoinSQL);
}
sub createUVBinSumTable {
    my ($joinTable) = @_;
    my $binSumSQL;
    if ($param{combineUVStrands}){
        $binSumSQL = "SELECT name2, bin, binsStart, binsEnd, 1 strand, sum(strand * uCount) uCount, sum(strand * nCount) nCount
                      FROM $joinTable
                      GROUP BY name2, bin, binsStart, binsEnd";       
    } else {
        $binSumSQL = "SELECT name2, bin, binsStart, binsEnd, strand, sum(strand * uCount) uCount, sum(strand * nCount) nCount
                      FROM $joinTable
                      GROUP BY name2, bin, binsStart, binsEnd, strand";                      
    }              
    return createUVTable('binSum', $binSumSQL);                               
}  
sub createUVTable {
    my ($tableSuffix, $sql) = @_;  
    status("    creating $tableSuffix table\n");
    my $outTable = "uv$tableSuffix$param{geneEnd}$mergedSample";
    dropTable($outTable);                 
    runSQL("CREATE TABLE $outTable AS $sql"); 
    return $outTable;      
}

1;

