#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#binSize is required but already defined by VAMP
defined $param{uvMinND} or $param{uvMinND} = 1;  #ninimum normalized density value for allowed genes
defined $param{uvPadding} or $param{uvPadding} = 10000;  #bp to pad onto the indicated end of allowed genes
defined $param{uvPlotDistance} or $param{uvPlotDistance} = 50000;  #bp within the gene to plot
defined $param{uvSpan} or $param{uvSpan} = 15000; #the anticipated span in bp of the tss uv enhancement region
defined $param{uvMinBinHits} or $param{uvMinBinHits} = 10; #the minimum number of combined bin hits required for fracUV calculation (NOT counts)
defined $command{uvMedian} or $command{uvMedian} = ['singleThread', '48:00:00', 2000, 0];

my ($nonUVSample, $uvSample, $geneEnd);

sub uvMedian {
    ($nonUVSample, $uvSample, $geneEnd) = @_;  
#    ($geneEnd eq 'Start' or $geneEnd eq 'End') or die "geneEnd must be eq Start or End\n";
#    my ($minBin, $maxBin) = (-$param{uvPadding}, $param{uvPlotDistance});
#    $geneEnd eq 'End' and ($minBin, $maxBin) = (-$param{uvPlotDistance}, $param{uvPadding});
#    $minBin = int(($minBin / $param{binSize}) + 0.5) * $param{binSize};
#    $maxBin = int(($maxBin / $param{binSize}) + 0.5) * $param{binSize};
#    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $uvHitCount) = getTSSScalars($nonUVSample, $uvSample); 
#    status("  nonUV scalar = $nonUVScalar\n  UV scalar = $uvScalar\n");
#    my $densityTable = createUVDensityTable($nonUVSample, $uvSample, $nonUVScalar, $uvScalar);
#    my $fracUVTable = "uv_fracUV_$geneEnd";
#    dropTable($fracUVTable);
#    runSQL("CREATE TABLE $fracUVTable (name2 varchar2(255), bin number, nRelDens number, uRelDens number, fracUV number)");
#    my $fracUVFile = "$fracUVTable.csv";
#    open my $fracUVFileH, ">", $fracUVFile or die "could not open $fracUVFile: $!"; 
#    
#    foreach my $chrom(1..nChrom()){
#    #foreach my $chrom(21..21){ 
#        
#        status("  chrom $chrom\n");
#        my $hitsTable = createUVHitsTable($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $chrom);
#        my $geneTable = createUVGeneTable($nonUVSample, $chrom);
#        my $joinTable = createUVJoinTable($hitsTable, $geneTable);   
#        my $binSumTable = createUVBinSumTable($joinTable);
#        status("    extracting fracUV data\n");
#        my $fracUVSQL = "SELECT b.name2, b.bin, b.binsStart, b.binsEnd,
#                                (b.nCount/$param{binSize})/d.nDensity nRelDens, 
#                                (b.uCount/$param{binSize})/d.uDensity uRelDens,
#                                case when (b.uCount + b.nCount) > $param{uvMinBinHits}
#                                     then b.uCount/(b.uCount + b.nCount)
#                                     else null end fracUV
#                         FROM $binSumTable b, $densityTable d
#                         WHERE b.name2 = d.name2 ";  
#        runSQL($fracUVSQL, \my($name2,$bin,$binsStart,$binsEnd,$nRelDens,$uRelDens,$fracUV)); 
#        my (%genes, %fracUV);
#        while (fetchRow()){ 
#            $genes{$name2} = [$binsStart,$binsEnd];
#            $fracUV{$name2}{$bin} = [$nRelDens,$uRelDens,$fracUV];    
#        }   
#        foreach my $name2(keys %genes){
#            my ($binsStart,$binsEnd) = @{$genes{$name2}};
#            for (my $bin = $minBin; $bin <= $maxBin; $bin += $param{binSize}){
#                my ($nRelDens,$uRelDens,$fracUV);
#                if($fracUV{$name2}{$bin}){
#                    ($nRelDens,$uRelDens,$fracUV) = @{$fracUV{$name2}{$bin}};                         
#                } else {
#                    my $inData = $bin <= $binsEnd;
#                    $geneEnd eq 'End' and $inData = $bin >= $binsStart;
#                    $inData or next;
#                    ($nRelDens,$uRelDens) = (0, 0);    
#                }
#                $fracUV or $fracUV = -1;
#                print $fracUVFileH "$name2,$bin,$nRelDens,$uRelDens,$fracUV\n";
#            }              
#        }    
#    }
#    status("  loading fracUV data\n");
#    close $fracUVFileH;
#    loadData($fracUVFile, $fracUVTable, ",", "NAME2, BIN, NRELDENS, URELDENS, FRACUV"); 
#    status("  calculating medians\n");  
#    my $medianSQL = "SELECT bin, count(*) nGenes,
#                     percentile_cont(0.05) WITHIN GROUP (ORDER BY nullif(fracUV,-1)) fracUV5, 
#                     percentile_cont(0.5) WITHIN GROUP (ORDER BY nullif(fracUV,-1)) fracUV50,
#                     percentile_cont(0.95) WITHIN GROUP (ORDER BY nullif(fracUV,-1)) fracUV95,
#                     percentile_cont(0.05) WITHIN GROUP (ORDER BY nRelDens) nRelDens5, 
#                     percentile_cont(0.5) WITHIN GROUP (ORDER BY nRelDens) nRelDens50,
#                     percentile_cont(0.95) WITHIN GROUP (ORDER BY nRelDens) nRelDens95,  
#                     percentile_cont(0.05) WITHIN GROUP (ORDER BY uRelDens) uRelDens5, 
#                     percentile_cont(0.5) WITHIN GROUP (ORDER BY uRelDens) uRelDens50,
#                     percentile_cont(0.95) WITHIN GROUP (ORDER BY uRelDens) uRelDens95   
#                     FROM $fracUVTable
#                     GROUP BY bin";            
#    runSQL($medianSQL, \my($bin,$nGenes,$fracUV5,$fracUV50,$fracUV95,$nRelDens5,$nRelDens50,$nRelDens95,$uRelDens5,$uRelDens50,$uRelDens95));
    my $file = "$param{inputPath}/uvMedian.output.$geneEnd.csv";
#    open my $fileH, ">", $file or die "could not open $file: $!\n";
#    print $fileH "bin,nGenes,fracUV5,fracUV50,fracUV95,nRelDens5,nRelDens50,nRelDens95,uRelDens5,uRelDens50,uRelDens95\n";
#    while(fetchRow()){ 
#        print $fileH "$bin,$nGenes,$fracUV5,$fracUV50,$fracUV95,$nRelDens5,$nRelDens50,$nRelDens95,$uRelDens5,$uRelDens50,$uRelDens95\n" 
#    }
#    close $fileH;            
    my ($xMin, $xMax) = (-$param{uvPadding}, $param{uvPlotDistance});
#    $geneEnd eq 'End' and ($xMin, $xMax) = (-$param{uvPlotDistance}, $param{uvPadding});
#    status("  making plots\n"); 
    system("Rscript $param{vampPath}/bin/bruSeq/uvMedianPlot.R $file $xMin $xMax $geneEnd");
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
    my $hitJoinSQL = "SELECT nvl(n.position,u.position) position,
                             $nAgg * $nonUVScalar nCount, $uAgg * $uvScalar uCount
                      FROM ($nonUVSQL) n FULL OUTER JOIN ($uvSQL) u
                        ON (n.position = u.position)";   
    return createUVTable('hits', $hitJoinSQL);
}
sub createUVGeneTable {
    my ($nonUVSample, $chrom) = @_;  
    my $geneSQL = "SELECT name2, strand,
                          decode(strand,1,corrstart_,end_) refPos,    
                          decode(strand,1,corrstart_-$param{uvPadding},corrstart_) geneStart, 
                          decode(strand,1,end_,end_+$param{uvPadding}) geneEnd,
                          decode(strand,1,corrstart_-$param{uvPadding},greatest(corrstart_,end_-$param{uvPlotDistance})) dataStart, 
                          decode(strand,1,least(end_,corrstart_+$param{uvPlotDistance}),end_+$param{uvPadding}) dataEnd,
                          decode(strand,1,1,-1) multiplier
                   FROM refGene_unq_$param{refSeq}
                   WHERE chromosome = $chrom ";      
    $geneEnd eq 'End' and 
    $geneSQL =    "SELECT name2, strand,
                          decode(strand,1,end_,corrstart_) refPos,  
                          decode(strand,1,corrstart_,corrstart_-$param{uvPadding}) geneStart, 
                          decode(strand,1,end_+$param{uvPadding},end_) geneEnd,
                          decode(strand,1,greatest(corrstart_+$param{uvSpan},end_-$param{uvPlotDistance}),corrstart_-$param{uvPadding}) dataStart, 
                          decode(strand,1,end_+$param{uvPadding},least(end_-$param{uvSpan},corrstart_+$param{uvPlotDistance})) dataEnd,
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
    my $joinSQL = "SELECT g.name2, g.refPos, g.dataStart, g.dataEnd, g.multiplier,
                          decode(strand, 1,
                                 (round((g.dataStart-g.refPos) / $param{binSize}) * $param{binSize}) * g.multiplier,
                                 (round((g.dataEnd-g.refPos) / $param{binSize}) * $param{binSize}) * g.multiplier) binsStart,   
                          decode(strand, 1,
                                 (round((g.dataEnd-g.refPos) / $param{binSize}) * $param{binSize}) * g.multiplier,
                                 (round((g.dataStart-g.refPos) / $param{binSize}) * $param{binSize}) * g.multiplier) binsEnd
                   FROM ($geneSQL) g, gmap_$nonUVSample\_unq m
                   WHERE g.name2 = m.name2
                     AND m.normalizedDensity >= $param{uvMinND}
                     AND g.name2 NOT IN ($ovrSQL)";  
    return createUVTable('genes', $joinSQL);
}
sub createUVJoinTable {
    my ($hitsTable, $geneTable) = @_; 
    my $geneJoinSQL = "SELECT (round((h.position-g.refPos) / $param{binSize}) * $param{binSize}) * g.multiplier bin,   
                              h.nCount, h.uCount, g.name2, g.binsStart, g.binsEnd   
                       FROM $hitsTable h, $geneTable g
                       WHERE h.position BETWEEN g.dataStart - $param{binSize} 
                                            AND g.dataEnd   + $param{binSize} ";   
    return createUVTable('join', $geneJoinSQL);
}
sub createUVBinSumTable {
    my ($joinTable) = @_;
    my $binSumSQL = "SELECT name2, bin, binsStart, binsEnd, sum(uCount) uCount, sum(nCount) nCount
                     FROM $joinTable
                     GROUP BY name2, bin, binsStart, binsEnd";      
    return createUVTable('binSum', $binSumSQL);                               
}  
sub createUVTable {
    my ($tableSuffix, $sql) = @_;  
    status("    creating $tableSuffix table\n");
    my $outTable = "uv_$tableSuffix\_$geneEnd";                                                                  
    dropTable($outTable);                 
    runSQL("CREATE TABLE $outTable AS $sql"); 
    return $outTable;      
}

1;

