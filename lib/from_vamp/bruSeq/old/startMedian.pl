#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#binSize is required but already defined by VAMP

defined $param{startMinND} or $param{startMinND} = 1; 
defined $param{startPadding} or $param{startPadding} = 10000;
defined $command{startMedian} or $command{startMedian} = ['singleThread', '48:00:00', 2000, 0];

sub startMedian {
    my ($nonUVSample, $uvSample) = @_;  
    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $uvHitCount) = getTSSScalars($nonUVSample, $uvSample); 
    status("  nonUV scalar = $nonUVScalar\n  UV scalar = $uvScalar\n");
    my $fracUVTable = 'start_fracUV';
    my $binSumTable = 'start_binSum';
    my $densityTable = createStartDensityTable($nonUVSample, $uvSample, $nonUVScalar, $uvScalar);
    dropTable($fracUVTable);
    runSQL("CREATE TABLE $fracUVTable (name2 varchar2(255), bin number, nRelDens number, uRelDens number, fracUV number)");
    my $fracUVFile = "$fracUVTable.csv";
    open my $fracUVFileH, ">", $fracUVFile or die "could not open $fracUVFile: $!"; 
    
    #foreach my $chrom(1..nChrom()){
    foreach my $chrom(21..21){ 
        
        status("  chrom $chrom\n");
        my $hitsTable = createTssHitsTable($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $chrom);
        my $geneTable = createStartGeneTable($nonUVSample, $chrom);
        my $joinTable = createStartJoinTable($hitsTable, $geneTable);                
        my $binSumSQL = "SELECT name2, bin, sum(uCount) uCount, sum(nCount) nCount
                         FROM $joinTable
                         GROUP BY name2, bin"; 
        dropTable($binSumTable);          
        runSQL("CREATE TABLE $binSumTable AS $binSumSQL");
        my $fracUVSQL = "SELECT b.name2, b.bin, 
                                (b.nCount/$param{binSize})/d.nDensity nRelDens, 
                                (b.uCount/$param{binSize})/d.uDensity uRelDens,
                                 b.uCount/(b.uCount + b.nCount) fracUV
                         FROM $binSumTable b, $densityTable d
                         WHERE b.name2 = d.name2 ";  
        runSQL($fracUVSQL, \my($name2,$bin,$nRelDens,$uRelDens,$fracUV));
        while (fetchRow()){ print $fracUVFileH "$name2,$bin,$nRelDens,$uRelDens,$fracUV\n" }
    }
    close $fracUVFileH;
    loadData($fracUVFile, $fracUVTable, ",", "NAME2, BIN, NRELDENS, URELDENS, FRACUV");   
    my $medianSQL = "SELECT bin, count(*) nGenes,
                     percentile_cont(0.05) WITHIN GROUP (ORDER BY fracUV) fracUV5, 
                     percentile_cont(0.5) WITHIN GROUP (ORDER BY fracUV) fracUV50,
                     percentile_cont(0.95) WITHIN GROUP (ORDER BY fracUV) fracUV95,
                     percentile_cont(0.05) WITHIN GROUP (ORDER BY nRelDens) nRelDens5, 
                     percentile_cont(0.5) WITHIN GROUP (ORDER BY nRelDens) nRelDens50,
                     percentile_cont(0.95) WITHIN GROUP (ORDER BY nRelDens) nRelDens95,  
                     percentile_cont(0.05) WITHIN GROUP (ORDER BY uRelDens) uRelDens5, 
                     percentile_cont(0.5) WITHIN GROUP (ORDER BY uRelDens) uRelDens50,
                     percentile_cont(0.95) WITHIN GROUP (ORDER BY uRelDens) uRelDens95   
                     FROM $fracUVTable
                     GROUP BY bin";            
    runSQL($medianSQL, \my($bin,$nGenes,$fracUV5,$fracUV50,$fracUV95,$nRelDens5,$nRelDens50,$nRelDens95,$uRelDens5,$uRelDens50,$uRelDens95));
    my $file = "$param{inputPath}/startMedian.output.csv";
    open my $fileH, ">", $file or die "could not open $file: $!\n";
    print $fileH "bin,nGenes,fracUV5,fracUV50,fracUV95,nRelDens5,nRelDens50,nRelDens95,uRelDens5,uRelDens50,uRelDens95\n";
    while(fetchRow()){ 
        print $fileH "$bin,$nGenes,$fracUV5,$fracUV50,$fracUV95,$nRelDens5,$nRelDens50,$nRelDens95,$uRelDens5,$uRelDens50,$uRelDens95\n" 
    }
    close $fileH;         
    system("Rscript $param{vampPath}/bin/bruSeq/startFracUV_2.R $file $param{maxGeneDistance} $param{startPadding}");
    system("Rscript $param{vampPath}/bin/bruSeq/startCounts_2.R $file $param{maxGeneDistance} $param{startPadding}");
}

sub createStartDensityTable {
    my ($nonUVSample, $uvSample, $nonUVScalar, $uvScalar) = @_;  
    my $densSQL = "SELECT n.name2, n.density*$nonUVScalar nDensity, u.density*$uvScalar uDensity
                   FROM gmap_$nonUVSample\_unq n, gmap_$uvSample\_unq u
                   WHERE n.name2 = u.name2";
    my $outTable = 'start_density';                                   
    dropTable($outTable);
    runSQL("CREATE TABLE $outTable AS $densSQL");
    return $outTable;  
}
sub createStartGeneTable {
    my ($nonUVSample, $chrom) = @_;         
    my $geneSQL = "SELECT name2, 
                          decode(strand,1,corrstart_,end_) refPos,    
                          decode(strand,1,corrstart_-$param{startPadding},corrstart_) geneStart, 
                          decode(strand,1,end_,end_+$param{startPadding}) geneEnd,
                          decode(strand,1,corrstart_-$param{startPadding},greatest(corrstart_,end_-$param{maxGeneDistance})) dataStart, 
                          decode(strand,1,least(end_,corrstart_+$param{maxGeneDistance}),end_+$param{startPadding}) dataEnd,
                          decode(strand,1,1,-1) multiplier
                   FROM refGene_unq_$param{refSeq}
                   WHERE chromosome = $chrom ";            
    my $ovrSQL = "SELECT g1.name2
                  FROM ($geneSQL) g1, ($geneSQL) g2
                  WHERE g1.geneStart <= g2.geneEnd
                    AND g2.geneStart <= g1.geneEnd
                    AND g1.name2 != g2.name2";    
    my $joinSQL = "SELECT g.name2, g.refPos, g.dataStart, g.dataEnd, g.multiplier
                   FROM ($geneSQL) g, gmap_$nonUVSample\_unq m
                   WHERE g.name2 = m.name2
                     AND m.normalizedDensity >= $param{startMinND}
                     AND g.name2 NOT IN ($ovrSQL)";  
    my $outTable = 'start_genes';                                   
    dropTable($outTable);
    runSQL("CREATE TABLE $outTable AS $joinSQL");
    return $outTable;
}
sub createStartJoinTable {
    my ($tssHitsTable, $tssGeneTable) = @_; 
    my $geneJoinSQL = "SELECT (round((h.position-g.refPos) / $param{binSize}) * $param{binSize}) * multiplier bin,   
                              h.nCount, h.uCount, g.name2
                       FROM $tssHitsTable h, $tssGeneTable g
                       WHERE h.position BETWEEN g.dataStart - $param{binSize} 
                                            AND g.dataEnd   + $param{binSize} ";   
    my $outTable = 'start_join';                                                                    
    dropTable($outTable);                 
    runSQL("CREATE TABLE $outTable AS $geneJoinSQL"); 
    return $outTable;     
}


1;


