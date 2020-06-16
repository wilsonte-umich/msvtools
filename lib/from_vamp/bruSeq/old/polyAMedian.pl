#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#binSize is required but already defined by VAMP
defined $command{polyAMedian} or $command{polyAMedian} = ['singleThread', '48:00:00', 2000, 0];

sub polyAMedian {
    my ($nonUVSample, $uvSample) = @_;  
    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $uvHitCount) = getTSSScalars($nonUVSample, $uvSample); 
    status("  nonUV scalar = $nonUVScalar\n  UV scalar = $uvScalar\n");
    my $fracUVTable = 'tss_fracUV';
    my $binSumTable = 'tss_binSum';
    my $geneSumTable = 'tss_geneCount';
    dropTable($fracUVTable);
    runSQL("CREATE TABLE $fracUVTable (name2 varchar2(255), strand number, bin number, nCount number, uCount number, fracUV number)");
    my $fracUVFile = "$fracUVTable.csv";
    open my $fracUVFileH, ">", $fracUVFile or die "could not open $fracUVFile: $!";  
    my $tssDensityTable = createPolyADensityTable($nonUVSample, $uvSample);
    
    foreach my $chrom(1..nChrom()){
    #foreach my $chrom(21..22){
    
        status("  chrom $chrom\n");
        my $tssDensityTable = createPolyADensityTable($nonUVSample, $uvSample, $chrom);
        my $tssHitsTable = createTssHitsTable($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $chrom);
        my $tssGeneTable = createPolyAGeneTable($chrom);
        my $tssJoinTable = createPolyAJoinTable($tssHitsTable, $tssGeneTable);                                                  
        my $binSumSQL = "SELECT name2, strand, bin, sum(uCount) uCount, sum(nCount) nCount
                         FROM $tssJoinTable
                         GROUP BY name2, strand, bin";
        dropTable($binSumTable);          
        runSQL("CREATE TABLE $binSumTable AS $binSumSQL");
        my $fracUVSQL = "SELECT b.name2, b.strand, b.bin, 
                                 (b.nCount/$param{binSize})/d.nDensity nCount, 
                                 (b.uCount/$param{binSize})/d.uDensity uCount,
                                 b.uCount/(b.uCount + b.nCount) fracUV
                         FROM $binSumTable b, $tssDensityTable d
                         WHERE b.name2 = d.name2
                           AND (b.uCount > $param{minBinHits} OR b.nCount > $param{minBinHits}) ";  
        runSQL($fracUVSQL, \my($name2,$strand,$bin,$nCount,$uCount,$fracUV));
        while (fetchRow()){ print $fracUVFileH "$name2,$strand,$bin,$nCount,$uCount,$fracUV\n" }
    }
    close $fracUVFileH;
    loadData($fracUVFile, $fracUVTable, ",", "NAME2, STRAND, BIN, NCOUNT, UCOUNT, FRACUV");   
    my $medianSQL = "SELECT bin, count(*) nGenes,
                     percentile_cont(0.05) WITHIN GROUP (ORDER BY fracUV) fracUV5, 
                     percentile_cont(0.5) WITHIN GROUP (ORDER BY fracUV) fracUV50,
                     percentile_cont(0.95) WITHIN GROUP (ORDER BY fracUV) fracUV95,
                     percentile_cont(0.5) WITHIN GROUP (ORDER BY nCount) nCount50,
                     percentile_cont(0.5) WITHIN GROUP (ORDER BY uCount) uCount50
                     FROM $fracUVTable
                     GROUP BY bin";            
    runSQL($medianSQL, \my($bin,$nGenes,$fracUV5,$fracUV50,$fracUV95,$nCount50,$uCount50));
    my $file = "$param{inputPath}/polyAMedian.output.csv";
    open my $fileH, ">", $file or die "could not open $file: $!\n";
    print $fileH "bin,nGenes,fracUV5,fracUV50,fracUV95,nCount50,uCount50\n";
    while(fetchRow()){ 
        print $fileH "$bin,$nGenes,$fracUV5,$fracUV50,$fracUV95,$nCount50,$uCount50\n" 
    }
    close $fileH;         
    system("Rscript $param{vampPath}/bin/bruSeq/polyAFracUV_2.R $file $param{maxGeneDistance}");
    system("Rscript $param{vampPath}/bin/bruSeq/polyACounts_2.R $file $param{maxGeneDistance}"); 
}

sub createPolyAGeneTable {
    my ($chrom) = @_;
    my $tssGeneTable = 'tss_genes';                 
    my $nvrTable = "refGene_nvr_$param{refSeq}";  
    my $allTable = "refGene_all_$param{refSeq}";  
    my $nvrSQL = "SELECT * FROM $nvrTable WHERE chromosome = $chrom"; 
    my $allSQL = "SELECT * FROM $allTable WHERE chromosome = $chrom"; 
    my $polyASQL = "SELECT name2, strand, 
                          decode(strand,1,end_,corrstart_) polyA,
                          decode(strand,1,greatest(end_-$param{maxGeneDistance},corrstart_+$param{antisenseDistance}),
                                          end_-$param{maxGeneDistance}) start_,
                          decode(strand,1,end_+$param{maxGeneDistance},
                                          least(end_+$param{maxGeneDistance},corrstart_-$param{antisenseDistance})) end_,
                          decode(strand,1,1,-1) signMultiplier
                   FROM ($nvrSQL)";    
    my $overlapSQL = "SELECT p.name2
                      FROM ($polyASQL) p, ($allSQL) a
                      WHERE p.start_ <= a.end_
                        AND a.corrstart_ <= p.end_
                        AND p.name2 != a.name2";        
    my $geneSQL = "SELECT * FROM ($polyASQL) WHERE name2 NOT IN ($overlapSQL)";                                                  
    dropTable($tssGeneTable);
    runSQL("CREATE TABLE $tssGeneTable AS $geneSQL");
    return $tssGeneTable;
}
sub createPolyAJoinTable { 
    my ($tssHitsTable, $tssGeneTable) = @_; 
    my $tssJoinTable = "tss_join";              
    my $geneJoinSQL = "SELECT (round((h.position-g.polyA) / $param{binSize}) * $param{binSize})*signMultiplier bin,   
                              h.nCount, h.uCount, g.name2, g.strand
                       FROM $tssHitsTable h, $tssGeneTable g
                       WHERE h.position BETWEEN g.start_ AND g.end_ ";                                                                         
    dropTable($tssJoinTable);                 
    runSQL("CREATE TABLE $tssJoinTable AS $geneJoinSQL"); 
    return $tssJoinTable;     
}

sub createPolyADensityTable {
    my ($nonUVSample, $uvSample)  = @_;
    my $tssDensityTable = "tss_density";   
    my $unionSQL = "SELECT nvl(n.name2, u.name2) name2, n.density nDensity, u.Density uDensity 
                        FROM gmap_$nonUVSample\_UNQ n FULL OUTER JOIN gmap_$uvSample\_UNQ u
                          ON n.name2 = u.name2 ";
    dropTable($tssDensityTable);                 
    runSQL("CREATE TABLE $tssDensityTable AS $unionSQL"); 
    return $tssDensityTable;    
}

1;


