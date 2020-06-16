#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
defined $param{asPadding} or $param{asPadding} = 1000;  #required gene free region at each end of gene
defined $param{asMinND} or $param{asMinND} = 0.01;  #minimum normalized density of genes considered (must have some basal expression...)
defined $param{asMinUVHits} or $param{asMinUVHits} = 100; #minimum UV hits in gene, as senseHits + antiSenseHits
defined $command{antiSense} or $command{antiSense} = ['singleThread', '48:00:00', 2000, 0];

my ($mergedSample);

sub antiSense {
    my ($nonUVSample, $uvSample) = @_;  
    $mergedSample = "$nonUVSample$uvSample";    
    $param{stranded} or die "param 'stranded' must be true for command 'antiSense'\n";
    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $uvHitCount) = getTSSScalars($nonUVSample, $uvSample); 
    status("  nonUV scalar = $nonUVScalar\n  UV scalar = $uvScalar\n");
    my $outFile = "antiSense_$mergedSample.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile $!"; 
    print $outFileH "name2,normD,nCount_sense,nCount_antisense,uCount_sense,uCount_antisense,uFracSense\n";

    foreach my $chrom(1..nChrom()){
    #foreach my $chrom(19..19){ 

        status("  chrom $chrom\n");
        my $hitsTable = createASHitsTable($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $chrom);
        my $geneTable = createASGeneTable($nonUVSample, $chrom);
        my $sumSQL = getASSumSQL($hitsTable, $geneTable);   
        runSQL($sumSQL, \my($name2,$normD,$nCount_sense,$nCount_antisense,$uCount_sense,$uCount_antisense));
        while (fetchRow()){ 
            my $uFracSense = "";
            my $denom = ($uCount_sense + $uCount_antisense);
            $denom >= $param{asMinUVHits} and $uFracSense = $uCount_sense / $denom;
            print $outFileH "$name2,$normD,$nCount_sense,$nCount_antisense,$uCount_sense,$uCount_antisense,$uFracSense\n" 
        }
        dropTable($hitsTable, $geneTable);    
    }
    close $outFileH;            
    status("  making plots\n"); 
    system("Rscript $param{vampPath}/bin/bruSeq/antiSensePlot.R $outFile");
}

sub createASHitsTable {
    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $chrom) = @_;
    my ($nAgg, $uAgg) = ("decode(n.count_,null,0,1)", "decode(u.count_,null,0,1)");
    $param{keepDups} and ($nAgg, $uAgg) = ("nvl(n.count_,0)", "nvl(u.count_,0)");
    my $nonUVSQL = "SELECT * FROM $nonUVHitsTable WHERE chromosome = $chrom";
    my $uvSQL =    "SELECT * FROM $uvHitsTable    WHERE chromosome = $chrom";
    my $hitJoinSQL =  "SELECT nvl(n.position,u.position) position, nvl(n.strand,u.strand) strand,
                              $nAgg * $nonUVScalar nCount, $uAgg * $uvScalar uCount
                       FROM ($nonUVSQL) n FULL OUTER JOIN ($uvSQL) u
                         ON (n.position = u.position AND n.strand = u.strand)";                    
    return createASTable('hits', $hitJoinSQL);
}
sub createASGeneTable {
    my ($nonUVSample, $chrom) = @_;  
    my $geneSQL = "SELECT name2, strand, corrstart_, end_
                   FROM refGene_unq_$param{refSeq}
                   WHERE chromosome = $chrom ";  
    my $ovrSQL = "SELECT g1.name2
                  FROM ($geneSQL) g1, ($geneSQL) g2
                  WHERE g1.corrstart_ - $param{asPadding} <= g2.end_
                    AND g2.corrstart_ <= g1.end_ + $param{asPadding}
                    AND g1.name2 != g2.name2";    
    my $gmapSQL = "SELECT * FROM gmap_$nonUVSample\_unq WHERE normalizeddensity >= $param{asMinND}";
    my $nvrSQL = "SELECT g.name2, g.strand, g.corrstart_, g.end_, m.normalizeddensity
                   FROM ($geneSQL) g, ($gmapSQL) m
                   WHERE g.name2 = m.name2
                     AND g.name2 NOT IN ($ovrSQL)";  
    return createASTable('genes', $nvrSQL);
}
sub getASSumSQL {
    my ($hitsTable, $geneTable) = @_;  
    my $strand = "case when g.strand = h.strand then 1 else -1 end";
    my $joinSQL = "SELECT g.name2, g.normalizeddensity, $strand strand, h.nCount, h.uCount    
                   FROM $hitsTable h, $geneTable g
                   WHERE h.position BETWEEN g.corrstart_ AND g.end_ ";  
    my $sumSQL = "SELECT name2, normalizeddensity,
                         sum(decode(strand, 1,nCount,0)) nCount_sense,
                         sum(decode(strand,-1,nCount,0)) nCount_antisense,
                         sum(decode(strand, 1,uCount,0)) uCount_sense,
                         sum(decode(strand,-1,uCount,0)) uCount_antisense
                  FROM ($joinSQL)
                  GROUP BY name2, normalizeddensity";       
    return $sumSQL;                                                                                      
}
sub createASTable {
    my ($tableSuffix, $sql) = @_;  
    status("    creating $tableSuffix table\n");
    my $outTable = "as$tableSuffix$mergedSample";                                                                  
    dropTable($outTable);                 
    runSQL("CREATE TABLE $outTable AS $sql"); 
    return $outTable;      
}

1;

