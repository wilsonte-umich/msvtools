#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#binSize is required but already defined by VAMP
defined $param{minBinHits} or $param{minBinHits} = 5; #3 hits/bin = 95% chance that a bin will have >= 1 hit
defined $param{maxGeneDistance} or $param{maxGeneDistance} = 100000; 
defined $param{antisenseDistance} or $param{antisenseDistance} = 20000; 
defined $command{tssMedian} or $command{tssMedian} = ['singleThread', '48:00:00', 2000, 0];

sub tssMedian {
    my ($nonUVSample, $uvSample) = @_;  
    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $uvHitCount) = getTSSScalars($nonUVSample, $uvSample); 
    status("  nonUV scalar = $nonUVScalar\n  UV scalar = $uvScalar\n");
    my $fracUVTable = 'tss_fracUV';
    foreach my $senseType('sense','antisense'){
        status("$senseType\n");
        dropTable($fracUVTable);
        runSQL("CREATE TABLE $fracUVTable (name2 varchar2(255), strand number, bin number, fracUV number)");
        my $fracUVFile = "$fracUVTable.csv";
        open my $fracUVFileH, ">", $fracUVFile or die "could not open $fracUVFile: $!"; 
        foreach my $chrom(1..nChrom()){
            status("  chrom $chrom\n");
            my $tssHitsTable = createTssHitsTable($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $chrom);
            my $tssGeneTable = createTssGeneTable($chrom);
            my $tssJoinTable = createTssJoinTable($tssHitsTable, $tssGeneTable, $senseType);                                               
            my $binSumSQL = "SELECT name2, strand, bin, sum(uCount) uCount, sum(nCount) nCount
                             FROM $tssJoinTable
                             GROUP BY name2, strand, bin";
            my $fracUVSQL = "SELECT name2, strand, bin, uCount/(uCount + nCount) fracUV
                             FROM ($binSumSQL)
                             WHERE (uCount > $param{minBinHits} OR nCount > $param{minBinHits}) ";  
            runSQL($fracUVSQL, \my($name2, $strand, $bin, $fracUV));
            while (fetchRow()){ print $fracUVFileH "$name2,$strand,$bin,$fracUV\n" }
        }
        close $fracUVFileH;
        loadData($fracUVFile, $fracUVTable, ",", "NAME2, STRAND, BIN, FRACUV");               
        # my $medianSQL = "SELECT strand, bin, count(*) nGenes,
                         # percentile_cont(0.05) WITHIN GROUP (ORDER BY fracUV) fracUV5, 
                         # percentile_cont(0.5) WITHIN GROUP (ORDER BY fracUV) fracUV50,
                         # percentile_cont(0.95) WITHIN GROUP (ORDER BY fracUV) fracUV95  
                         # FROM $fracUVTable
                         # GROUP BY strand, bin";            
        # my $seriesSQL = "SELECT bin,
                                # sum(decode(strand,1,nGenes,'')) forwardNGenes,
                                # sum(decode(strand,2,nGenes,'')) reverseNGenes,
                                # sum(decode(strand,1,fracUV5,'')) forward5,
                                # sum(decode(strand,2,fracUV5,'')) reverse5,
                                # sum(decode(strand,1,fracUV50,'')) forward50,
                                # sum(decode(strand,2,fracUV50,'')) reverse50,
                                # sum(decode(strand,1,fracUV95,'')) forward95,
                                # sum(decode(strand,2,fracUV95,'')) reverse95
                         # FROM ($medianSQL)  
                         # GROUP BY bin
                         # ORDER BY bin";      
        # runSQL($seriesSQL, \my($bin,$forwardNGenes,$reverseNGenes,$forward5,$reverse5,$forward50,$reverse50,$forward95,$reverse95));
        # my $file = "$param{inputPath}/tssMedian.output.$senseType.csv";
        # open my $fileH, ">", $file or die "could not open $file: $!\n";
        # print $fileH "bin,forwardNGenes,reverseNGenes,forward5,reverse5,forward50,reverse50,forward95,reverse95\n";
        # while(fetchRow()){ 
            # $forwardNGenes or $forwardNGenes = '';
            # $reverseNGenes or $reverseNGenes = '';
            # $forward5 or $forward5 = '';
            # $reverse5 or $reverse5 = '';
            # $forward50 or $forward50 = '';
            # $reverse50 or $reverse50 = '';  
            # $forward95 or $forward95 = '';
            # $reverse95 or $reverse95 = '';  
            # print $fileH "$bin,$forwardNGenes,$reverseNGenes,$forward5,$reverse5,$forward50,$reverse50,$forward95,$reverse95\n" 
        # }
        # close $fileH;         
        # system("Rscript $param{vampPath}/bin/bruSeq/tssMedian.R $file $param{maxGeneDistance}"); 
        my $medianSQL = "SELECT bin, count(*) nGenes,
                         percentile_cont(0.05) WITHIN GROUP (ORDER BY fracUV) fracUV5, 
                         percentile_cont(0.5) WITHIN GROUP (ORDER BY fracUV) fracUV50,
                         percentile_cont(0.95) WITHIN GROUP (ORDER BY fracUV) fracUV95  
                         FROM $fracUVTable
                         GROUP BY bin";            
        runSQL($medianSQL, \my($bin,$nGenes,$fracUV5,$fracUV50,$fracUV95));
        my $file = "$param{inputPath}/tssMedian.output.$senseType.csv";
        open my $fileH, ">", $file or die "could not open $file: $!\n";
        print $fileH "bin,nGenes,fracUV5,fracUV50,fracUV95\n";
        while(fetchRow()){ 
            print $fileH "$bin,$nGenes,$fracUV5,$fracUV50,$fracUV95\n" 
        }
        close $fileH;         
        system("Rscript $param{vampPath}/bin/bruSeq/tssMedian2.R $file $param{maxGeneDistance}");
    }
}

sub createTssHitsTable {
    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $chrom) = @_;
    my $tssHitsTable = 'tss_hits';
    my ($nAgg, $uAgg) = ("decode(n.count_,null,0,1)", "decode(u.count_,null,0,1)");
    $param{keepDups} and ($nAgg, $uAgg) = ("nvl(n.count_,0)", "nvl(u.count_,0)");
    my $nonUVSQL = "SELECT * FROM $nonUVHitsTable WHERE chromosome = $chrom";
    my $uvSQL = "SELECT * FROM $uvHitsTable WHERE chromosome = $chrom";
    my $hitJoinSQL = "SELECT nvl(n.position,u.position) position,
                             $nAgg * $nonUVScalar nCount, $uAgg * $uvScalar uCount
                      FROM ($nonUVSQL) n FULL OUTER JOIN ($uvSQL) u
                        ON (n.position = u.position)";                                  
    dropTable($tssHitsTable);
    runSQL("CREATE TABLE $tssHitsTable AS $hitJoinSQL");
    return $tssHitsTable;
}
sub createTssGeneTable {
    my ($chrom) = @_;
    my $tssGeneTable = 'tss_genes';                       
    my $geneTable = "refGene_nvr_$param{refSeq}";     
    my $chromSQL = "SELECT * FROM $geneTable WHERE chromosome = $chrom"; 
    my $geneSQL = "SELECT name2, strand, 
                          decode(strand,1,corrstart_,end_) senseTSS,  
                          decode(strand,1,end_,corrstart_) antisenseTSS, 
                          decode(strand,1,corrstart_,greatest(corrstart_+$param{antisenseDistance}, end_-$param{maxGeneDistance})) senseStart, 
                          decode(strand,1,least(end_-$param{antisenseDistance},corrstart_+$param{maxGeneDistance}),end_) senseEnd,
                          decode(strand,1,greatest(corrstart_+$param{antisenseDistance}, end_-$param{maxGeneDistance}),corrstart_) antisenseStart, 
                          decode(strand,1,end_,least(end_-$param{antisenseDistance},corrstart_+$param{maxGeneDistance})) antisenseEnd
                   FROM ($chromSQL)";                        
    dropTable($tssGeneTable);
    runSQL("CREATE TABLE $tssGeneTable AS $geneSQL");
    return $tssGeneTable;
}
# sub createTssJoinTable {
    # my ($tssHitsTable, $tssGeneTable, $senseType) = @_; 
    # my $tssJoinTable = "tss_join";       
    # my $geneJoinSQL = "SELECT round((h.position-g.$senseType"."TSS) / $param{binSize}) * $param{binSize} bin,   
                              # h.nCount, h.uCount, g.name2, g.strand
                       # FROM $tssHitsTable h, $tssGeneTable g
                       # WHERE h.position BETWEEN g.$senseType"."Start AND g.$senseType"."End";                              
    # dropTable($tssJoinTable);                 
    # runSQL("CREATE TABLE $tssJoinTable AS $geneJoinSQL"); 
    # return $tssJoinTable;     
# }
sub createTssJoinTable {
    my ($tssHitsTable, $tssGeneTable, $senseType) = @_; 
    my $tssJoinTable = "tss_join";   
    my $geneJoinSQL = "SELECT abs(round((h.position-g.$senseType"."TSS) / $param{binSize}) * $param{binSize}) bin,   
                              h.nCount, h.uCount, g.name2, g.strand
                       FROM $tssHitsTable h, $tssGeneTable g
                       WHERE h.position BETWEEN g.$senseType"."Start AND g.$senseType"."End";                              
    dropTable($tssJoinTable);                 
    runSQL("CREATE TABLE $tssJoinTable AS $geneJoinSQL"); 
    return $tssJoinTable;     
}

1;


