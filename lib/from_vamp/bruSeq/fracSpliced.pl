#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
defined $param{minND} or $param{minND} = 1;  #minimum normalized density value for allowed genes
defined $command{fracSpliced} or $command{fracSpliced} = ['singleThread', '24:00:00', 2000, 0];

sub fracSpliced {
    my ($earlySample, $lateSample) = @_;  
    foreach my $sample ($earlySample, $lateSample) {
        my $outFile = "$param{inputPath}/$sample.fracSpliced.csv";
        open my $outFileH, ">", $outFile or die "could not open $outFile: $!"; 
        print $outFileH "name2,eDensity,iDensity,splicedFrac,geneSize\n";
        
        foreach my $chrom(1..nChrom()){
        #foreach my $chrom(21..21){ 
        
            my ($geneTable, $exonTable, $intronTable) = createGeneTablesFracSpliced($earlySample, $chrom);
            my $exonDensTable = createFracSplicedTable('eDens', getMrnaDensitySQL($sample, $chrom, $exonTable)); 
            my $intronDensTable = createFracSplicedTable('iDens', getIntronDensitySQL($sample, $chrom, $intronTable)); 
            my $splicedFracTable = createSplicedFracTable($exonDensTable, $intronDensTable, $geneTable);
            runSQL("SELECT name2, eDensity, iDensity, splicedFrac, geneSize FROM $splicedFracTable", \my($name2, $eDensity, $iDensity, $splicedFrac, $geneSize));
            while (fetchRow()){ print $outFileH "$name2,$eDensity,$iDensity,$splicedFrac,$geneSize\n" }
        }
        close $outFileH;
        system("Rscript $param{vampPath}/bin/bruSeq/fracSplicedPlot.R $outFile");   
    } 
}
sub createGeneTablesFracSpliced {
    my ($earlySample, $chrom) = @_;
    my $geneSQL = "SELECT * FROM refGene_nvr_$param{refSeq} WHERE chromosome = $chrom";
    my $mapSQL =  "SELECT * FROM gmap_$earlySample\_unq     WHERE normalizeddensity >= $param{minND}";
    my $geneTable = createFracSplicedTable('genes', "SELECT g.name2, g.geneSize
                                                     FROM  ($geneSQL) g, ($mapSQL) m
                                                     WHERE g.name2 = m.name2"); 
    my $refGeneExonsUniqueTable = "refGeneExons_Unq_$param{refSeq}"; #single-entry list of EXONS
    my $exonTable = createFracSplicedTable('exons', "SELECT g.name2, $chrom chromosome, e.corrstart_, e.end_
                                                     FROM $refGeneExonsUniqueTable e, $geneTable g
                                                     WHERE e.name2 = g.name2"); 
    my $refGeneIntronsUniqueTable = "refGeneIntrons_Unq_$param{refSeq}"; #single-entry list of INTRONS
    my $intronTable = createFracSplicedTable('introns', "SELECT g.name2, $chrom chromosome, i.start_, i.corrend_
                                                       FROM $refGeneIntronsUniqueTable i, $geneTable g
                                                       WHERE i.name2 = g.name2");                                                                                                   
    return ($geneTable, $exonTable, $intronTable);                                                                                                   
}
sub createSplicedFracTable {
    my ($exonDensTable, $intronDensTable, $geneTable) = @_;
    my $fracSQL = "SELECT e.name2, e.density eDensity, i.density iDensity, i.density/e.density splicedFrac
                   FROM $exonDensTable e, $intronDensTable i
                   WHERE e.name2 = i.name2";
    my $joinSQL = "SELECT f.name2, f.eDensity, f.iDensity, f.splicedFrac, g.geneSize
                   FROM ($fracSQL) f,  $geneTable g
                   WHERE f.name2 = g.name2";
    return createFracSplicedTable('spliced', $joinSQL);
}
sub createFracSplicedTable {
    my ($tableSuffix, $sql) = @_;  
    status("    creating $tableSuffix table\n");
    my $outTable = "frSp_$tableSuffix";                                                                  
    dropTable($outTable);                 
    runSQL("CREATE TABLE $outTable AS $sql"); 
    return $outTable;      
}

1;

