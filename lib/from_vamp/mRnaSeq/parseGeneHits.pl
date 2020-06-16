#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{xx} or $command{xx} = ['multiThread', '24:00:00', 5000, 0];
defined $param{geneType} or $param{geneType} = "Unq"; 
defined $param{isMrna} or $param{isMrna} = 0; 
defined $command{parseGeneHits} or $command{parseGeneHits} = ['multiThread', '48:00:00', 2000, 0];

sub parseGeneHits {
    #GMAP will restrict hits to gene strand if stranded is true
    my ($sample) = @_;
    getIsMrna(\my%isMrna);
    my $isMrna = $isMrna{$sample};
    my $gMapTable = newTable('GMap', "$sample\_$param{geneType}");  
    calculateGeneDensity($sample, $gMapTable, $isMrna);
    #my $gMapTable = getTableName('GMap', "$sample\_$param{geneType}");  
    calculateNormalizedGeneDensity($sample, $gMapTable, $isMrna);
    my $histogramTable = newTable('h', $gMapTable);
    createHitMapHistogram($histogramTable, $gMapTable);
}
sub getIsMrna {
    my ($isMrna) = @_;
    $param{isMrna} or return;
    my @isMrna = split(",", $param{isMrna});
    foreach my $sample(@isMrna){ $$isMrna{$sample} = 1 }
}
sub calculateGeneDensity {
    my ($sample, $gMapTable, $isMrna) = @_;
    status("calculating hit density per gene...\n    Chr: ");
    foreach my $chrom (1..nChrom()){
        status("$chrom ");
        my $densitySql;
        if ($isMrna) {
            $densitySql = getMrnaDensitySQL($sample, $chrom);
        } else {
            $densitySql = getGeneDensitySQL($sample, $chrom);
        } 
        my $gMapFile = "$gMapTable.csv";
        open my $gMapFileH, ">", $gMapFile;
        runSQL($densitySql, \my($name,$coverage,$density));
        while(fetchRow()){ print $gMapFileH "$name,$coverage,$density\n"}
        close $gMapFileH;
        loadData($gMapFile, $gMapTable, ",", "NAME2, COVERAGE, DENSITY"); 
    } 
    status("\n");  
}
sub calculateNormalizedGeneDensity {
    my ($sample, $gMapTable, $isMrna) = @_;
    status("calculating normalized hit density per gene...\n");
    my $genomeHitCount = getGenomeGeneHitCount($sample);
    my $genomeSize = getGenomeGeneSize($isMrna);
    my $genomeDensity = $genomeHitCount / $genomeSize;
    status("    $sample has $genomeDensity average gene hit density\n");
    runSQL("UPDATE $gMapTable SET NORMALIZEDDENSITY = DENSITY / $genomeDensity");
}
sub getGenomeGeneHitCount {
    my ($sample) = @_;  
    my $gmapTable = getTableName('GMap', "$sample\_$param{geneType}");
    runSQL("SELECT sum(coverage) N from $gmapTable", \my$count);
    fetchRow();
    status("    $sample has $count gene hits\n");
    return $count;    
}
sub getGenomeGeneSize {
    my ($isMrna) = @_;
    my $refGeneGenesTable = getRefGeneTableName();
    my $sizeField = "geneSize";
    $isMrna and $sizeField = "mRnaSize";
    runSQL("SELECT sum($sizeField) N from $refGeneGenesTable", \my$size);
    fetchRow();
    status("    $param{refSeq} has $size bp in genes\n");
    return $size; 
}

1;


