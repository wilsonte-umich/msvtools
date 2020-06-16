#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{xx} or $command{xx} = ['multiThread', '24:00:00', 5000, 0];
defined $command{parseIntronHits} or $command{parseIntronHits} = ['multiThread', '12:00:00', 2000, 0];

sub parseIntronHits {
    #IMAP will restrict hits to gene strand if stranded is true
    my ($sample) = @_;
    my $iMapTable = newTable('IMap', $sample);
    calculateIntronDensity($sample, $iMapTable);
    calculateNormalizedIntronDensity($sample, $iMapTable);
    my $histogramTable = newTable('h', $iMapTable);
    createHitMapHistogram($histogramTable, $iMapTable);
}
sub calculateIntronDensity {
    my ($sample, $iMapTable) = @_;
    status("calculating hit density per intron...\n    Chr: ");
    foreach my $chrom (1..nChrom()){
        status("$chrom ");
        my $intronHitsSQL = getInIntronHitsSQL($sample, $chrom);
        my $densitySql = getPerFeatureDensitySQL($intronHitsSQL);
        my $iMapFile = "$iMapTable.csv";
        open my $iMapFileH, ">", $iMapFile or die "could not open $iMapFile\n";
        runSQL($densitySql, \my($name,$chrom,$start,$end,$coverage,$density));
        while(fetchRow()){ print $iMapFileH "$name,$chrom,$start,$end,$coverage,$density\n" }
        close $iMapFileH;
        loadData($iMapFile, $iMapTable, ",", "NAME2, CHROMOSOME, START_, CORREND_, COVERAGE, DENSITY"); 
    } 
    status("\n");  
}
sub calculateNormalizedIntronDensity {
    my ($sample, $iMapTable) = @_;
    status("calculating normalized hit density per intron...\n");
    my $genomeHitCount = getGenomeIntronHitCount($sample, $iMapTable);
    my $genomeSize = getGenomeIntronSize();
    my $genomeDensity = $genomeHitCount / $genomeSize;
    status("    $sample has $genomeDensity average intron hit density\n");
    runSQL("UPDATE $iMapTable SET NORMALIZEDDENSITY = DENSITY / $genomeDensity");
}
sub getGenomeIntronHitCount {
    my ($sample, $iMapTable) = @_;
    runSQL("SELECT sum(coverage) N from $iMapTable", \my$count);
    fetchRow();
    status("    $sample has $count intron hits\n");
    return $count;    
}
sub getGenomeIntronSize {
    my $refGeneIntronsUniqueTable = "refGeneIntrons_Unq_$param{refSeq}"; #single-entry list of INTRONS
    runSQL("SELECT sum(corrend_ - start_ - $param{readLength}) N from $refGeneIntronsUniqueTable", \my$size);
    fetchRow();
    status("    $param{refSeq} has $size bp in introns\n");
    return $size; 
}

1;


