#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{xx} or $command{xx} = ['multiThread', '24:00:00', 5000, 0];
defined $command{parseExonHits} or $command{parseExonHits} = ['multiThread', '48:00:00', 2000, 0];
defined $command{parseIntronHits} or $command{parseIntronHits} = ['multiThread', '48:00:00', 2000, 0];

sub parseExonHits {
    my ($sample) = @_;
    parseFeatureHits($sample, 'EMap', 'exon', \&getInMrnaHitsSQL, \&getGenomeExonSize);
}
sub parseIntronHits {
    my ($sample) = @_;
    parseFeatureHits($sample, 'IMap', 'intron', \&getInIntronHitsSQL, \&getGenomeIntronSize);
}
sub parseFeatureHits {
    #xMAP will restrict hits to gene strand if stranded is true
    my ($sample, $mapType, $featureName, $hitsSqlSub, $getSizeSub) = @_;
    my $mapTable = newTable($mapType, $sample);
    calculateFeatureDensity($sample, $mapTable, $featureName, $hitsSqlSub);
    calculateNormalizedFeatureDensity($sample, $mapTable, $featureName, $getSizeSub); 
    my $histogramTable = newTable('h', $mapTable);
    createHitMapHistogram($histogramTable, $mapTable);
}
sub calculateFeatureDensity {
    my ($sample, $mapTable, $featureName, $hitsSqlSub) = @_;
    status("calculating hit density per $featureName...\n    Chr: ");
    foreach my $chrom (1..nChrom()){
        status("$chrom ");
        my $hitsSQL = &$hitsSqlSub($sample, $chrom);
        my $densitySql = getPerFeatureDensitySQL($hitsSQL);
        my $mapFile = "$mapTable.csv";
        open my $mapFileH, ">", $mapFile or die "could not open $mapFile\n";
        runSQL($densitySql, \my($name,$chrom,$start,$end,$coverage,$density));
        while(fetchRow()){ print $mapFileH "$name,$chrom,$start,$end,$coverage,$density\n" }
        close $mapFileH;
        loadData($mapFile, $mapTable, ",", "NAME2, CHROMOSOME, START_, END_, COVERAGE, DENSITY"); 
    } 
    status("\n");  
}
sub calculateNormalizedFeatureDensity {
    my ($sample, $mapTable, $featureName, $getSizeSub) = @_;
    status("calculating normalized hit density per $featureName...\n");
    my $genomeHitCount = getGenomeFeatureHitCount($sample, $mapTable, $featureName);
    my $genomeSize = &$getSizeSub();
    my $genomeDensity = $genomeHitCount / $genomeSize;
    status("    $sample has $genomeDensity average $featureName hit density\n");
    runSQL("UPDATE $mapTable SET NORMALIZEDDENSITY = DENSITY / $genomeDensity");
}
sub getGenomeFeatureHitCount {
    my ($sample, $mapTable, $featureName) = @_;
    runSQL("SELECT sum(coverage) N from $mapTable", \my$count);
    fetchRow();
    status("    $sample has $count $featureName hits\n");
    return $count;    
}
sub getGenomeExonSize {
    my $refGeneTable = "refGeneExons_Unq_$param{refSeq}"; #single-entry list of EXONS
    runSQL("SELECT sum(end_ - (corrstart_ - $param{readLength}) + 1) N from $refGeneTable", \my$size);
    fetchRow();
    status("    $param{refSeq} has $size bp in exons\n");
    return $size; 
}
sub getGenomeIntronSize {
    my $refGeneTable = "refGeneIntrons_Unq_$param{refSeq}"; #single-entry list of INTRONS
    runSQL("SELECT sum((corrend_ - $param{readLength}) - start_ + 1) N from $refGeneTable", \my$size);
    fetchRow();
    status("    $param{refSeq} has $size bp in introns\n");
    return $size; 
}

1;

