#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#binSize is required but already defined by VAMP
defined $param{minTssND} or $param{minTssND} = 1;
defined $command{corrTssGMap} or $command{corrTssGMap} = ['singleThread', '48:00:00', 2000, 0];

sub corrTssGMap {
    my ($nonUVsample1, $uvSample1, $nonUVsample2, $uvSample2) = @_;
    my $refGeneTable = getRefGeneTableName();
    my $tssTable = getTableName('TSS', "$uvSample1\_$uvSample2");
    unless(tableExists($tssTable)){
        ($uvSample1, $uvSample2) = ($uvSample2, $uvSample1);
        $tssTable = getTableName('TSS', "$uvSample1\_$uvSample2");
        tableExists($tssTable) or die "could not find TSS table for $uvSample1, $uvSample2";
    }
    my $ctTable = "gmap_$nonUVsample1\_$nonUVsample2\_$param{geneType}_ct";
    unless(tableExists($ctTable)){
        ($nonUVsample1, $nonUVsample2) = ($nonUVsample2, $nonUVsample1);
        $ctTable = "gmap_$nonUVsample1\_$nonUVsample2\_$param{geneType}_ct";
        tableExists($ctTable) or die "could not find TSS table for $nonUVsample1, $nonUVsample2";
    }
    my $enhTable = getTableName('Enh', "$uvSample1\_$uvSample2");

    my $etssSql = "SELECT chromosome, round(start_+(end_-start_)/2) tss,
                          $uvSample1, fraction_$uvSample1
                   FROM $tssTable
                   WHERE ingene = 0
                     AND greatest($uvSample1, $uvSample2) > $param{minTssND}";
    my $gtssSql = "SELECT g.name2, g.chromosome, decode(g.strand,1,g.corrstart_,g.end_) tss,
                          ct.$nonUVsample1, ct.fraction_$nonUVsample1
                   FROM $ctTable ct, $refGeneTable g
                   WHERE ct.name2 = g.name2";
    my $distSql = "SELECT e.chromosome, e.tss etss, e.$uvSample1, e.fraction_$uvSample1,
                           g.name2, g.tss gtss, g.$nonUVsample1, g.fraction_$nonUVsample1,
                           abs(e.tss - g.tss) dist
                   FROM ($etssSql) e, ($gtssSql) g
                   WHERE e.chromosome = g.chromosome";
    my $distTable = substr("tssgmap$nonUVsample1$uvSample1$nonUVsample2$uvSample2", 0, 30);
    dropTable($distTable);
    runSQL("CREATE TABLE $distTable AS $distSql");
    my $minDistSql = "SELECT chromosome, etss, min(dist) minDist
                      FROM $distTable
                      GROUP BY chromosome, etss";
    my $enhSql = "SELECT d.*
                   FROM $distTable d, ($minDistSql) md
                   WHERE d.chromosome = md.chromosome
                     AND d.etss = md.etss
                     AND d.dist = md.minDist";
                     
    dropTable($enhTable);
    runSQL("CREATE TABLE $enhTable AS $enhSql");
    
    dropTable($distTable);


    my @fractionAxis = (0, 1, undef);
    my $fractionBinSize = 0.05;
    plotCorrelation($enhTable, undef, [], $enhTable,
                "fraction_$uvSample1", undef, "TSS fraction_$uvSample1", @fractionAxis,
                "fraction_$nonUVsample1", undef, "nearest gene fraction_$nonUVsample1", @fractionAxis,
                1, [], []);


    plotCorrelation($enhTable, undef, [], $enhTable,
                $nonUVsample1, undef, undef, undef, undef, undef,
                $uvSample1, undef, undef, undef, undef, undef,
                1, [], []);
}


1;

