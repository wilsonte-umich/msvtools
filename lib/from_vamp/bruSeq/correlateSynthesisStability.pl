#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));
    
#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{xx} or $command{xx} = ['multiThread', '24:00:00', 5000, 0];
defined $command{correlateSynStab} or $command{correlateSynStab} = ['singleThread', '1:00:00', 2000, 0];

sub correlateSynStab {
    my ($sampleEarly, $sampleLate) = @_;
    my $synTable = getTableName('GMap', "$sampleEarly\_$param{geneType}");
    my $stabSql = getNormalizedStabilitySql($sampleEarly, $sampleLate);
    my @gmaps = (['synthesis', $synTable], ['stability', $stabSql]);
    my ($crosstabSql, $groupByName) = getSampleCrosstabSql("name2", "normalizedDensity", undef, undef, @gmaps);
    my ($ratioSql, $ratio1Name, $ratio2Name, $fractionName) = getTwoSampleRatioSql($crosstabSql, $groupByName, 'synthesis', 'stability');
    my $refGeneGenesTable = getRefGeneTableName();
    my $geneSQL = "SELECT g.chromosome, g.corrstart_, g.end_, g.strand, r.*
                   FROM ($ratioSql) r, $refGeneGenesTable g
                   WHERE r.name2 = g.name2";
    my $crosstabTable = "synStab_$sampleLate\_$sampleEarly";
    updateCrosstabTable($crosstabTable, $geneSQL);
    my $include = ['name2'];
    my $ndMax = 50;
    my @ndAxis = (1/$ndMax, $ndMax, 1);
    plotCorrelation($crosstabTable, undef, $include, undef,
                    'synthesis', undef, undef, @ndAxis,
                    'stability', undef, undef, @ndAxis,
                    1, [], []);
}

sub getNormalizedStabilitySql {
    my ($sampleEarly, $sampleLate) = @_;
    my $gmapTableLate = getTableName('GMap', "$sampleLate\_$param{geneType}");
    my $gmapTableEarly = getTableName('GMap', "$sampleEarly\_$param{geneType}");
    my @gmaps = ([$sampleLate, $gmapTableLate], [$sampleEarly, $gmapTableEarly]);
    my ($crosstabSql, $groupByName) = getSampleCrosstabSql("name2", "normalizedDensity", undef, undef, @gmaps);
    my ($ratioSql, $ratio1Name, $ratio2Name, $fractionName) = getTwoSampleRatioSql($crosstabSql, $groupByName, $sampleLate, $sampleEarly);
    return "SELECT name2, $ratio1Name normalizedDensity FROM ($ratioSql)";
}



1;
    
