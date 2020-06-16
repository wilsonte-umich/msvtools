#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{xx} or $command{xx} = ['multiThread', '24:00:00', 5000, 0];

sub getInGeneHitsSQL { #list of hits in genes
    my ($sample, $chrom, $geneType) = @_; #$geneType allows overrides
    my $refGeneGenesTable = getRefGeneTableName($geneType); #undef geneType will reflex to param geneType
    return getHitJoinSQL($sample, $chrom, $refGeneGenesTable, "name2", "chromosome", "corrstart_", "end_", "strand");
}
#In subs below, hit positions within readLength of the start of an exon
#are considered more likely to be part of the exon than the intron
#the position is actually in.  This results from the fact that mapping
#is done to the genome first, so many reads crossing splice junctions
#map, with mismatches, to the intron when in fact they should have mapped
#to the upstream exon.
sub getInMrnaHitsSQL { #list of hits in exons
    my ($sample, $chrom, $refGeneExonsUniqueTable) = @_;
    $refGeneExonsUniqueTable or $refGeneExonsUniqueTable = "refGeneExons_Unq_$param{refSeq}"; #single-entry list of EXONS
    return getHitJoinSQL($sample, $chrom, $refGeneExonsUniqueTable, "name2", "chromosome", "corrstart_ - $param{readLength}", "end_", "strand");
}
sub getInIntronHitsSQL { #list of hits in introns
    my ($sample, $chrom, $refGeneIntronsUniqueTable) = @_;
    $refGeneIntronsUniqueTable or $refGeneIntronsUniqueTable = "refGeneIntrons_Unq_$param{refSeq}"; #single-entry list of INTRONS
    return getHitJoinSQL($sample, $chrom, $refGeneIntronsUniqueTable, "name2", "chromosome", "start_", "corrend_ - $param{readLength}", "strand");
}
sub getHitJoinSQL {
    #a hit could end up on this list more than once, once in each overlapping feature
    my ($sample, $chrom, $featureTable, $featName, $featChrom, $featStart, $featEnd, $featStrand) = @_;
    my $hitsTable = getTableName('Hits', $sample);
    my $hitSQL = "SELECT * FROM $hitsTable WHERE chromosome = $chrom";
    my $featureSQL = "SELECT * FROM $featureTable WHERE $featChrom = $chrom";
    my $strandFilter = "";
    $param{stranded} and $featStrand and $strandFilter = " AND h.strand = f.$featStrand ";
    my $joinSQL = "SELECT f.$featName name, f.$featStart start_, f.$featEnd end_, h.chromosome, h.position, h.strand, h.count_
                   FROM ($hitSQL) h, ($featureSQL) f
                   WHERE h.position BETWEEN f.$featStart AND f.$featEnd $strandFilter"; #DON'T group by, could purge desired duplicate reads in data
    return $joinSQL;
}
sub getNotInGeneHitsSQL { #list of hits outside of any annotated gene
    my ($sample, $chrom) = @_;
    return getHitAntiJoinSQL($sample, $chrom, getInGeneHitsSQL($sample, $chrom, "All"));
}
sub getHitAntiJoinSQL {
    my ($sample, $chrom, $inFeatureSQL) = @_;
    my $hitsTable = getTableName('Hits', $sample);
    my $hitSQL = "SELECT * FROM $hitsTable WHERE chromosome = $chrom";
    my $hitKey = " chromosome || ':' || position || ':' || strand ";
    my $inFeatureHitKeySQL = "SELECT $hitKey FROM ($inFeatureSQL)";
    my $antiJoinSQL = "SELECT chromosome, position, strand, count_
                        FROM ($hitSQL)
                        WHERE $hitKey NOT IN ($inFeatureHitKeySQL)";
    return $antiJoinSQL;
}

#these subs return the count over all features contained within a gene
#e.g. all exons in the gene are summed
sub getGeneCountSQL { #count of hits in gene spans
    my ($sample, $chrom) = @_;
    return getFeatureCountSQL($sample, $chrom, getInGeneHitsSQL($sample, $chrom));
}
sub getMrnaCountSQL { #count of hits in all gene exons
    my ($sample, $chrom, $refGeneExonsUniqueTable) = @_;
    return getFeatureCountSQL($sample, $chrom, getInMrnaHitsSQL($sample, $chrom, $refGeneExonsUniqueTable));
}
sub getIntronCountSQL { #count of hits in all gene introns
    my ($sample, $chrom, $refGeneIntronsUniqueTable) = @_;
    return getFeatureCountSQL($sample, $chrom, getInIntronHitsSQL($sample, $chrom, $refGeneIntronsUniqueTable));
}
sub getFeatureCountSQL {
    my ($sample, $chrom, $hitJoinSQL) = @_;
    my $agg = "count(*)";
    $param{keepDups} and $agg = "sum(count_)";
    my $hitCountSQL = "SELECT name, $agg nHits FROM ($hitJoinSQL) GROUP BY name"; #list will NOT include features with 0 hits
    return $hitCountSQL;
}

#these subs return the density over all features contained within a gene
#e.g. all exons in the gene are summed
sub getGeneDensitySQL {
    my ($sample, $chrom) = @_; 
    return getPerGeneDensitySQL($sample, $chrom, getGeneCountSQL($sample, $chrom), "geneSize");
}
sub getMrnaDensitySQL {
    my ($sample, $chrom, $refGeneExonsUniqueTable) = @_; 
    return getPerGeneDensitySQL($sample, $chrom, getMrnaCountSQL($sample, $chrom, $refGeneExonsUniqueTable), "mRnaSize");
}
sub getIntronDensitySQL {
    my ($sample, $chrom, $refGeneIntronsUniqueTable) = @_; 
    return getPerGeneDensitySQL($sample, $chrom, getIntronCountSQL($sample, $chrom, $refGeneIntronsUniqueTable), "(genesize - mRnaSize)");
}
sub getPerGeneDensitySQL {
    my ($sample, $chrom, $featureCountSQL, $featureSize) = @_;
    my $refGeneGenesTable = getRefGeneTableName();
    my $refGeneSQL = "SELECT * FROM $refGeneGenesTable WHERE chromosome = $chrom";
    my $featureDensitySQL = "SELECT f.name, f.nHits COVERAGE, f.nHits / $featureSize DENSITY
                             FROM ($featureCountSQL) f, ($refGeneSQL) g
                             WHERE f.name = g.name2";
    return $featureDensitySQL;
}

#this sub returns the density over every individual feature
#e.g. each intron is a separate item on the list
sub getPerFeatureDensitySQL {
    my ($hitJoinSQL) = @_;
    my $agg = "count(*)";
    $param{keepDups} and $agg = "sum(count_)";
    my $featureDensitySQL = "SELECT name, chromosome, start_, end_, $agg COVERAGE, $agg / (end_ - start_ + 1) DENSITY
                             FROM ($hitJoinSQL)
                             GROUP BY name, chromosome, start_, end_";
    return $featureDensitySQL;
}

1;

