#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{addMrnaMaps} or $command{addMrnaMaps} = ['multiThread', '24:00:00', 5000, 0];

sub getInGeneHitsSQL { #list of hits in genes
    #a hit could end up on this list more than once, once in each overlapping gene
    my ($sample, $chrom, $geneType) = @_; #$geneType allows overrides
    my $hitsTable = getTableName('Hits', $sample);
    my $refGeneGenesTable = getRefGeneTableName($geneType); #undef geneType will reflex to param geneType
    my $inGeneSQL = "SELECT g.name2, g.corrstart_, g.end_, h.chromosome, h.position, h.strand, h.count_
                     FROM $hitsTable h, $refGeneGenesTable g
                     WHERE h.chromosome = $chrom
                       AND g.chromosome = $chrom
                       AND h.position >= g.corrstart_
                       AND h.position <= g.end_"; #DON'T group by, could purge desired duplicate reads in data
    return $inGeneSQL;
}
sub getNotInGeneHitsSQL { #list of hits outside of any annotated gene
    my ($sample, $chrom) = @_;
    my $hitsTable = getTableName('Hits', $sample);  
    my $inGeneSQL = getInGeneHitsSQL($sample, $chrom, "All");  
    my $hitKey = " chromosome || ':' || position || ':' || strand ";
    my $inGeneHitKeySQL = "SELECT $hitKey FROM ($inGeneSQL)";
    my $notInGeneSQL = "SELECT chromosome, position, strand, count_
                        FROM $hitsTable 
                        WHERE $hitKey NOT IN ($inGeneHitKeySQL)";
    return $notInGeneSQL;
}
sub getInMrnaHitsSQL { #list of hits in exons
    my ($sample, $chrom) = @_;
    my $hitsTable = getTableName('Hits', $sample);
    my $refGeneExonsUniqueTable = "refGeneExons_Unq_$param{refSeq}"; #single-entry list of EXONS
    my $inMrnaSQL = "SELECT g.name2, g.corrstart_, g.end_, h.chromosome, h.position, h.strand, h.count_
                     FROM $hitsTable h, $refGeneExonsUniqueTable g
                     WHERE h.chromosome = $chrom
                       AND g.chromosome = $chrom
                       AND h.position >= g.corrstart_
                       AND h.position <= g.end_"; #DON'T group by, could purge desired duplicate reads in data
    return $inMrnaSQL;
}
sub getGeneCountSQL { #count of hits in gene spans
    my ($sample, $chrom) = @_;
    my $hitsSQL = getInGeneHitsSQL($sample, $chrom);
    my $agg = "count(*)";
    $param{keepDups} and $agg = "sum(count_)";
    my $geneCountSQL = "SELECT name2, $agg nHits
                        FROM ($hitsSQL)
                        GROUP BY name2"; #list will NOT include genes with 0 hits
    return $geneCountSQL;   
}
sub getMrnaCountSQL { #count of hits in all gene exons
    my ($sample, $chrom) = @_;
    my $hitsSQL = getInMrnaHitsSQL($sample, $chrom);
    my $agg = "count(*)";
    $param{keepDups} and $agg = "sum(count_)";
    my $mRnaCountSQL = "SELECT name2, $agg nHits
                        FROM ($hitsSQL)
                        GROUP BY name2"; #list will NOT include genes with 0 hits         
    return $mRnaCountSQL;   
}
sub getGeneDensitySQL {
    my ($sample, $chrom) = @_; 
    my $geneCountSQL = getGeneCountSQL($sample, $chrom);
    my $refGeneGenesTable = getRefGeneTableName();
    my $geneDensitySQL = "SELECT g.name2, c.nHits COVERAGE, c.nHits / g.geneSize DENSITY
                          FROM ($geneCountSQL) c, 
                               (SELECT * FROM $refGeneGenesTable WHERE chromosome = $chrom) g
                          WHERE c.name2 = g.name2";
    return $geneDensitySQL;
}
sub getMrnaDensitySQL {
    my ($sample, $chrom) = @_; 
    my $mRnaCountSQL = getMrnaCountSQL($sample, $chrom);
    my $refGeneGenesTable = getRefGeneTableName();
    my $mRnaDensitySQL = "SELECT g.name2, c.nHits COVERAGE, c.nHits / g.mRnaSize DENSITY
                          FROM ($mRnaCountSQL) c, 
                               (SELECT * FROM $refGeneGenesTable WHERE chromosome = $chrom) g
                          WHERE c.name2 = g.name2";       
    return $mRnaDensitySQL;
}





#sub getDensityRankSQL {
#    my ($sample, $isMrna) = @_;
#    my $densitySql;
#    if ($isMrna) {
#        $densitySql = getMrnaDensitySQL($sample);
#    } else {
#        $densitySql = getGeneDensitySQL($sample);
#    }  
#    my $rankSQL = "SELECT sum(1) OVER (ORDER BY hitDensity ROWS UNBOUNDED PRECEDING) rank,
#                          sample, name2, chromosome, strand, corrstart_, end_, nHits, hitDensity
#                   FROM ($densitySql)
#                   GROUP BY sample, name2, chromosome, strand, corrstart_, end_, nHits, hitDensity";
#    return $rankSQL;
#}

##normalize to ??all hits in all genes??  ??median density??

#sub getCrosstabDensitySQL {
#    my ($sample1, $sample2, $isMrna1, $isMrna2, $scalar1, $scalar2) = @_;
#    my $densitySQL1 = getDensityRankSQL($sample1, $isMrna1);
#    my $densitySQL2 = getDensityRankSQL($sample2, $isMrna2);
#    my $crosstabSQL = "SELECT name2, chromosome, strand, corrstart_, end_,
#                              sum(decode(sample,'$sample1',nHits,0))      * $scalar1 $sample1\_nHits,
#                              sum(decode(sample,'$sample2',nHits,0))      * $scalar2 $sample2\_nHits,
#                              sum(decode(sample,'$sample1',hitDensity,0)) * $scalar1 $sample1\_hitDensity,
#                              sum(decode(sample,'$sample2',hitDensity,0)) * $scalar2 $sample2\_hitDensity,   
#                              sum(decode(sample,'$sample1',rank,0))                  $sample1\_rank,
#                              sum(decode(sample,'$sample2',rank,0))                  $sample2\_rank      
#                       FROM ( $densitySQL1 UNION ALL $densitySQL2 ) 
#                       GROUP BY name2, chromosome, strand, corrstart_, end_";
#    return $crosstabSQL;          
#}

1;
