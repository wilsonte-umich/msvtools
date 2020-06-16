#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{updateRefGene} or $command{updateRefGene} = ['singleThread', '2:00:00', 500, 0];

sub updateRefGene { #this sub refreshes refGene from UCSC and regenerates all derivative tables
    
    loadUCSCData("refGene");
    loadUCSCData("refGeneExons"); 
    
    my $refGeneTable = "refGene_$param{refSeq}"; 
    my $refGeneExonsTable = "refGeneExons_$param{refSeq}"; 
    my $refGeneAllTable = getRefGeneTableName("All");
    my $refGeneUniqueTable = getRefGeneTableName("Unq");
    my $refGeneOverTable = getRefGeneTableName("Ovr");
    my $refGeneNonOverTable = getRefGeneTableName("Nvr");
    my $refGeneExonsUniqueTable = "refGeneExons_Unq_$param{refSeq}";
    
    updateRefGeneTable($refGeneAllTable,
                       "SELECT name2, chromosome, strand, 
                               min(corrstart_) corrstart_, max(end_) end_, 
                               max(end_) - min(corrstart_) + 1 geneSize, 0 mRnaSize
                        FROM $refGeneTable
                        GROUP BY name2, chromosome, strand"); #lose transcripts, just keep gene limits
                        
    #POSSIBLE CONFOUNDER - some genes may have overlapping exons, which will each be listed!!
    #i.e. some intron regions may sometimes not be spliced out of an exon on the list
    #leading to segment countable as both exon and intron 
    updateRefGeneTable($refGeneExonsUniqueTable,
                       "SELECT name2, chromosome, strand, corrstart_, end_
                          FROM $refGeneExonsTable
                          GROUP BY name2, chromosome, strand, corrstart_, end_"); #non-redundant gene exons                            

    runSQL("UPDATE $refGeneAllTable g SET mRnaSize = 
            (SELECT sum(e.end_ - e.corrstart_ + 1) mRnaSize 
             FROM $refGeneExonsUniqueTable e
             WHERE g.name2 = e.name2
             GROUP BY e.name2)"); #fill mRnaSize into genes table - as always this is slow, but only takes ~5 min
         
    my $countSQL = "SELECT name2, count(*) N FROM $refGeneAllTable GROUP BY name2";
    my $singleSQL = "SELECT name2 FROM ($countSQL) WHERE N = 1";   
    updateRefGeneTable($refGeneUniqueTable,
                       "SELECT * FROM $refGeneAllTable WHERE name2 IN ($singleSQL)"); #lose ambiguously named genes                   
    updateRefGeneTable($refGeneOverTable,
                       "SELECT g1.*
                        FROM $refGeneUniqueTable g1, $refGeneUniqueTable g2
                        WHERE g1.chromosome = g2.chromosome
                        AND g1.corrstart_ <= g2.end_
                        AND g2.corrstart_ <= g1.end_
                        AND g1.name2 != g2.name2"); #overlapping unique genes
    my $ovName2SQL = "SELECT name2 FROM $refGeneOverTable";                    
    updateRefGeneTable($refGeneNonOverTable,
                       "SELECT * FROM $refGeneUniqueTable WHERE name2 NOT IN ($ovName2SQL)"); #non-overlapping unique genes                 
    
                
}

sub getRefGeneTableName {
    my ($geneType) = @_;
    $geneType or $geneType = $param{geneType};
    return "refGene_$geneType\_$param{refSeq}"; 
}

sub updateRefGeneTable {
    my ($tableName, $sql) = @_;   
    status("creating derivative table $tableName\n"); 
    dropTable($tableName);
    runSQL("CREATE TABLE $tableName AS $sql");
}


1;
