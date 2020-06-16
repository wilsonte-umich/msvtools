#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{updateRefGene} or $command{updateRefGene} = ['singleThread', '4:00:00', 500, 0];

sub updateRefGene { #this sub refreshes refGene from UCSC and regenerates all derivative tables
    
    loadUCSCData("refGene");
    loadUCSCData("refGeneExons"); 
    loadUCSCData("chromInfo"); 
    loadUCSCData("cytoBand"); 
    loadUCSCData("rmsk"); 
    
    my $refGeneTable = "refGene_$param{refSeq}"; 
    my $refGeneExonsTable = "refGeneExons_$param{refSeq}"; 
    my $refGeneAllTable = getRefGeneTableName("All");
    my $refGeneUniqueTable = getRefGeneTableName("Unq");
    my $refGeneOverTable = getRefGeneTableName("Ovr");
    my $refGeneNonOverTable = getRefGeneTableName("Nvr");
    my $refGeneExonsUniqueTable = "refGeneExons_Unq_$param{refSeq}";
    my $refGeneIntronsUniqueTable = "refGeneIntrons_Unq_$param{refSeq}";
    
    updateRefGeneTable($refGeneAllTable,
                       "SELECT name2, chromosome, strand, 
                               min(corrstart_) corrstart_, max(end_) end_, 
                               max(end_) - min(corrstart_) + 1 geneSize, 0 mRnaSize
                        FROM $refGeneTable
                        GROUP BY name2, chromosome, strand"); #lose transcripts, just keep gene limits

    updateRefGeneTable($refGeneExonsUniqueTable,
                       "SELECT name2, chromosome, strand, corrstart_, end_
                          FROM $refGeneExonsTable
                          GROUP BY name2, chromosome, strand, corrstart_, end_"); #non-redundant gene exons          
    purgeOverlappingExons($refGeneExonsUniqueTable); #see comment below                  

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
                       
    parseRefGeneIntrons($refGeneExonsUniqueTable, $refGeneIntronsUniqueTable);
               
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

#POSSIBLE CONFOUNDER - to this point, some genes may have overlapping exons, which will each be listed
#i.e. some intron regions may sometimes not be spliced out of an exon on the list
#leading to segments countable as both exon and intron 
#this script purges overlapping exons, always labeling ambiguous regions as exons, i.e keeping the longest possible exon
sub purgeOverlappingExons {
    my ($refGeneExonsUniqueTable) = @_;
    my $outFile = "$refGeneExonsUniqueTable.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile\n";
    foreach my $chrom(1..nChrom()){
        foreach my $strand(1..2){
            runSQL("SELECT name2, corrstart_, end_
                    FROM $refGeneExonsUniqueTable
                    WHERE CHROMOSOME = $chrom
                      AND STRAND = $strand
                    ORDER BY name2, corrstart_, end_",
                    \my($name2, $corrStart, $end));
            my ($prevName2, $exonStart, $exonEnd) = ('');
            while (fetchRow()){
                if ($name2 eq $prevName2){
                   if($corrStart > $exonEnd){
                       print $outFileH "$prevName2,$chrom,$strand,$exonStart,$exonEnd\n";
                       $exonStart = $corrStart; 
                       $exonEnd = $end;
                   } else { #ovelrapping exons, take widest possible exon span
                       $exonEnd >= $end or $exonEnd = $end;
                   }
                } else {
                    defined $exonStart and print $outFileH "$prevName2,$chrom,$strand,$exonStart,$exonEnd\n";
                    $exonStart = $corrStart; 
                    $exonEnd = $end;
                }
                $prevName2 = $name2; 
            }
            print $outFileH "$prevName2,$chrom,$strand,$exonStart,$exonEnd\n";
        }        
    }
    close $outFileH;
    runSQL("DELETE FROM $refGeneExonsUniqueTable WHERE 1=1");
    loadData($outFile, $refGeneExonsUniqueTable, ",", "NAME2, CHROMOSOME, STRAND, CORRSTART_, END_");   
}

sub parseRefGeneIntrons { #create the unique intron table as the inverse of the unique, non-overlapping exon table
    my ($refGeneExonsUniqueTable, $refGeneIntronsUniqueTable) = @_;
    dropTable($refGeneIntronsUniqueTable); 
    runSQL("CREATE TABLE $refGeneIntronsUniqueTable 
             (NAME2 VARCHAR2(255), CHROMOSOME NUMBER, STRAND NUMBER, START_ NUMBER, CORREND_ NUMBER)");
    my $outFile = "$refGeneIntronsUniqueTable.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile\n";
    foreach my $chrom(1..nChrom()){
        foreach my $strand(1..2){
            runSQL("SELECT name2, corrstart_, end_
                    FROM $refGeneExonsUniqueTable
                    WHERE CHROMOSOME = $chrom
                      AND STRAND = $strand
                    ORDER BY name2, corrstart_, end_",
                    \my($name2, $corrStart, $end));
            my ($prevName2, $intronStart, $intronEnd) = ('');
            while (fetchRow()){
                if ($name2 eq $prevName2){
                   $intronEnd = $corrStart - 1;
                   $intronEnd > $intronStart and print $outFileH "$name2,$chrom,$strand,$intronStart,$intronEnd\n";
                }
                $intronStart = $end + 1; 
                $prevName2 = $name2;
            }
        }     
    }
    close $outFileH;
    loadData($outFile, $refGeneIntronsUniqueTable, ",", "NAME2, CHROMOSOME, STRAND, START_, CORREND_");
}


1;


