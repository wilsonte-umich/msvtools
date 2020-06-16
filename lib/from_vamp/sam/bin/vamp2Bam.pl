#!/usr/bin/perl
use strict;
use warnings;

my ($fragsTable, $normalOnly, $refSeq, $vampPath, $maxChrom, $dbLogin) = @ARGV;
our %refSeqs;
our %param = (dbLogin=>$dbLogin);
require "$vampPath/bin/sam/samSchema.pl";
require "$vampPath/bin/RefSeqs.pl";
require "$vampPath/bin/HandleDB_DBI.pl";
use vars(qw(%samFlags %samFields %reverseRefSeqs));

openOracle();
print  "\@HD\tVN:1.4\tSO:unsorted\n";    
foreach my $chrom(1..$maxChrom){
     my $sq = $chrom;
     fixBamSexChrom(\$sq);
     runSQL("SELECT end_ FROM chrominfo_$refSeq WHERE chromosome = $chrom", \my$chromEnd);
     fetchRow();
     print "\@SQ\tSN:$sq\tLN:$chromEnd\n"; 
}
my $readGroupSQL = "SELECT (PAIRID/1E8 - trunc(PAIRID/1E8))*1E8 readGroup FROM $fragsTable";
$readGroupSQL = "SELECT readGroup FROM ($readGroupSQL) GROUP BY readGroup";
runSQL($readGroupSQL, \my($rg));
while(fetchRow()){ print "\@RG\tID:$rg\n" }    
my $anomalyFilter = " OR (FRAGMENTTYPE > 0 AND FRAGMENTTYPE <= 16 AND NSETSPAIR = 1 AND NSETSFRAG = 1) ";
$normalOnly and $anomalyFilter = ""; 
my $fragsSQL = "SELECT trunc(PAIRID/1E8) QNAME, (PAIRID/1E8 - trunc(PAIRID/1E8))*1E8 readGroup,
                CHROMOSOME1 RNAME1, decode(CHROMOSOME2,0,CHROMOSOME1,CHROMOSOME2) RNEXT1, decode(CHROMOSOME2,0,CHROMOSOME1,CHROMOSOME2) RNAME2, CHROMOSOME1 RNEXT2, 
                POSITION1 POS1, POSITION2 PNEXT1, POSITION2 POS2, POSITION1 PNEXT2,
                LENGTH1 || 'M' CIGAR1, LENGTH2 || 'M' CIGAR2,
                decode(STRAND1,1,0,$samFlags{isReverse}) isReverse1, decode(STRAND2,1,0,$samFlags{nextIsReverse}) nextIsReverse1,
                decode(STRAND2,1,0,$samFlags{isReverse}) isReverse2, decode(STRAND1,1,0,$samFlags{nextIsReverse}) nextIsReverse2
                FROM $fragsTable
                WHERE (FRAGMENTTYPE = 0 $anomalyFilter)
                  AND CHROMOSOME1 <= $maxChrom AND CHROMOSOME2 <= $maxChrom";                           
runSQL($fragsSQL, \my($qname, $readGroup,
                      $rName1, $rNext1, $rName2, $rNext2,
                      $pos1, $pNext1, $pos2, $pNext2,
                      $cigar1, $cigar2,
                      $isReverse1, $nextIsReverse1, $isReverse2, $nextIsReverse2));     
#my $counter = 0;        
while (fetchRow()){
    my $tLen1 = fixBamRNames(\$rName1,\$rNext1,$pos1,$pNext1);
    my $tLen2 = fixBamRNames(\$rName2,\$rNext2,$pos2,$pNext2);
    my $flag1 = $samFlags{isMultiFrag} + $samFlags{isFirst} + $isReverse1 + $nextIsReverse1;
    my $flag2 = $samFlags{isMultiFrag} + $samFlags{isLast}  + $isReverse2 + $nextIsReverse2;
    print join("\t", $qname,$flag1,$rName1,$pos1,255,$cigar1,$rNext1,$pNext1,$tLen1,'*','*','RG:Z:'.$readGroup)."\n";
    print join("\t", $qname,$flag2,$rName2,$pos2,255,$cigar2,$rNext2,$pNext2,$tLen2,'*','*','RG:Z:'.$readGroup)."\n";
    #$counter++;
    #$counter == 10000 and last;
}
closeDB();


sub fixBamRNames {
    my($rName,$rNext,$pos,$pNext) = @_;
    fixBamSexChrom($rName);
    fixBamSexChrom($rNext);
    if($$rName eq $$rNext){
        $$rNext = '*';
        return $pNext - $pos;
    } else  {
        return 0;
    }
}
sub fixBamSexChrom {
    my ($rName) = @_;
    $reverseRefSeqs{$refSeq}{$$rName} =~ m/X/ and $$rName = 'X' and return;
    $reverseRefSeqs{$refSeq}{$$rName} =~ m/Y/ and $$rName = 'Y';
}

1;

