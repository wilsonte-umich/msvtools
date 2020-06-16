#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#collateDiscrepancies.pl pulls catenated discrepancies out of Frags tables
#and splits them into single position entries in Discs table
#only accepting discrepancies from Normal Fragments, SingleReads and Fragments uniquely placed into Sets
#Counts are maintained for identical discrepancies, matching ambiguous reads,
#and the total number of crossing reads.
#Functional consequence is determined using gene tables if available.
#When comparing two samples, LOD scores are calculated.
#########################################################################

use vars(qw(%param %types %fields %refSeqs %geneticCode));

my %allowedTypes;
my %mismatchTypes;

sub processDiscrepancies1{
    my ($sample) = @_;
    initializeTypes();
    my %discCounts; #number of times a dicrepancy was detected
    my %consistentCounts; #number of times a discrepancy OR matching ambiguous were detected
    my %readCounts; #number of reads crossing a discrepancy position (mutant or not)    
    #process sample and collect discrepancy counts
    fillDiscCounts($sample, \%discCounts);
    #calculate consistent counts
    status("  calculating consistent counts...\n");
    foreach my $disc (keys %discCounts) {
        $consistentCounts{$disc} = $discCounts{$disc}; #initialize sample consistentCount to discCount
        my ($chrom, $pos, $discType, $extra) = split(":", $disc);  #retrieve discrepancy information
        if (defined $mismatchTypes{$discType}){ #include ambiguous positions in consistent dicrepancy count for mismatches   
                my $ambiguous = getDisc($chrom, $pos, $types{Discs}{N}, $extra);
                $consistentCounts{$disc} += ($discCounts{$ambiguous} or 0);    
        }
        my $discPos = "$chrom:$pos";  #initialize readCounts to 0 for discrepancy positions
        $readCounts{$discPos} = 0;   
    }         
    #scan for number of reads crossing discrepancy positions
    fillReadCounts($sample, \%readCounts);   
    #commit discrepancies
    status("  committing discrepancies...\n");
    my %chromSeqs;
    loadChromosomes(\%chromSeqs);
    my $cdsTable = "\UCDS_$param{refSeq}";
    tableExists($cdsTable) or $cdsTable = 0;
    my $nDisc;
    my $timeStamp = getTimeStamp();
    my $discsTable = newTable('Discs', $sample);
    my $discsFile = "$discsTable.csv";
    open my $discsFileH, ">", $discsFile;
    foreach my $disc(keys %discCounts){
        my ($chrom, $pos, $discType, $extra) = split(":", $disc);  #retrieve discrepancy information
        if ($consistentCounts{$disc} > 1){        
            $nDisc++;
            my $discID = ($nDisc * 1E8) + $timeStamp;
            my $discPos = "$chrom:$pos";
            my $s1_sameCount = $discCounts{$disc};
            my $s2_sameCount = 0;
            my $s1_consistentCount = $consistentCounts{$disc};
            $s1_consistentCount <= $readCounts{$discPos} or $s1_consistentCount = $readCounts{$discPos};  #fixes rare error in consistent count (I think propagated from Pass)
            my $s2_consistentCount = 0;
            my $sameLOD = 0;
            my $consistentLOD = 0;
            my $consequence = getConsequence($chrom, $pos, $discType, \%chromSeqs, $cdsTable);
            print $discsFileH join(",", ($discID,
                                         $chrom, $pos, $discType, $extra,
                                         $s1_sameCount, $s1_consistentCount, $readCounts{$discPos},
                                         $s2_sameCount, $s2_consistentCount, 0,
                                         $sameLOD, $consistentLOD, $consequence))."\n"; 
        }  
    }
    close $discsFileH;
    loadData($discsFile, $discsTable, ",", "DISCREPANCYID, 
                                            CHROMOSOME, POSITION, DISCREPANCYTYPE, EXTRA,
                                            S1_SAMECOUNT, S1_CONSISTENTCOUNT, S1_READCOUNT, 
                                            S2_SAMECOUNT, S2_CONSISTENTCOUNT, S2_READCOUNT,
                                            SAMELOD, CONSISTENTLOD, CONSEQUENCE");
    purgeAmbiguous($discsTable);   
}

sub processDiscrepancies2{
    my ($sample1, $sample2) = @_;
    initializeTypes();
    my %discCounts = (1 => {}, 2 => {}); #number of times a dicrepancy was detected
    my %consistentCounts = (1 => {}, 2 => {}); #number of times a discrepancy OR matching ambiguous were detected
    my %readCounts = (1 => {}, 2 => {}); #number of reads crossing a discrepancy position (mutant or not)    
    #process each sample and collect discrepancy counts
    fillDiscCounts($sample1, $discCounts{1});
    fillDiscCounts($sample2, $discCounts{2});    
    #collect combined discrepancies and calculate consistent counts
    status("  calculating consistent counts...\n");
    my %discs;
    foreach my $sampleN(1..2){
        foreach my $disc (keys %{$discCounts{$sampleN}}) {
            $discs{$disc} = 0;  #record the discrepancy
            ${$consistentCounts{$sampleN}}{$disc} = ${$discCounts{$sampleN}}{$disc}; #initialize sample consistentCount to discCount
            my ($chrom, $pos, $discType, $extra) = split(":", $disc);  #retrieve discrepancy information
            if (defined $mismatchTypes{$discType}){ #include ambiguous positions in consistent dicrepancy count for mismatches   
                    my $ambiguous = getDisc($chrom, $pos, $types{Discs}{N}, $extra);
                    ${$consistentCounts{$sampleN}}{$disc} += (${$discCounts{$sampleN}}{$ambiguous} or 0);    
            }
            my $discPos = "$chrom:$pos";  #initialize readCounts to 0 for discrepancy positions
            ${$readCounts{1}}{$discPos} = 0;
            ${$readCounts{2}}{$discPos} = 0;
        }         
    }
    #increment consistent counts when discrepancy positively identified in only one of two samples
    DISC: foreach my $disc(keys %discs){    
        my ($chrom, $pos, $discType, $extra) = split(":", $disc);  #retrieve discrepancy information
        defined $mismatchTypes{$discType} or next;     
        foreach my $refSampleN(1..2){   
            my $testSampleN = ($refSampleN % 2) + 1;
            if ($discCounts{$refSampleN}{$disc} and !$consistentCounts{$testSampleN}{$disc}){    
                my $ambiguous = getDisc($chrom, $pos, $types{Discs}{N}, $extra);     
                ${$consistentCounts{$testSampleN}}{$disc} = (${$discCounts{$testSampleN}}{$ambiguous} or 0);     
                next DISC;
            }        
        }
    }            
    #scan for number of reads crossing discrepancy positions
    fillReadCounts($sample1, $readCounts{1});
    fillReadCounts($sample2, $readCounts{2});    
    #commit discrepancies
    status("  committing discrepancies...\n");
    my %chromSeqs;
    loadChromosomes(\%chromSeqs);
    my $cdsTable = "\UCDS_$param{refSeq}";
    tableExists($cdsTable) or $cdsTable = 0;
    my $nDisc;
    my $timeStamp = getTimeStamp();
    my $discsTable = newTable('Discs', "$sample1\_$sample2");
    my $discsFile = "$discsTable.csv";
    open my $discsFileH, ">", $discsFile;
    foreach my $disc(keys %discs){
        my ($chrom, $pos, $discType, $extra) = split(":", $disc);  #retrieve discrepancy information
        my $s1_consistentCount = (${$consistentCounts{1}}{$disc} or 0);
        my $s2_consistentCount = (${$consistentCounts{2}}{$disc} or 0); 
        if (($s1_consistentCount + $s2_consistentCount) > 1){
            my $s1_sameCount = (${$discCounts{1}}{$disc} or 0);
            my $s2_sameCount = (${$discCounts{2}}{$disc} or 0);        
            $nDisc++;        
            my $discID = ($nDisc * 1E8) + $timeStamp;
            my $discPos = "$chrom:$pos";  
            $s1_consistentCount <= ${$readCounts{1}}{$discPos} or $s1_consistentCount = ${$readCounts{1}}{$discPos};  #fixes rare error in consistent count (I think propagated from Pass)
            $s2_consistentCount <= ${$readCounts{2}}{$discPos} or $s2_consistentCount = ${$readCounts{2}}{$discPos};  #fixes rare error in consistent count (I think propagated from Pass)
            my $sameLOD = calculateLOD( (${$readCounts{1}}{$discPos} - $s1_sameCount), $s1_sameCount,
                                        (${$readCounts{2}}{$discPos} - $s2_sameCount), $s2_sameCount );
            my $consistentLOD = calculateLOD( (${$readCounts{1}}{$discPos} - $s1_consistentCount), $s1_consistentCount,
                                              (${$readCounts{2}}{$discPos} - $s2_consistentCount), $s2_consistentCount );        
            my $consequence = getConsequence($chrom, $pos, $discType, \%chromSeqs, $cdsTable);
            print $discsFileH join(",", ($discID,
                                         $chrom, $pos, $discType, $extra,
                                         $s1_sameCount, $s1_consistentCount, ${$readCounts{1}}{$discPos},
                                         $s2_sameCount, $s2_consistentCount, ${$readCounts{2}}{$discPos},
                                         $sameLOD, $consistentLOD, $consequence))."\n"; 
        }  
    }
    close $discsFileH;
    loadData($discsFile, $discsTable, ",", "DISCREPANCYID, 
                                            CHROMOSOME, POSITION, DISCREPANCYTYPE, EXTRA,
                                            S1_SAMECOUNT, S1_CONSISTENTCOUNT, S1_READCOUNT, 
                                            S2_SAMECOUNT, S2_CONSISTENTCOUNT, S2_READCOUNT,
                                            SAMELOD, CONSISTENTLOD, CONSEQUENCE");    
    purgeAmbiguous($discsTable);                                       
}
 
sub initializeTypes{ #cannot be done at compile
    %allowedTypes = ($types{Discs}{A}=>0, $types{Discs}{C}=>0,
                        $types{Discs}{G}=>0, $types{Discs}{T}=>0,
                        $types{Discs}{Missing}=>0, $types{Discs}{Extra}=>0);
    %mismatchTypes = ($types{Discs}{A}=>0, $types{Discs}{C}=>0,
                        $types{Discs}{G}=>0, $types{Discs}{T}=>0);
}
 
sub fillDiscCounts{
    my ($sample, $discCountsRef) = @_;
    my $fragsTable = getTableName('Frags', $sample);   
    fillDiscCounts_Read($fragsTable, 1, $discCountsRef);
    fillDiscCounts_Read($fragsTable, 2, $discCountsRef); 
}

sub fillDiscCounts_Read{
    my ($fragsTable, $read, $discCountsRef) = @_;
    status("  filling discrepancy counts, $fragsTable, read $read...\n");
    #DISCREPANCIES2 = 0 for SingleRead, non-existent read2 will not be retrieved
    runSQL("SELECT CHROMOSOME, POSITION, DISCREPANCIES
            FROM 
            (SELECT POSITION$read POSITION, (DISCREPANCIES$read||'x') DISCREPANCIES,
                     (CASE $read
                        WHEN 1 THEN CHROMOSOME1
                        WHEN 2 THEN (CASE CHROMOSOME2
                                    WHEN 0 THEN CHROMOSOME1
                                    ELSE CHROMOSOME2 END) END) CHROMOSOME
             FROM $fragsTable
             WHERE (FRAGMENTTYPE = $types{Frags}{Normal}
                     OR FRAGMENTTYPE = $types{Frags}{ReverseNormal}
                     OR FRAGMENTTYPE = $types{Frags}{SingleRead}
                     OR (NSETSPAIR = 1 AND NSETSFRAG = 1)) 
                 AND DISCREPANCIES$read > 0 )
            WHERE CHROMOSOME <= $refSeqs{$param{refSeqBase}}{nChrom}",
            \my($chrom, $position, $discs) ); #'x' forces Perl to use Discrepancies as a string
    while (fetchRow()){ #parse out, collect, and count all discrepancies    
        $discs =~ m/(.*)x$/; #strip the 'x' appended to force Perl to use $discreps as string
        while ($1){
            #$1 = discreps remaining; $2 = relPos, $3 = discType, $4 = extra
            $1 =~ m/(.*)(..)(.)(.)$/ or $1 =~ m/(.*)(.)(.)(.)$/;
            if($3){
                my $pos = $position + $2 - 1;
                my $disc = getDisc($chrom, $pos, $3, $4);
                $$discCountsRef{$disc}++; #count the number of times this discrepancy postion and value occurred            
            }
        }
    }
}

sub fillReadCounts{
    my ($sample, $readCountsRef) = @_;
    my $fragsTable = getTableName('Frags', $sample);   
    fillReadCounts_Read($fragsTable, 1, $readCountsRef);
    fillReadCounts_Read($fragsTable, 2, $readCountsRef); 
}

sub fillReadCounts_Read{
    my ($fragsTable, $read, $readCountsRef) = @_;
    status("  filling read counts, $fragsTable, read $read...\n");
    foreach my $chrom (1..$refSeqs{$param{refSeqBase}}{nChrom}){
        #POSITION2 = 0 for SingleRead, non-existent read2 will not be retrieved
        runSQL("SELECT POSITION, LENGTH_
                FROM 
                (SELECT POSITION$read POSITION, LENGTH$read LENGTH_, 
                         (CASE $read
                            WHEN 1 THEN CHROMOSOME1
                            WHEN 2 THEN (CASE CHROMOSOME2
                                        WHEN 0 THEN CHROMOSOME1
                                        ELSE CHROMOSOME2 END) END) CHROMOSOME
                 FROM $fragsTable
                 WHERE (FRAGMENTTYPE = $types{Frags}{Normal}
                         OR FRAGMENTTYPE = $types{Frags}{ReverseNormal} 
                         OR FRAGMENTTYPE = $types{Frags}{SingleRead}                
                         OR (NSETSPAIR = 1 AND NSETSFRAG = 1)) )
                WHERE CHROMOSOME = $chrom
                  AND POSITION > 0",
                \my($position, $length) ); #'x' forces Perl to use Discrepancies as a string
        while (fetchRow()){
            for my $pos ($position..($position + $length - 1)){
                my $discPos = "$chrom:$pos";
                if (defined $$readCountsRef{$discPos}){$$readCountsRef{$discPos}++}
            }
        }          
    }
}

sub getDisc{
    my ($chrom, $pos, $discType, $extra) = @_;
    return join(":", $chrom, $pos, $discType, $extra);
}

sub calculateLOD{
    my ($wtInwt, $mutInwt, $wtInmut, $mutInmut) = @_;
    my $NR = $wtInwt + $mutInmut;
    my $R = $mutInwt + $wtInmut;
    if (!($NR + $R)){return -1} #shouldn't happen, but just in case passed zeros
    my $theta = $R / ($NR + $R);
    my $pWithLinkage = ((1 - $theta)**$NR) * ($theta**$R);
    my $pWithoutLinkage = 0.5**($NR + $R);
    if (!$pWithoutLinkage){return -1} #happens when NR and R get too big and pWithoutLinkage overflows variable
    my $LOD = log10($pWithLinkage/$pWithoutLinkage);
    return $LOD;
}

sub log10 {
    my ($n) = @_; 
    my $returnValue = -1;
    eval{$returnValue = log($n)/log(10)};
    return $returnValue;
}

sub loadChromosomes{
    my ($chromSeqsRef) = @_;
    my $refSeqFile = "$param{refSeqPath}/$param{refSeq}/$param{refSeq}.fa";
    open my $refSeqFileH, "<", $refSeqFile;
    my $chrom = 0;
    while (<$refSeqFileH>){
        chomp $_;
        if ($_ =~ m/^>(.+)/){
            $chrom = $refSeqs{$param{refSeq}}{$1};  
        } elsif ($chrom){
            $$chromSeqsRef{$chrom} .= $_; 
        }   
    }
    close $refSeqFileH;
}

##############################################
sub fixConsequences{
    my ($sample1, $sample2) = @_;
    initializeTypes();
    my %chromSeqs;
    loadChromosomes(\%chromSeqs);       
    my $cdsTable = "\UCDS_$param{refSeq}";    
    my $discsTable = getTableName('Discs', "$sample1\_$sample2");
    my $discsTableTMP = newTable('Discs', "$sample1\_$sample2\_TMP");  
    my $discsFileTMP = "$discsTableTMP.csv";
    open my $discsTMPH, ">", $discsFileTMP;
    runSQL("SELECT * FROM $discsTable");   
    my @discRefs;
    while (my $__ = fetchRowHashRef()){ push @discRefs, $__ }
    foreach my $__(@discRefs){
        my $consequence = getConsequence($$__{CHROMOSOME}, $$__{POSITION}, $$__{DISCREPANCYTYPE}, \%chromSeqs, $cdsTable);     
        print $discsTMPH join(",", $$__{DISCREPANCYID}, 
                                  $$__{CHROMOSOME}, $$__{POSITION}, $$__{DISCREPANCYTYPE}, $$__{EXTRA},
                                  $$__{S1_SAMECOUNT}, $$__{S1_CONSISTENTCOUNT}, $$__{S1_READCOUNT}, 
                                  $$__{S2_SAMECOUNT}, $$__{S2_CONSISTENTCOUNT}, $$__{S2_READCOUNT},
                                  $$__{SAMELOD}, $$__{CONSISTENTLOD}, $consequence )."\n"; 
    }
    close $discsTMPH;
    loadData($discsFileTMP, $discsTableTMP, ",", "DISCREPANCYID, 
                                            CHROMOSOME, POSITION, DISCREPANCYTYPE, EXTRA,
                                            S1_SAMECOUNT, S1_CONSISTENTCOUNT, S1_READCOUNT, 
                                            S2_SAMECOUNT, S2_CONSISTENTCOUNT, S2_READCOUNT,
                                            SAMELOD, CONSISTENTLOD, CONSEQUENCE");   
}
##############################################


sub getConsequence {
    my ($chrom, $pos, $discType, $chromSeqsRef, $cdsTable) = @_;
    $cdsTable or return -1; #flag indicating no gene table was available
    $discType == $types{Discs}{N} and return $types{Consequences}{Silent}; #cannot determine consequence of ambiguous bases
    unless (defined $mismatchTypes{$discType}){return $types{Consequences}{InDel}}  #non-mismatch in cds returns InDel
    runSQL("SELECT START_, CARRYOVER, STRAND  FROM $cdsTable
           WHERE CHROMOSOME = $chrom
            AND START_ <= $pos
            AND END_ >= $pos",
            \my($start, $carryover, $strand));
    while (fetchRow()){  #loop through crossing cds to see if any change coding
        # +1 in the next line accounts for fact that UCSC gene table
        # is zero-referenced, i.e. the 1st base is called 0
        my $start = $start - $carryover + 1  ;
        my $framePos = ($pos - $start) % 3;
        my $codonStart = $pos - $framePos;
        # -1 in the following line needed to get substr to
        #return the correct base.  not sure why, is this an error in reading GFF or correct use of substr...
        my $wtCodon = substr($$chromSeqsRef{$chrom}, $codonStart - 1, 3);
        my $wtAA = translateCodon($wtCodon, $strand);
        $wtCodon =~ m/^(.{$framePos})(.)(.*)/; #replace mutation in codon
        my $mutCodon = "$1$types{RevDiscs}{$discType}$3";
        my $mutAA = translateCodon($mutCodon, $strand);
        unless ($mutAA eq $wtAA){ #do nothing yet for silent mutations
            if ($mutAA eq 'x'){ #new value is a stop codon
                return $types{Consequences}{Nonsense} 
            } else {  #new value is not a stop codon, i.e. another aa                
                return $types{Consequences}{Missense} 
            }
        } 
    }        
    return $types{Consequences}{Silent};  #either no crossing cds OR no non-silent mutation encountered in cds
}

sub translateCodon{
    my ($codon, $strand) = @_;
    if ($strand - 1){$codon = reverseComplement($codon)}
    return ($geneticCode{$codon} or 'X')
}

sub reverseComplement{
    my ($sequence) = @_;
    $sequence = reverse "\U$sequence";
    $sequence =~ tr/ACGT/TGCA/;
    return $sequence;
}

sub purgeAmbiguous{ #only keep ambiguous discrepancies if a mismatch was not also positiviely identified at that position
    my ($discsTable) = @_;
    my $positionIndex = "CHROMOSOME || 'x' || POSITION";
    my $mismatchIndexes = "SELECT $positionIndex
                           FROM $discsTable
                           WHERE DISCREPANCYTYPE = $types{Discs}{A}
                              OR DISCREPANCYTYPE = $types{Discs}{C}
                              OR DISCREPANCYTYPE = $types{Discs}{G}
                              OR DISCREPANCYTYPE = $types{Discs}{T}";
    runSQL("DELETE FROM $discsTable
            WHERE DISCREPANCYTYPE = $types{Discs}{N}
              AND $positionIndex IN ($mismatchIndexes)");
}

1;













