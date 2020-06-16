#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %gffFields %geneticCode));

my (%stats1, %stats2);

sub findHomozygous{
    my ($childRefSeq, $sample1, $sample2) = @_;

    status("getting parameters...\n");
        my $sample = $sample1;
        $sample2 and $sample .= "_$sample2"; 
        my $discsTable = getTableName('Discs', $sample);
        tableExists($discsTable) or die "could not find table $discsTable";
        getStatistics(getTableName('Stats', $sample1), \%stats1);        
        $sample2 and getStatistics(getTableName('Stats', $sample2), \%stats2);  
        my $lineSize = getLineSize(); 

    status("retrieving homozygous discrepancies...\n"); 
        my $mutsRef = getHomoMut($discsTable, $lineSize, $sample2);          
       
    status("finding missing/extra bases consistent with homozygosity...\n");        
        my $chromLinesRef = loadChromLines($param{refSeq});        
        checkOtherME($discsTable, $sample1, $sample2, $mutsRef, $chromLinesRef, $lineSize, \&checkMEHomozygous, \%stats1, \%stats2);
          
    status("creating mutated copy of $param{refSeq}...\n");
        mutateChromSeqs($childRefSeq, $mutsRef);  
        rectifyRefSeqLines($childRefSeq, $lineSize);               
        my $childChromLinesRef = loadChromLines($childRefSeq); 
        printChromFiles($childChromLinesRef, $childRefSeq); 

    status("printing table of mutations by position...\n");
        my $mutsByGeneRef = printMutations($childRefSeq, 'Homozygous', $mutsRef);    

    status("printing table of mutations by gene and generating CDS table...\n");
        printGenes($childRefSeq, $lineSize, $chromLinesRef, $childChromLinesRef, 'Homozygous', $mutsByGeneRef);    
}
        
sub getLineSize{   
    my $refSeqFile = "$param{refSeqPath}/$param{refSeq}/$param{refSeq}.fa";    
    open my $refSeqH, "<", $refSeqFile;
    my $chrom = 0;
    while (<$refSeqH>){
        if ($_ =~ m/^>.+\n/){
            $chrom = 1;
        } elsif ($chrom){
            chomp $_;
            my $length = length $_;
            if ($length){return $length}
        }   
    }
    close $refSeqH;
}

sub getHomoMut{
    my ($discsTable, $lineSize, $sample2) = @_;
    my $sampleFilters = getSampleFilters(\%stats1, \%stats2);
    my $sample2Filter = '';
    $sample2 and $sample2Filter = ' AND S2_CONSISTENTCOUNT > 0 ';
    #homozygosity logic = homozygous in combined sample but only if detected in all samples
    return getMutations($lineSize,
                        "SELECT CHROMOSOME, POSITION, DISCREPANCYTYPE, EXTRA,
                                S1_CONSISTENTCOUNT, S1_READCOUNT, S2_CONSISTENTCOUNT, S2_READCOUNT, CONSISTENTLOD
                         FROM $discsTable
                         WHERE ((S1_CONSISTENTCOUNT + S2_CONSISTENTCOUNT)/(S1_READCOUNT + S2_READCOUNT)) >= $param{minHomoF}
                           AND S1_CONSISTENTCOUNT > 0 $sample2Filter
                           AND $sampleFilters");    
}

sub getMutations{
    my ($lineSize, $sql) = @_; 
    runSQL($sql, \my($chrom, $pos, $discType, $extra, $s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD));
    my %newMuts;
    while (fetchRow()){ 
        fillMutation($chrom, $pos, $discType, $extra, $s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD, $lineSize, \%newMuts) 
    }
    return \%newMuts;
}

sub fillMutation{
    my ($chrom, $pos, $discType, $extra, $s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD, $lineSize, $mutsRef) = @_;
    my ($zrLine, $zrLinePos) = getzrValues($pos, $lineSize);
    $$mutsRef{$chrom}{$zrLine}{$zrLinePos} = {Position => $pos, DiscType => $discType, Extra => $extra,
                                              S1Consistent => $s1Consistent, S1Reads => $s1Reads, S2Consistent => $s2Consistent, S2Reads => $s2Reads, LOD => $LOD}    
}

sub getzrValues{
    my ($pos, $lineSize) = @_;
    my $zrPos = $pos - 1;
    my $zrLine = int($zrPos / $lineSize);
    my $zrLinePos = $zrPos % $lineSize;
    return ($zrLine, $zrLinePos);
}

sub checkOtherME{
    my ($discsTable, $sample1, $sample2, $mutsRef, $chromLinesRef, $lineSize, $checkMESubRef, $stats1Ref, $stats2Ref) = @_;
    runSQL(getOtherMESQL($discsTable, $sample1, $sample2, $stats1Ref, $stats2Ref));
    my $readRef = fetchRowHashRef();     
    my @reads;   
    my $prevReadRef = $readRef;                                                   
    push @reads, $readRef;
    while (my $readRef = fetchRowHashRef()){
        unless ($$prevReadRef{DISCREPANCYID} eq $$readRef{DISCREPANCYID}){
            my ($chrom, $pos, $discType, $extra) = ($$prevReadRef{CHROMOSOME}, $$prevReadRef{POSITION}, $$prevReadRef{DISCREPANCYTYPE}, $$prevReadRef{EXTRA});
            my %nReads = ($sample1 => 0); $sample2 and $nReads{$sample2} = 0;
            my %consistent = ($sample1 => 0); $sample2 and $consistent{$sample2} = 0;
            READ_ : foreach my $readRef(@reads){
                $nReads{$$readRef{SAMPLE}}++;  
                my %discs;             
                $$readRef{DISCREPANCIES} =~ m/(.*)x$/; #strip the 'x' appended to force Perl to use $discreps as string 
                while ($1){ #$1 = discreps remaining; $2 = relPos, $3 = discType, $4 = extra 
                    $1 =~ m/(.*)(..)(.)(.)$/ or $1 =~ m/(.*)(.)(.)(.)$/;  
                    if ($3){
                        my $zrRelPos = $2 - 1;           
                        if($pos == ($$readRef{LOWPOSITION} + $zrRelPos) and $discType == $3 and $extra == $4){
                            $consistent{$$readRef{SAMPLE}}++;  
                            next READ_;
                        } else {
                            $discs{$zrRelPos} = {DiscType => $3, Extra => $4};
                            $3 == $types{Discs}{Missing} and $$readRef{HIGHPOSITION}++;; 
                            $3 == $types{Discs}{Extra} and $$readRef{HIGHPOSITION}--; 
                        }                    
                    }   
                }  
                my $read = getDNASegment($chrom, $$readRef{LOWPOSITION}, $$readRef{HIGHPOSITION}, $lineSize, $chromLinesRef);
                $read or ($nReads{$$readRef{SAMPLE}}-- and next READ_);
                my $mutRead = mutateDNASegmentMulti($read, \%discs);       
                my ($left, $right, $readEnd, $refEnd, $length) = ($pos - $$readRef{LOWPOSITION} + 1, $$readRef{HIGHPOSITION} - $pos + 1);  
                if ($left <= $right){ #left 
                    $length = $left;             
                    $readEnd = substr($mutRead, 0, $length);          
                    if ($discType == $types{Discs}{Missing}){ #Missing                       
                        $refEnd = getDNASegment($chrom, $pos - $length, $pos - 1, $lineSize, $chromLinesRef);
                    } else { #Extra
                        if($length >= 2){
                            $refEnd = getDNASegment($chrom, $pos - $length + 2, $pos, $lineSize, $chromLinesRef).$types{RevDiscs}{$extra}; 
                        } else {
                            $refEnd = $types{RevDiscs}{$extra};
                        }                 
                    }   
                } else { #right   
                    $length = $right;                           
                    $readEnd = substr($mutRead, -$length);
                    if ($discType == $types{Discs}{Missing}){   #Missing                      
                        $refEnd = getDNASegment($chrom, $pos + 1, $pos + $length, $lineSize, $chromLinesRef);                     
                    } else { #Extra
                        if($length >= 2){
                            $refEnd = getDNASegment($chrom, $pos, $pos + $length - 2, $lineSize, $chromLinesRef);
                        } else {
                            $refEnd = getDNASegment($chrom, $pos, $pos, $lineSize, $chromLinesRef);
                        }
                        $refEnd =~ m/^(.{1})(.*)/;
                        $refEnd = "$1$types{RevDiscs}{$extra}$2";
                    }
                }   
                ($readEnd and length($readEnd) == $length) or ($nReads{$$readRef{SAMPLE}}-- and next READ_);                
                ($refEnd and length($refEnd) == $length) or ($nReads{$$readRef{SAMPLE}}-- and next READ_);
                my @readEnd = split(//, $readEnd);                   
                my @refEnd = split(//, $refEnd); 
                foreach my $i(0..($length - 1)){ unless($readEnd[$i] eq $refEnd[$i] or $readEnd[$i] eq 'N'){ next READ_ } }
                $consistent{$$readRef{SAMPLE}}++;  
            }     
            my ($passed, $s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD) = &$checkMESubRef($sample1, $sample2, \%nReads, \%consistent);
            $passed and fillMutation($chrom, $pos, $discType, $extra, $s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD, $lineSize, $mutsRef);
            @reads = ();                     
        }   
        $prevReadRef = $readRef;  
        push @reads, $readRef;
    }   
}

sub getOtherMESQL{
    my ($discsTable, $sample1, $sample2, $stats1Ref, $stats2Ref) = @_;
    my $discsSQL = getMEDiscsSQL($discsTable, $stats1Ref, $stats2Ref);  
    my $readsSQL = getMEReadsSQL($sample1);                               
    $sample2 and $readsSQL .= " UNION ALL " . getMEReadsSQL($sample2);                              
    return "SELECT DISCREPANCYID, CHROMOSOME, POSITION, DISCREPANCYTYPE, EXTRA,
                   SAMPLE, LOWPOSITION, HIGHPOSITION, DISCREPANCIES
            FROM
            (SELECT d.DISCREPANCYID, d.CHROMOSOME, d.POSITION, d.DISCREPANCYTYPE, d.EXTRA,
                    r.SAMPLE, r.LOWPOSITION, r.HIGHPOSITION, r.DISCREPANCIES
            FROM ($discsSQL) d, ($readsSQL) r
            WHERE d.CHROMOSOME = r.CHROMOSOME
              AND d.POSITION >= r.LOWPOSITION
              AND d.POSITION <= r.HIGHPOSITION )
            ORDER BY DISCREPANCYID" ;  
}

sub getMEDiscsSQL{
    my ($discsTable, $stats1Ref, $stats2Ref) = @_;
    my $sampleFilters = getSampleFilters($stats1Ref, $stats2Ref);
    return "SELECT DISCREPANCYID, CHROMOSOME, POSITION, DISCREPANCYTYPE, EXTRA
            FROM $discsTable
            WHERE (DISCREPANCYTYPE = $types{Discs}{Extra} OR DISCREPANCYTYPE = $types{Discs}{Missing})
               AND ((S1_CONSISTENTCOUNT + S2_CONSISTENTCOUNT)/(S1_READCOUNT + S2_READCOUNT)) < $param{minHomoF}
               AND $sampleFilters ";
}

sub getMEReadsSQL{
    my ($sample) = @_;
    my $fragsTable = getTableName('Frags', $sample);
    my $read1SQL = getMEReadSQL($fragsTable, 1, $sample);
    my $read2SQL = getMEReadSQL($fragsTable, 2, $sample);
    return " $read1SQL UNION ALL $read2SQL ";
}

sub getMEReadSQL{
    my ($fragsTable, $read, $sample) = @_;
    my $chromLabel = " CHROMOSOME1 ";
    $read == 2 and $chromLabel = " CASE CHROMOSOME2 WHEN 0 THEN CHROMOSOME1 ELSE CHROMOSOME2 END ";
    return "SELECT '$sample' SAMPLE, $chromLabel CHROMOSOME, 
                   POSITION$read LOWPOSITION, (POSITION$read + LENGTH$read - 1) HIGHPOSITION, 
                   DISCREPANCIES$read || 'x' DISCREPANCIES
             FROM $fragsTable
             WHERE (FRAGMENTTYPE = $types{Frags}{Normal}
                OR FRAGMENTTYPE = $types{Frags}{ReverseNormal} 
                OR FRAGMENTTYPE = $types{Frags}{SingleRead}                
                OR (NSETSPAIR = 1 AND NSETSFRAG = 1))
               AND POSITION$read > 0 "; #last filter prevents recovery of non-existent read2 for SingleRead Fragments
}

sub getSampleFilters{
    my ($stats1Ref, $stats2Ref, $suppressMin1) = @_;
    my $sampleFilters = getSampleFilter(1, $stats1Ref, $suppressMin1);
    $$stats2Ref{minCoverageRMap} and $sampleFilters .= " AND " . getSampleFilter(2, $stats2Ref);
    return $sampleFilters;
}

sub getSampleFilter{
    my ($sampleN, $statsRef, $suppressMin) = @_;
    my $minFilter = " S$sampleN\_READCOUNT >= $$statsRef{minCoverageRMap} AND ";
    $suppressMin and $minFilter = '';
    return " $minFilter S$sampleN\_READCOUNT <= $$statsRef{maxCoverageRMap} "
}

sub checkMEHomozygous{
    my ($sample1, $sample2, $nReadsRef, $consistentRef) = @_;
    my ($consistent, $nReads, $s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD) =  
        getMECounts($sample1, $sample2, $nReadsRef, $consistentRef);
    $nReads or return undef;        
    my $inBothSamples = ($s1Consistent and (!$sample2 or $s2Consistent));          
    my $countConsistent = (($consistent / $nReads) >= $param{minHomoF}); 
    my $passed = ($inBothSamples and $countConsistent);
    return ($passed, $s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD, $consistent, $nReads);
}   

sub getMECounts{
    my ($sample1, $sample2, $nReadsRef, $consistentRef) = @_;
    my ($s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD) = ($$consistentRef{$sample1}, $$nReadsRef{$sample1}, 0, 0, 0);
    if($sample2){
        ($s2Consistent, $s2Reads) = ($$consistentRef{$sample2}, $$nReadsRef{$sample2});
        my ($wtInwt, $mutInwt, $wtInmut, $mutInmut) = ($s1Reads - $s1Consistent, $s1Consistent, $s2Reads - $s2Consistent, $s2Consistent);
        $LOD = calculateLOD($wtInwt, $mutInwt, $wtInmut, $mutInmut);        
    }
    my ($consistent, $nReads) = ($s1Consistent + $s2Consistent, $s1Reads + $s2Reads);
    return ($consistent, $nReads, $s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD);
}

sub loadChromLines{
    my ($refSeq, $getChrom) = @_; #inputs optional, getChrom restricts return to a specific chromosome
    $refSeq or $refSeq = $param{refSeq};
    my $refSeqFile = "$param{refSeqPath}/$refSeq/$refSeq.fa";
    open my $refSeqFileH, "<", $refSeqFile;
    my ($chrom, $zrLine, %chromLines) = (0, 0);
    while (<$refSeqFileH>){
        my $line = $_;
        $line =~ m/\S/ or next; #skip blank lines                  
        chomp $line;
        if ($line =~ m/^>(.*)/){
            $chrom = $refSeqs{$param{refSeqBase}}{$1};  
            defined $getChrom and ($getChrom == $chrom or $chrom = 0);
            $zrLine = 0;
        } elsif ($chrom){
            $chromLines{$chrom}{$zrLine} = $line; 
            $zrLine++;  
        }   
    }
    close $refSeqFileH;
    return \%chromLines;
}

sub getDNASegment{
    my ($chrom, $minPos, $maxPos, $lineSize, $chromLinesRef) = @_;
    ($minPos > 0 and $maxPos > 0) or return undef;
    my $expectedLength = $maxPos - $minPos + 1;
    ($expectedLength and $expectedLength > 0) or return undef;
    my $segment;      
    my ($minZRLine, $zrLinePos) = getzrValues($minPos, $lineSize); 
    my ($maxZRLine) = getzrValues($maxPos, $lineSize);   
    foreach my $zrLine($minZRLine..$maxZRLine){ $$chromLinesRef{$chrom}{$zrLine} and $segment .= $$chromLinesRef{$chrom}{$zrLine} }
    $segment or return undef;
    my $segmentLength = length($segment); 
    $segmentLength >= $zrLinePos + $expectedLength or return undef;
    $segment = substr($segment, $zrLinePos, $expectedLength); 
    return $segment;
}

sub mutateChromSeqs{
    my ($childRefSeq, $mutsRef) = @_;  
    status("  introducing mutations...\n");  
    my $refSeqInFile = "$param{refSeqPath}/$param{refSeq}/$param{refSeq}.fa";   
    my $refSeqOutPath = "$param{refSeqPath}/$childRefSeq"; 
    -d $refSeqOutPath or mkdir $refSeqOutPath;
    my $refSeqOutFile = "$refSeqOutPath/$childRefSeq.fa"; 
    open my $refSeqInH, "<", $refSeqInFile or die "could not open $refSeqInFile";
    open my $refSeqOutH, ">", $refSeqOutFile or die "could not open $refSeqOutFile";
    my ($chrom, $zrLine) = (0, 0);
    while (<$refSeqInH>){
        my $line = $_;
        $line =~ m/\S/ or next; #skip blank lines        
        chomp $line;
        if ($line =~ m/^>(.*)/){
            $chrom = $refSeqs{$param{refSeqBase}}{$1};        
            $zrLine = 0;       
        } elsif ($chrom){
            if (defined $$mutsRef{$chrom}{$zrLine}){ $line = mutateDNASegmentMulti($line, $$mutsRef{$chrom}{$zrLine}) }
            $zrLine++;  
        }        
        print $refSeqOutH "$line\n";
    }
    close $refSeqInH;
    close $refSeqOutH;
}

sub mutateDNASegmentMulti{
    my ($segment, $discsRef) = @_; 
    foreach my $zrLinePos(sort {$b <=> $a} keys %$discsRef){ #scroll BACKWARDS through dicrepancies!    
        $segment = mutateDNASegment($segment, $zrLinePos, $$discsRef{$zrLinePos}{DiscType}, $$discsRef{$zrLinePos}{Extra})
    } 
    return $segment;
}

sub mutateDNASegment{
    my ($segment, $zrLinePos, $discType, $extra) = @_;
    unless($segment =~ m/^(.{$zrLinePos})(.)(.*)/){
        print "  mutateSegment error on segment $segment, $zrLinePos, $discType, $extra\n"; 
        return $segment;
    }
    my $mutatedBase = $2;
    if ($types{RevDiscs}{$discType}){                       
        $mutatedBase = $types{RevDiscs}{$discType};  
    } elsif ($discType == $types{Discs}{Extra}){
        $mutatedBase = "$2$types{RevDiscs}{$extra}";                     
    } elsif ($discType == $types{Discs}{Missing}){
        $mutatedBase = '';
    } 
    return "$1$mutatedBase$3"; 
}

sub rectifyRefSeqLines {
    my ($refSeq, $lineSize) = @_; 
    status("  rectifying refSeq file to $lineSize-base lines...\n");  
    my $refSeqFile = "$param{refSeqPath}/$refSeq/$refSeq.fa";      
    my $refSeqTmpFile = "$refSeqFile.tmp";       
    system("mv $refSeqFile $refSeqTmpFile");  
    open my $refSeqTmpH, "<", $refSeqTmpFile or die "could not open $refSeqTmpFile";
    open my $refSeqH, ">", $refSeqFile or die "could not open $refSeqFile";
    my $working = '';
    while (<$refSeqTmpH>){     
        my $line = $_;
        $line =~ m/\S/ or next; #skip blank lines   
        chomp $line; 
        if ($line =~ m/^>/){ #pass name lines     
            $working and print $refSeqH "$working\n";
            $working = '';
            print $refSeqH "$line\n";
            next;
        }
        $working .= $line;
        while($working =~ m/^(.{$lineSize})(.*)/){
            print $refSeqH "$1\n";
            $working = $2;
        }
    }
    close $refSeqTmpH;
    close $refSeqH;
    unlink $refSeqTmpFile;  
}

sub printChromFiles{
    my ($chromLinesRef, $refSeq) = @_;
    my $outputPath = "$param{refSeqPath}/$refSeq";
    -d $outputPath or mkdir $outputPath;
    foreach my $chrom(keys %$chromLinesRef){
        my $chr = $reverseRefSeqs{$param{refSeqBase}}{$chrom};
        $chr or next;
        my $chromFile = "$outputPath/$chr.fa";
        open my $chromFileH, ">", $chromFile;
        print $chromFileH ">$chr\n";
        foreach my $zrLine(sort {$a <=> $b} keys %{$$chromLinesRef{$chrom}}){ 
            print $chromFileH "$$chromLinesRef{$chrom}{$zrLine}\n" 
        }
        close $chromFileH;
    }
}

sub printMutations{
    my ($childRefSeq, $resultsFileName, $mutsRef) = @_;
    my %discLabels = (1=>'Mismatch', 5=>'Ambiguous', 6=>'Missing', 7=>'Extra');
    $discLabels{2} = $discLabels{1}; $discLabels{3} = $discLabels{1}; $discLabels{4} = $discLabels{1};
    my $mutsFile = "$param{refSeqPath}/$childRefSeq/$childRefSeq\_$resultsFileName\_Mutations.xls";
    open my $mutsFileH, ">", $mutsFile;    
    my @labels = qw(Candidate Chromosome Position Mutation Added_Base Counts(S1Mut/Tot//S2Mut/Tot) LOD Affected_Genes);
    print $mutsFileH join("\t", @labels)."\n";    
    my (%mutsByGene, %discCounts);
    foreach my $chrom(sort {$a <=> $b} keys %$mutsRef){     
        my $chr = $reverseRefSeqs{$param{refSeqBase}}{$chrom};
        print $mutsFileH "\t$chr\n";   
        foreach my $zrLine(sort {$a <=> $b} keys %{$$mutsRef{$chrom}}){
            foreach my $zrLinePos(sort {$a <=> $b} keys %{$$mutsRef{$chrom}{$zrLine}}){  
                my $candidate = ' ';
                $$mutsRef{$chrom}{$zrLine}{$zrLinePos}{isCandidate} and $candidate = "CANDIDATE";
                my $pos = $$mutsRef{$chrom}{$zrLine}{$zrLinePos}{Position};
                my $discType = $$mutsRef{$chrom}{$zrLine}{$zrLinePos}{DiscType};
                $discCounts{$discLabels{$discType}}++;        
                $types{RevDiscsFull}{$discType} and $discType = $types{RevDiscsFull}{$discType};
                my $extra = $$mutsRef{$chrom}{$zrLine}{$zrLinePos}{Extra};    
                $types{RevDiscsFull}{$extra} and $extra = $types{RevDiscsFull}{$extra}; 
                my ($s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD) = 
                        ($$mutsRef{$chrom}{$zrLine}{$zrLinePos}{S1Consistent},
                         $$mutsRef{$chrom}{$zrLine}{$zrLinePos}{S1Reads},
                         $$mutsRef{$chrom}{$zrLine}{$zrLinePos}{S2Consistent},
                         $$mutsRef{$chrom}{$zrLine}{$zrLinePos}{S2Reads},
                         $$mutsRef{$chrom}{$zrLine}{$zrLinePos}{LOD}); 
                my $counts = "$s1Consistent/$s1Reads // $s2Consistent/$s2Reads";   
                my @mutInfo = ($candidate, $chrom, $pos, $discType, $extra, $counts, $LOD);          
                my @affectedGenes;                
                my $genesRef = getGenes($chrom, $pos);               
                foreach my $geneRef(@$genesRef){
                    my $url = "http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=$$geneRef{NAME1}";               
                    my $hyperlink = "=HYPERLINK(\"$url\",\"$$geneRef{NAME2}\")";                  
                    push @affectedGenes, $hyperlink;
                    push @{$mutsByGene{$$geneRef{NAME1}}}, @mutInfo;
                } 
                print $mutsFileH join("\t", @mutInfo, @affectedGenes)."\n";              
            }
        }
    }
    close $mutsFileH;
    print "  mutations committed to $childRefSeq\:\n";
    foreach my $discLabel(keys %discCounts){print "    $discLabel -> $discCounts{$discLabel}\n"}      
    print "  see $mutsFile\n";
    return \%mutsByGene;
}

sub getGenes{
    my ($chrom, $pos) = @_;
    my $cdsTable = "CDS_$param{refSeq}";
    runSQL("SELECT * FROM $cdsTable
           WHERE CHROMOSOME = $chrom
            AND START_ <= $pos
            AND END_ >= $pos");
    return fetchAllHashRef();
}

sub printGenes{
    my ($childRefSeq, $lineSize, $chromLinesRef, $childChromLinesRef, $resultsFileName, $mutsByGeneRef) = @_;
    my $cdsTableIn = getTableName('CDS', $param{refSeq});
    my $cdsTableOut = newTable('CDS', $childRefSeq);
    my $cdsFileOut = "$cdsTableOut.csv";
    open my $cdsFileOutH, ">", $cdsFileOut or die "could not open $cdsFileOut";
    my $genesFile = "$param{refSeqPath}/$childRefSeq/$childRefSeq\_$resultsFileName\_Genes.xls";
    open my $genesFileH, ">", $genesFile;
    my @labels = qw(Candidate Name1 Name2 Exon Chromosome Chromosome Start End DNA_Changes Protein_Changes DNA Protein Mutations(Candidate/Chr/Pos/Mut/Added/Counts/LOD...));
    print $genesFileH join("\t", @labels)."\n";
    my ($cdsDNADelta, $cdsProtDelta) = (0, 0);
    runSQL("SELECT CDSID, NAME1, NAME2, CHROMOSOME, STRAND, EXON, START_, CARRYOVER, END_
            FROM $cdsTableIn",
            \my($cdsID, $name1, $name2, $chrom, $strand, $exon, $start, $carryover, $end));   
    while (fetchRow()){
        my ($nDNADiscs, $nProtDiscs, $parentStart, $parentEnd, $childStart, $childEnd);
        my $parentDNA = getDNASegment($chrom, $start + 1, $end, $lineSize, $chromLinesRef);
        $parentDNA or next;
        ($nDNADiscs, $parentStart, $parentEnd, $childStart, $childEnd) = blastParent2Child($childRefSeq, 'blastn', $chrom, $parentDNA);        
        defined $nDNADiscs or next;  #too short for blast           
        print $cdsFileOutH join(",", $cdsID, $name1, $name2, $chrom, $strand, $exon, $childStart - 1, $carryover, $childEnd)."\n"; 
        if($nDNADiscs){
            $cdsDNADelta++;
            my $parentProtein = translateGene($parentDNA, $strand);        
            $parentProtein =~ m/(.*)x$/ and $parentProtein = $1;
            my $childDNA = getDNASegment($chrom, $childStart, $childEnd, $lineSize, $childChromLinesRef);
            $childDNA or next;
            my $childProtein = translateGene($childDNA, $strand);              
            $childProtein =~ m/(.*)x$/ and $childProtein = $1;
            ($nProtDiscs) = blastParent2Child($childRefSeq, 'blastp', $chrom, $parentProtein, $childProtein); 
            $nProtDiscs and $cdsProtDelta++;
            my $chr = $reverseRefSeqs{$param{refSeqBase}}{$chrom};
            my $candidate = ' ';
            foreach my $value(@{$$mutsByGeneRef{$name1}}){if($value eq "CANDIDATE"){$candidate = "CANDIDATE";last}}
            my $mutations = ' ';
            $$mutsByGeneRef{$name1} and $mutations = join("\t", @{$$mutsByGeneRef{$name1}});            
            my $url = "http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=$name1";    
            $name1 = "=HYPERLINK(\"$url\",\"$name1\")";         
            $name2 = "=HYPERLINK(\"$url\",\"$name2\")";    
            print $genesFileH join("\t", $candidate, $name1, $name2, $exon, $chr, $chrom, $start, $end, 
                                         $nDNADiscs, $nProtDiscs, $childDNA, $childProtein, $mutations)."\n";
        }                  
    }
    close $cdsFileOutH;
    close $genesFileH;    
    loadData($cdsFileOut, $cdsTableOut, ",", "CDSID, NAME1, NAME2, CHROMOSOME, STRAND, EXON, START_, CARRYOVER, END_");
    print "  $cdsDNADelta coding exons have $resultsFileName DNA changes\n";
    print "  $cdsProtDelta coding exons have protein coding changes\n";
    print "  see $genesFile\n";
}

sub translateGene{
    my ($sequence, $strand) = @_;
    if ($strand - 1){$sequence = reverseComplement($sequence)}
    my $translation;
    while ($sequence =~ m/^(.{3})(.*)/){
        my $codon = "\U$1";
        $translation .= ($geneticCode{$codon} or 'X');
        $sequence = $2;
    }
    return $translation;  
}

sub blastParent2Child{
    my ($childRefSeq, $blastCommand, $chrom, $parentSeq, $childSeq) = @_;
    my $timeStamp = getTimeStamp();
    my $queryFile = "$param{blastPath}$blastCommand\_$childRefSeq\_$timeStamp\_query.seq";     
    open my $queryFileH, ">", $queryFile;    
    print $queryFileH $parentSeq;
    close $queryFileH;     
    my $subjectFile = "$param{refSeqPath}/$childRefSeq/$reverseRefSeqs{$param{refSeqBase}}{$chrom}.fa";      
    if($childSeq){
        $subjectFile = "$param{blastPath}$blastCommand\_$childRefSeq\_$timeStamp\_subject.seq"; 
        open my $subjectFileH, ">", $subjectFile;    
        print $subjectFileH $childSeq;
        close $subjectFileH; 
    }
    my ($mismatch, $gapOpen, $pStart, $pEnd, $chStart, $chEnd) = 
        blast2Seqs($blastCommand, $queryFile, $subjectFile);
    unlink $queryFile;
    $childSeq and unlink $subjectFile;
    defined $mismatch or return undef;
    my $before = $pStart - 1;
    my $after = length($parentSeq) - $pEnd;     
    $chStart -= $before;
    $chEnd += $after;      
    my $nDiscs = ($mismatch + $gapOpen + $before + $after);      
    return ($nDiscs, $pStart, $pEnd, $chStart, $chEnd);  
}

sub blast2Seqs{
    my ($blastCommand, $queryFile, $subjectFile, $subjectStart, $subjectEnd, $options) = @_;
    my $resultsFile = "$queryFile\_blastResults.txt"; 
    my $blast = "$param{blastPath}$blastCommand -max_target_seqs 1 ";
    $blast .= "-query $queryFile -subject $subjectFile -out $resultsFile ";
    $blast .= "-outfmt \"10 mismatch gapopen qstart qend sstart send\" "; 
    $subjectEnd and $blast .= "-subject_loc $subjectStart-$subjectEnd ";
    $options and $blast .= "$options ";
    system($blast);
    open my $resultsFileH, "<", $resultsFile or return undef;
    my $bestMatch = <$resultsFileH>;
    close $resultsFileH; 
    unlink $resultsFile;
    $bestMatch or return undef;   
    chomp $bestMatch;     
    my ($mismatch, $gapOpen, $qstart, $qend, $sstart, $send) = split(",", $bestMatch); 
    return ($mismatch, $gapOpen, $qstart, $qend, $sstart, $send);
}

1;



