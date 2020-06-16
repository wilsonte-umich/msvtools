#!/usr/bin/perl -w
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{upstreamBins} or $param{upstreamBins} = 5; #number of bins preceding the gene start
defined $param{charPerLine} or $param{charPerLine} = 100; 
defined $command{mapRegion} or $command{mapRegion} = ['singleThread', '24:00:00', 5000, 0];
defined $command{mapRegionRooted} or $command{mapRegionRooted} = ['singleThread', '1:00:00', 2000, 0];
defined $command{mapRegionRooted2} or $command{mapRegionRooted2} = ['singleThread', '1:00:00', 2000, 0];

sub mapRegion { 
    my ($inSample, $targetName, $chrom, $start, $end) = @_;
    
    status("creating refSeq for $targetName...\n");
        createRegionRefSeq($targetName, $chrom, $start, $end);
        createDatabase($targetName);
        $refSeqs{$targetName} = {nChrom => 1,$targetName => 1};
        $reverseRefSeqs{$targetName} = {1 => $targetName};  
        $param{refSeq} = $targetName;
        $param{refSeqBase} = $targetName;
    
    status("mapping $inSample to $targetName...\n");
        my $outSample = "$inSample\_$targetName";
        $param{mapType} = 'pass';
        getDirectories($outSample, \my%outDirs);
        -d $outDirs{sample} or mkdir $outDirs{sample};
        -d $outDirs{read1} or mkdir $outDirs{read1};
        -d $outDirs{read2} or mkdir $outDirs{read2};
        mapReads($inSample, $outSample, 'read1');  
        mapReads($inSample, $outSample, 'read2'); 

    status("extracting pairs...\n");   
        $param{extractDataOnly} = 1;
        extractPairs($outSample);
        getMapRootedStats($inSample, $outSample);
 
    status("parsing fragments...\n");  
        $param{parseFragmentsOnly} = 1;
        parseFragments($outSample);
}

sub createRegionRefSeq {
    my ($targetName, $chrom, $start, $end) = @_;
    my $lineSize = getLineSize(); 
    my $chromLines = loadChromLines(); 
    my $regionDNA = getDNASegment($chrom, $start, $end, $lineSize, $chromLines);  
    my $faFile = "$param{refSeqPath}/$targetName.fa";
    open my $faFileH, ">", $faFile or die "could not open $faFile\n";
    print $faFileH ">$targetName\n$regionDNA\n";
    close $faFileH;
}

sub mapRegionRooted {
    my ($inSample, $targetName, $chrom, $start, $end) = @_;
    my $targetSeq = getRegionRefSeq($targetName);
    my $targetLength = length($$targetSeq);
    my $rootedSample = "$inSample\_$param{rooted}";
    my $regionSample = "$inSample\_$targetName";
    getRegionsReads($regionSample, \my%regionReads, $start, $types{Frags}{Normal});
    getRegionsReads($regionSample, \%regionReads,   $start, $types{Frags}{SingleRead});
    getRegionRootReads($rootedSample, \%regionReads, \my%rootReads);
    loadRootSequences(\my%rootSeqs);
    mapRootReadsToRoots($inSample, $targetName, \%rootSeqs, \%rootReads);
    printRegionMap($inSample, $targetName, \%regionReads, \%rootReads, $start, $targetSeq, $targetLength);
}

sub mapRootReadsToRoots {
    my ($inSample, $targetName, $rootSeqs, $rootReads) = @_;
    my $outFile = "$param{inputPath}/$inSample/$inSample\_$targetName\_$param{rooted}.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile\n";
    print $outFileH "GENOME STRAND,GENOME POSITION,PAIRID,ROOT READ,PARTNER READ,";
    foreach my $rootName(sort {$a cmp $b} keys %$rootSeqs){
        foreach my $rootStrand(1..2){
            print $outFileH "$rootName\_$rootStrand,";   
        }
    } 
    print $outFileH "\n";
    foreach my $genomeStrand(sort {$a <=> $b} keys %$rootReads){
        foreach my $pos(sort {$a <=> $b} keys %{$$rootReads{$genomeStrand}}){
            foreach my $rootRead(@{$$rootReads{$genomeStrand}{$pos}}){
                my ($pairID, $rootRead, $partnerRead) = @$rootRead;
                print $outFileH "$genomeStrand,$pos,$pairID,$rootRead,$partnerRead,";
                foreach my $rootName(sort {$a cmp $b} keys %$rootSeqs){
                    foreach my $rootStrand(1..2){
                        print $outFileH mapRootReadToRoot($rootRead, $rootName, $rootSeqs, $rootStrand).",";
                    }   
                }
                print $outFileH "\n";
            }
        }
    }
    close $outFileH;
}

sub mapRootReadToRoot { 
    my ($rootRead, $rootName, $rootSeqs, $strand) = @_;
    $strand == 2 and $rootRead = reverseComplement($rootRead);
    my $readLength = length($rootRead);
    my $targetLength = length($$rootSeqs{$rootName});
    my $maxPos = $targetLength - $readLength;
    my @maps;
    foreach my $pos(0..$maxPos){
        my $target = substr($$rootSeqs{$rootName},$pos,$readLength);
        my $nDiscs = $readLength - (($rootRead ^ $target) =~ tr/\0/\0/);
        $nDiscs <= $param{maxDisc} and push @maps, $pos;
    }
    return (join(":",@maps), scalar(@maps));
}

sub printRegionMap {
    my ($inSample, $targetName, $regionReads, $rootReads, $start, $targetSeq, $targetLength) = @_;
    my $fragType = $types{Frags}{SingleRead};
    my %srLines;
    $srLines{1} = '.' x $targetLength; 
    $srLines{2} = '.' x $targetLength; 
    foreach my $strand(keys %$rootReads){
        my $symbol = ">";
        $strand == 2 and $symbol = "<";
        foreach my $pos(keys %{$$rootReads{$strand}}){
            my $before = $pos - $start - 1;
            $srLines{$strand} =~ m/^(.{$before})(.{$param{readLength}})(.*)/;
            $srLines{$strand} = $1.$symbol x $param{readLength}.$3;
            length($srLines{$strand}) == $targetLength or die "srLine length error\n";
        }
    }
    $fragType = $types{Frags}{Normal};
    my $normalLine = ' ' x $targetLength; 
    foreach my $pairID(keys %{$$regionReads{$fragType}}){
        foreach my $frag(@{$$regionReads{$fragType}{$pairID}}){
            my ($pos1, $strand1, $pos2, $strand2) = @$frag;
            my $before = $pos1 - $start - 1;
            my $length = $pos2 - $pos1 + 1;
            $normalLine =~ m/^(.{$before})(.{$length})(.*)/;
            $normalLine = $1.'-' x $length.$3;
            length($normalLine) == $targetLength or die "normalLine length error\n";
        }
    }
    interleaveLines ($inSample, $targetName, $targetLength, $srLines{1}, $srLines{2}, $$targetSeq, $normalLine);
}

sub interleaveLines {
    my ($inSample, $targetName, $targetLength, @inputs) = @_;
    my $outFile = "$param{inputPath}/$inSample/$inSample\_$targetName\_$param{rooted}.txt";
    open my $outFileH, ">", $outFile or die "could not open $outFile\n";
    my $nLines = int(($targetLength / $param{charPerLine})+1);
    foreach my $i(0..($nLines-1)){
        foreach my $input(@inputs){
            my $line = substr($input,$i * $param{charPerLine}, $param{charPerLine});
            $line =~ m/\S/ or next;
            print $outFileH "$line\n";  
        }
        print $outFileH ' ' x $param{charPerLine}."\n";
    } 
    close $outFileH;
}

sub getRegionsReads { #collect singleReads within region; strands already corrected for isCircles
    my ($regionSample, $regionReads, $start, $fragType) = @_;
    my $fragsTable = getTableName('Frags', $regionSample);
    my $sql = "SELECT trunc(pairid/1E8) pairid, 
                      $start + position1 pos1, strand1,
                      $start + position2 pos2, strand2
               FROM $fragsTable
               WHERE fragmenttype = $fragType
               GROUP BY pairid, position1, strand1, position2, strand2";
    runSQL($sql, \my($pairID, $pos1, $strand1, $pos2, $strand2));
    while(fetchRow()){ 
        push @{$$regionReads{$fragType}{$pairID}}, [$pos1, $strand1, $pos2, $strand2] }
}

sub getRegionRootReads {
    my ($rootedSample, $regionReads, $rootReads) = @_;
    my $srType = $types{Frags}{SingleRead};
    my $pairedReadsFile = "$param{inputPath}/$rootedSample/$rootedSample\_$param{rooted}\_pairedReads.csv"; 
    open my $fileH, "<", $pairedReadsFile or die "could not open $pairedReadsFile\n";
    while(my $line = <$fileH>){
        chomp $line;
        $line or next;
        my ($pairID, $rootRead, $partnerRead) = split(",", $line);
        $rootRead = "\U$rootRead";
        $partnerRead = "\U$partnerRead";
        if ($param{isCircles}){
            $rootRead = reverseComplement($rootRead);
            $partnerRead = reverseComplement($partnerRead);
        }
        $pairID or next;
        $$regionReads{$srType}{$pairID} or next;
        foreach my $regionRead(@{$$regionReads{$srType}{$pairID}}){
            my ($pos1, $strand1, $pos2, $strand2) = @$regionRead;
            push @{$$rootReads{$strand1}{$pos1}}, [$pairID, $rootRead, $partnerRead];
        }
    }
    close $fileH;
}

sub getRegionRefSeq {
    my ($targetName) = @_;
    my $faFile = "$param{refSeqPath}/$targetName.fa";
    open my $faFileH, "<", $faFile or die "could not open $faFile\n";
    my $nameLine = <$faFileH>;
    my $seqLine = <$faFileH>;
    chomp $seqLine;
    return \$seqLine;
}

sub loadRootSequences {
    my ($rootSeqs, $faFile) = @_;
    $faFile or $faFile = "$param{refSeqPath}/$param{rooted}.fa";
    open my $faFileH, "<", $faFile or die "could not open $faFile\n";
    my ($name, $seq) = (undef, undef);
    while(my $line = <$faFileH>){
        chomp $line;
        $line =~ s/\r//g;
        $line or next;
        if ($line =~ m/^>(.*)/){
            if ($name) {
                $$rootSeqs{$name} = "\U$seq";
                ($name, $seq) = (undef, undef);
            }
            $name = $1;
        } else {
            $seq .= $line;
        }
    }
    $name and $$rootSeqs{$name} = "\U$seq";
}

sub mapRegionRooted2 {
    my ($inSample, $targetName, $faFile, $chrom, $start, $end) = @_;
    my $rootedSample = "$inSample\_$param{rooted}";
    my $regionSample = "$inSample\_$targetName";
    getRegionsReads($regionSample, \my%regionReads, $start, $types{Frags}{SingleRead});
    getRegionRootReads($rootedSample, \%regionReads, \my%rootReads);
    loadRootSequences(\my%rootSeqs, "$param{inputPath}/$inSample/$faFile.fa");
    my @rootNames = keys %rootSeqs;
    my $rootName = $rootNames[0];
    my $rootSeq = $rootSeqs{$rootName};
    my $rootLength = length($rootSeq);
    my $fragType = $types{Frags}{SingleRead};
    my %srLines;
    $srLines{1}{symbol} = ' ' x $rootLength; 
    $srLines{2}{symbol} = ' ' x $rootLength; 
    $srLines{1}{seq} = ' ' x $rootLength; 
    $srLines{2}{seq} = ' ' x $rootLength;  
    my %symbols = (1=>'>',2=>'<');
    foreach my $genomeStrand(sort {$a <=> $b} keys %rootReads){
        my $rootStrand = ($genomeStrand % 2) + 1;
        foreach my $pos(sort {$a <=> $b} keys %{$rootReads{$genomeStrand}}){
            foreach my $readPair(@{$rootReads{$genomeStrand}{$pos}}){
                my ($pairID, $rootRead, $partnerRead) = @$readPair;
                my ($rootPos, $rootHits) = mapRootReadToRoot($rootRead, $rootName, \%rootSeqs,  $rootStrand);
                $rootHits == 1 or next;
                my ($partnerPos, $partnerHits) = mapRootReadToRoot($partnerRead, $rootName, \%rootSeqs, $genomeStrand);
                $partnerHits == 1 or next;  
                my $before = $rootPos;
                $srLines{$genomeStrand}{symbol} =~ m/^(.{$before})(.{$param{readLength}})(.*)/;
                $srLines{$genomeStrand}{symbol} = $1.$symbols{$rootStrand} x $param{readLength}.$3;
                $srLines{$genomeStrand}{seq} =~ m/^(.{$before})(.{$param{readLength}})(.*)/;
                $rootStrand == 2 and $rootRead = reverseComplement($rootRead);
                $srLines{$genomeStrand}{seq} = $1.$rootRead.$3;
                $before = $partnerPos;
                $srLines{$genomeStrand}{symbol} =~ m/^(.{$before})(.{$param{readLength}})(.*)/;
                $srLines{$genomeStrand}{symbol} = $1.$symbols{$genomeStrand} x $param{readLength}.$3;
                $srLines{$genomeStrand}{seq} =~ m/^(.{$before})(.{$param{readLength}})(.*)/;
                $genomeStrand == 2 and $partnerRead = reverseComplement($partnerRead);
                $srLines{$genomeStrand}{seq} = $1.$partnerRead.$3;
            }
        }
    } 
    $srLines{1}{symbol} = highlightMismatchedBases($rootLength, $srLines{1}{seq}, $rootSeq, $srLines{1}{symbol});
    $srLines{2}{symbol} = highlightMismatchedBases($rootLength, $srLines{2}{seq}, $rootSeq, $srLines{2}{symbol});
    interleaveLines ($inSample, $rootName, $rootLength, 
                     $srLines{1}{seq}, $srLines{1}{symbol}, $rootSeq, $srLines{2}{symbol}, $srLines{2}{seq});
}

sub highlightMismatchedBases {
    my ($rootLength, $readLine, $seqLine, $symbolLine) = @_;
    foreach my $i(0..($rootLength-1)){
        my $readChar = substr($readLine,$i,1);
        ($readChar eq ' ' or $readChar eq "\n") and next;
        my $seqChar = substr($seqLine,$i,1);
        unless("\U$readChar" eq "\U$seqChar"){
            $symbolLine =~ m/(.{$i})(.)(.*)/;
            $symbolLine = $1.'*'.$3;
        }
    }
    return $symbolLine;
}



1;





