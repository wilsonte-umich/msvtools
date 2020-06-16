#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{xx} or $command{xx} = ['multiThread', '24:00:00', 5000, 0];
defined $param{minPolyALength} or $param{minPolyALength} = 6; 
defined $param{minMappableLength} or $param{minMappableLength} = 20; 
defined $param{padding} or $param{padding} = 1000;

#CAUTION: only run one polyA from an input path at a time!!!
#not handling not_aligned.fa correctly!!
defined $command{mapPolyA} or $command{mapPolyA} = ['multiThread', '12:00:00', 20000, 0];

#known limitation: currently only considering read1
#if ever need to run when unpaired=FALSE
#merge read2.fa into read1.fa and renumber in read1.fa.p
    
sub mapPolyA {
    my ($inSample) = @_;
    my $outSample = $inSample."A"; #force rename by simply appending 'A' to signify 'polyA'
    my $sample = "$inSample\::$outSample"; 
    
    checkPolyAReadLength($inSample);
    
#    status("finding candidate poly A reads\n");
#    prepareReads($sample, \&getPolyAReads_full);
#    
#    status("mapping to purge genome matches\n");
#    purgePolyAGenomeMatches($outSample);
#    
#    status("stripping poly A portion of unmapped reads\n");
#    $param{readType} = 'fasta';
#    prepareReads($outSample, \&getPolyAReads_root);
#    
#    status("remapping non-poly A portion of reads\n");
#    runPolyAMap($outSample);
    
    status("parsing poly A addition sites\n");
    my $hitsTable = newTable('Hits', $outSample);
    parsePolyAHits($outSample, $hitsTable);
#    my $hitsTable = getTableName('Hits', $outSample);
    calcHitStats($outSample);
}
sub checkPolyAReadLength {
    my ($inSample) = @_;
    getSampleReadLength($inSample);
    my $maxGapLength = $param{readLength} - $param{minPolyALength} - $param{minMappableLength};
    $maxGapLength >= 10 or die "readLength=$param{readLength} not long enough to provide useful polyA mapping for ".
                               "minPolyALength=$param{minPolyALength} and minMappableLength=$param{minMappableLength}\n";
}
sub getPolyAReads_full {
    my ($seq) = @_;
    chomp $seq;   
    $seq = "\U$seq";
    $seq =~ m/A{$param{minPolyALength},}$/ and return $seq;
    $seq = reverseComplement($seq);
    $seq =~ m/A{$param{minPolyALength},}$/ and return $seq;
    return undef; #skip all other reads 
}
sub purgePolyAGenomeMatches {
    my ($sample) = @_;
    $param{unpaired} = 1;
    $param{maxHits} = 1; 
    $param{mapType} = 'pass';
    mapReads($sample, $sample, 'read1');    
    getReadFiles($sample, 'purgedFasta', \my%purgedFiles);
    system("mv $purgedFiles{read1} $purgedFiles{read1}.bu");
    getReadFiles($sample, 'fasta', \my%fastaFiles);
    system("mv not_aligned.fa $fastaFiles{read1}");
}
sub getPolyAReads_root {
    my ($seq) = @_;
    chomp $seq;   
    $seq =~ m/(A{$param{minPolyALength},})$/ or return undef;
    my $polyALength = length($1);
    $seq =~ m/(.{$param{minMappableLength},})A{$polyALength}$/ or return undef;
    return $1;
}
sub runPolyAMap {
    my ($sample) = @_;
    $param{mapType} = 'bowtie';
    mapReads($sample, $sample, 'read1');
}
sub parsePolyAHits {
    my ($sample, $hitsTable) = @_;
    my $lineSize = getLineSize(); 
    my $chromLines = loadChromLines(); 
    getMapFiles($sample, \my%mapFiles); 
    open my $mapFileH, "<", $mapFiles{read1} or die "could not open $mapFiles{read1} for reading: $!\n";
    my %hitCount;
    while (my $line = <$mapFileH>){ 
        my ($pairID, $map) = getLine_bowtie(\$line);
        my ($chrom, $pos, $length, $strand, $discs) = split(":", $$map);
        my $polyASite = getPolyASite($chrom, $pos, $length, $strand, $lineSize, $chromLines, $line);
        defined $polyASite or next;
        $hitCount{"$chrom,$polyASite,$strand"}++;
    }
    close $mapFileH; 
    my $outFile = "$hitsTable.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile for writing : $!\n";
    foreach my $hitKey(keys %hitCount){ print $outFileH "$hitKey,$hitCount{$hitKey}\n" }
    close $outFileH;
    loadData($outFile, $hitsTable, ',', "CHROMOSOME,POSITION,STRAND,COUNT_");
}
sub getPolyASite {
    my ($chrom, $pos, $length, $strand, $lineSize, $chromLines, $line) = @_;
    if($strand == 1){
    
        return getPolyASite_top($chrom, $pos, $length, $lineSize, $chromLines, $line);
        
    } else {
        return getPolyASite_bottom($chrom, $pos, $length, $lineSize, $chromLines);
    }
}
sub getPolyASite_top {
    my ($chrom, $pos, $length, $lineSize, $chromLines, $line) = @_;
    my $flankingLength = $param{readLength} - $length;
    
    my $genomic = getDNASegment($chrom, $pos, $pos + $param{readLength} - 1, $lineSize, $chromLines); 
    $genomic = "\U$genomic";
    my @line = split("\t", $line);
    my $readRoot = $line[4];
    my $polyA = 'A' x $flankingLength;
    my $read = $readRoot.$polyA;
    $read = "\U$read";
    my $matches;
    foreach my $i(0..($param{readLength} - 1)){
        my $match = " ";
        substr($genomic,$i,1) eq substr($read,$i,1) and $match = "|";
        $matches .= $match;
    }
    print "+\t$genomic\n\t$matches\nr\t$read\n\n";
    
    


    my $minPos = $pos + $length;
    my $maxPos = $minPos + $flankingLength - 1;
    my $flanking = getDNASegment($chrom, $minPos, $maxPos, $lineSize, $chromLines); 
    $flanking = "\U$flanking";
    
    #yes, collect last non-reference base
    #however, procees entire thing and only return if
    #are at least some number of mismatches in the polyA terminus
    #possibly even consider using disc information in 3' terminal bases of mapped portion?
    #
    #should start by getting a simple printout of ref to read match
    #knowing that the read can easiliy be reconstructed from Discs and adding the missing poly A
    
    for my $i(0..($flankingLength-1)){
        my $base = substr($flanking,$i,1);
        if($base eq 'A'){
            $length++;
        } else {
            return $pos + $length - 1; #the last non-reference base
        }
    }
    return undef;
}
sub getPolyASite_bottom {
    my ($chrom, $pos, $length, $lineSize, $chromLines) = @_;
    my $flankingLength = $param{readLength} - $length;
    my $maxPos = $pos - 1;
    my $minPos = $maxPos - $flankingLength + 1;
    my $flanking = getDNASegment($chrom, $minPos, $maxPos, $lineSize, $chromLines); 
    $flanking = reverseComplement($flanking);
    $flanking = "\U$flanking";
    for my $i(0..($flankingLength-1)){
        my $base = substr($flanking,$i,1);
        if($base eq 'A'){
            $pos--;
        } else {
            return $pos; #the last non-reference base
        }
    }
    return undef;
}

1;

