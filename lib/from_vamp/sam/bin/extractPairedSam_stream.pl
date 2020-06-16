#!/usr/bin/perl
use strict;
use warnings;

my ($loadFile, $paramPass, $timeStamp) = @ARGV; #passed by the stream caller

#recover vamp's operating parameters
our %param = ();
my @paramJoins = split("##", $paramPass);
foreach my $paramJoin (@paramJoins) {
    my ($param, $value) = split("#", $paramJoin);
    $param{$param} = $value;
}

#initialize schema and nested requires
our (%types, %fields, %fieldNames, %partitions, %refSeqs, $uid, %archiveTables, %samFlags);
require "$param{vampPath}/bin/Schema.pl"; 
require "$param{vampPath}/bin/sam/samSchema.pl"; 

#loop through the stream input
open my $loadFileH, ">", $loadFile or die "could not open $loadFile: $!\n";
my ($nPairs, $prevQName, %maps) = (0, 0);
while (my $line = <STDIN>) {
    $line =~ m/^\@/ and next; #skip header lines
    my ($qName,$flag,$rName,$pos,$mapQ,$cigar,$rNext,$pNext,$tLen,$seq,$qual,@optional) = split("\t", $line);
    getSamFlagBit($flag, 'failedQuality') and next; 
    getSamFlagBit($flag, 'isUnmapped') and next; #vamp does not store unmapped read information
    $seq eq "*" and next; #best practice if alignment deemed secondary by mapper; vamp obeys mapper in this case
    $qName =~ m/(.+)\/\d$/ and $qName = $1;  #strip read identifiers (i.e. trailing /1 etc.) if present
    $qName eq $prevQName or commitSamPairs();
    my $chrom = getSamChrom($rName); 
    $chrom or next;
    my $length = length($seq);
    my $strand = getSamStrand($flag, 'isReverse');
    if ($strand == 2){
        $seq = reverseComplement($seq);
    } else {
        $seq = "\U$seq";
    }
    
    ##############################################
    #replace this with actual vamp-format extraction of the SAM discrepancy information
    #just comment these lines out, don't delete them
    my $discs = getSamOption(\@optional, 'NM');
    $discs or $discs = 0;
    ##############################################
    
    push @{$maps{$seq}}, [$chrom, $pos, $length, $strand, $discs];
    $prevQName = $qName;
}
close $loadFileH;

sub commitSamPairs { 
    my $pairID = ($nPairs * 1E8) + $timeStamp;
    my ($read1, $read2) = keys %maps;
    if($read1){
        if($read2){
            commitPairedMaps($loadFileH, $pairID, $maps{$read1}, $maps{$read2});
        } else {
            commitSingletonMaps($loadFileH, $pairID, $maps{$read1});
        }   
    }
    %maps = ();
    $nPairs++;
}

1;

