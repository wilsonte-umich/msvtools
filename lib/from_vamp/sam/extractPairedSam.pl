#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %fieldNames %refSeqs %reverseRefSeqs %samFlags %samFields));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{extractPairedSam} or $command{extractPairedSam} = ['multiThread', '24:00:00', 10000, 0];
defined $command{extractPairedBam} or $command{extractPairedBam} = ['multiThread', '24:00:00', 10000, 0];
defined $command{extractPairedSamFile} or $command{extractPairedSamFile} = ['singleThread', '24:00:00', 10000, 0];
defined $command{extractPairedBamFile} or $command{extractPairedBamFile} = ['singleThread', '24:00:00', 10000, 0];

sub extractPairedSam {
    my ($sample) = @_;
    my ($inputSample, $outputSample) = splitSampleName($sample);
    getMapFiles($inputSample, \my%mapFiles);
    extractPairedSamFile($outputSample, $mapFiles{sample});
}
sub extractPairedBam {
    my ($sample) = @_;
    my ($inputSample, $outputSample) = splitSampleName($sample);
    getMapFiles($inputSample, \my%mapFiles);
    extractPairedBamFile($outputSample, $mapFiles{sample});
}
sub extractPairedBamFile { 
    my ($sample, $bamFile) = @_;
    extractPairedSamFile($sample, $bamFile, 1);
}
sub extractPairedSamFile { 
    my ($sample, $samFile, $isBam) = @_;
    my $isSamFlag = getIsSamFlag($isBam);
    my $samCommand = "$param{samPath}samtools view $isSamFlag $samFile";
    unless (checkBamSort($samFile, $isSamFlag)){
        $isBam or $samFile = sam2Bam($samFile); #samtools only sorts bam apparently
        $samCommand = "$param{samPath}samtools sort -no -m 10000000000 $samFile $samFile.queryName | ";
        $samCommand .= "$param{samPath}samtools view -";
    }
    my $pairsTable = newTable('Pairs', $sample); 
    my $loadFile = "$pairsTable.csv";
    my $paramPass = catenateParameters_();  
    my $timeStamp = getTimeStamp();
    my $vampCommand = "perl $param{vampPath}/bin/sam/bin/extractPairedSam_stream.pl $loadFile $paramPass $timeStamp";

    status("extracting mapped pairs from $samFile\n");
        my $stream = "$samCommand | $vampCommand";
        system($stream);
    
    status("uploading pairs and discrepancies into $pairsTable\n");
        loadData($loadFile, $pairsTable, ",", $fieldNames{Pairs});

    unless($param{extractDataOnly}){
        status("creating pair size histogram...\n");
            my $histogramTable = newTable('h', $pairsTable);
            createPairSizeHistogram($histogramTable, $pairsTable);

        status("calculating pair statistics...\n");
            my $statsTable = newTable('Stats', $sample);    
            calculateNormalStats($histogramTable, $statsTable);
    }
    
    return $pairsTable;
}

sub getIsSamFlag {
    my ($isBam) = @_;
    my $isSamFlag = "-S";
    $isBam and $isSamFlag = "";
    return $isSamFlag;
}
sub checkBamSort { #returns true only if the bam is sorted by query name, as required for pair extraction
    my ($samFile, $isSamFlag) = @_;
    status("checking if $samFile is sorted by query name\n");
    my $samtoolsCommand = "$param{samPath}samtools view $isSamFlag -H $samFile";
    my $header = qx/$samtoolsCommand/;
    if($header =~ m/SO:/){
        #sort specified, must be sorted by queryname
        return ($header =~ m/SO:queryname/);
    } else {
        #sort unspecified, assume file is mapping output sorted by read pair
        return 1; 
    }
}
sub sam2Bam { 
    my ($samFile) = @_;
    my $bamFile = "$samFile.bam";
    status("converting $samFile to $bamFile\n");
    my $samtoolsCommand = "$param{samPath}samtools view -Shb -o $bamFile $samFile";
    my $header = qx/$samtoolsCommand/;
    return $bamFile;
}

sub catenateParameters_{
    my @joinedParams;
    foreach my $param (keys %param) {push @joinedParams, join("#", ($param, $param{$param})) }
    my $paramPass = join("##", @joinedParams);
    return $paramPass;
}

1;

