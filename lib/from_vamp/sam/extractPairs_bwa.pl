#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %fieldNames %refSeqs %reverseRefSeqs %samFlags %samFields));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{extractPairs_bwa} or $command{extractPairs_bwa} = ['multiThread', '24:00:00', 10000, 0];

#extractPairs_bwa streams bwa sampe output directly into vamp pair extraction
#saves disk write of bam file and associated conversions
#however, bam file is of course no longer available for separate manipulation

sub extractPairs_bwa{
    my ($sample) = @_;
    my ($inputSample, $outputSample) = splitSampleName($sample);
    getReadFiles($inputSample, 'solexa_fastq', \my%fastqFiles);
    getMapFiles($outputSample, \my%mapFiles);
    my $fastqFile1 = $fastqFiles{read1};
    my $fastqFile2 = $fastqFiles{read2};
    my $saiFile1 = $mapFiles{read1};
    my $saiFile2 = $mapFiles{read2};
    my $bwaCommand = "$param{bwaPath}bwa sampe -r '\@RG\tID:$outputSample' ";
    $bwaCommand .= "$param{refSeqPath}/$param{refSeq}/$param{refSeq}.fa $saiFile1 $saiFile2 $fastqFile1 $fastqFile2";
    my $pairsTable = newTable('Pairs', $sample);
    my $loadFile = "$pairsTable.csv";
    my $paramPass = catenateParameters_();
    my $timeStamp = getTimeStamp();
    my $vampCommand = "perl $param{vampPath}/bin/sam/bin/extractPairedSam_stream.pl $loadFile $paramPass $timeStamp";

    status("extracting $sample mapped pairs\n");
        my $stream = "$bwaCommand | $vampCommand";
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

1;

