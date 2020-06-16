#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command));

#callable parameters and commands in this script
defined $param{ftpURL} or $param{ftpURL} = 0;
defined $param{ftpAccount} or $param{ftpAccount} = 0;
defined $param{ftpPassword} or $param{ftpPassword} = 0;
defined $param{runFolder} or $param{runFolder} = 0;
defined $command{downloadIllumina} or $command{downloadIllumina} = ['singleThread', '24:00:00', 1000, 0];

sub downloadIllumina {
    my (@targets) = @_;
    my $fileExt = checkFTPParams(@targets);
    my $site = "$param{ftpURL}/$param{runFolder}";
    my $login = "--ftp-user $param{ftpAccount} --ftp-password $param{ftpPassword}";
    my $recursive = "-r -l 20 -A $fileExt";
    my $wgetCommand = "wget -nv -nd -np $recursive $login $site";
    status("$wgetCommand\n");
    my $wgetResults = qx/$wgetCommand/;
    $wgetResults and status("$wgetResults\n");
    foreach my $target(@targets){
        my ($lane, $sample) = splitSampleName($target);
        getDirectories($sample, \my%dirs);
        -d $dirs{sample} or mkdir $dirs{sample};
        foreach my $readN(1..3){ #3 reads includes possible barcode read
           my $read = "read$readN";
           -d $dirs{$read} or mkdir $dirs{$read};
           my $mvKey = $read;
           if ($param{readType} eq 'qseq'){
               my $qseq = "qseq$readN";
               -d $dirs{$qseq} or mkdir $dirs{$qseq};
               $mvKey = $qseq;
           } 
           #s_8_3_2101_qseq.txt
           #s_8_3_sequence.txt
           my $mvCommand = "mv $lane\_$readN".'*'."$fileExt $dirs{$mvKey}";
           if($param{unpaired} and $param{readType} eq 'solexa_fastq'){
               $readN == 1 or next;
               $mvCommand = "mv $lane".'*'."$fileExt $dirs{$mvKey}";
           } 
           status("$mvCommand\n");
           my $mvResults = qx/$mvCommand/;
           $mvResults and status("$mvResults\n");
        }
    }
}

sub checkFTPParams {
    my (@targets) = @_;
    foreach my $target(@targets){ #check this before downloading...
        my ($lane, $sample) = split('=', $target);
        ($lane and $sample) or die "bad target $target:  must be in format 'lane=sample'\n";
    }
    foreach my $param(qw(ftpURL ftpAccount ftpPassword runFolder)){
        $param{$param} or die "must specify parameter '$param'\n";
    }
    ($param{ftpURL} =~ m/^ftp/ or $param{ftpURL} =~ m/^FTP/) or $param{ftpURL} = "ftp://$param{ftpURL}";
    my $fileExt;
    $param{readType} eq 'solexa_fastq' and $fileExt = 'sequence.txt';
    $param{readType} eq 'qseq' and $fileExt = 'qseq.txt';
    $fileExt or die "parameter 'readType' must be either 'solexa_fastq' or 'qseq'\n";
    return $fileExt;
}

1;



