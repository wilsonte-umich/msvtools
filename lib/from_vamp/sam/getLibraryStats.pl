#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %gffFields %bowtieFields));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{bam2Fasta} or $command{bam2Fasta} = ['multiThread', '24:00:00', 10000, 0];

#bam2Fasta actually took 16.5 Gb to process 7 lanes of a 500 bp library
#could reduce memory footprint a bit by trimming redundant parts of read names

sub getLibraryStats { #retrieves reads from bam file and create renamed fasta file
    my ($sample) = @_;    
    getDirectories($sample, \my%dirs); #sample directory must exist and contain $sample.bam
    -d $dirs{read1} or mkdir $dirs{read1};
    -d $dirs{read2} or mkdir $dirs{read2};
    getReadFiles($sample, 'fasta', \my%fastaFiles);
    getReadFiles($sample, 'sam', \my%samFiles);
    my $rgFile = "$dirs{sample}/$sample.rg"; #filter to read groups specific by user in rg file
    my $command = "$param{samPath}samtools view ";
    -e $rgFile and $command .= "-R $rgFile "; 
    $command .= "$samFiles{bam} ";
    $command .= "| perl $param{vampPath}/bin/sam/bin/bam2Fasta.pl ";
    $command .= "$fastaFiles{read1} $fastaFiles{read2} $param{baseOffset} $param{readLength}";
    my $return = qx/$command/;
    $return and status("$return\n");
}

1;
