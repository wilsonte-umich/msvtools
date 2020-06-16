#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %samFlags %samFields));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $param{bamNormalOnly} or $param{bamNormalOnly} = 0; 
defined $command{vamp2Bam} or $command{vamp2Bam} = ['multiThread', '24:00:00', 2000, 0];

sub vamp2Bam { #converts vamp Oracle Frags table to SAM or BAM
    my ($sample) = @_;
    my $fragsTable = getTableName('Frags', $sample);    
    getDirectories($sample, \my%dirs);  
    my $outRoot = "$dirs{sample}/$fragsTable";  
    my $unsortedBam = "$outRoot.unsorted.bam"; 
    my $sortedBam = "$outRoot.bam"; 
    my $maxChrom = nChrom(); 
    #my $pipe = "perl $param{vampPath}/bin/sam/bin/vamp2Bam.pl $fragsTable $param{bamNormalOnly} $param{refSeq} $param{vampPath} $maxChrom $param{dbLogin}";
    my $pipe = "perl $param{vampPath}/bin/sam/bin/vamp2Bam.pl $fragsTable $param{bamNormalOnly} $param{refSeq} $param{vampPath} $maxChrom $param{dbLogin} | ";
    $pipe .= "$param{samPath}samtools view -S -b -h -o $unsortedBam -";
    system($pipe);
    system("$param{samPath}samtools sort -m 1000000000 $unsortedBam $outRoot");
    unlink $unsortedBam;
    system("$param{samPath}samtools index $sortedBam");
}



1;
