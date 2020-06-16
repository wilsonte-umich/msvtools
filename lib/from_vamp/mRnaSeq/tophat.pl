#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx;
defined $command{tophat} or $command{tophat} = ['multiThread', '48:00:00', 2000, 0];

sub tophat {
    my ($sample) = @_;
    getDirectories($sample, \my%dirs);
    my $tophatDir = "$dirs{sample}/tophat";
    -d $tophatDir or mkdir $tophatDir;
    getReadFiles($sample, $param{readType}, \my%readFiles);
    my $libType = $param{stranded} ? $param{reverseStrands} ? "fr-secondstrand" : "fr-firststrand" : "fr-unstranded";
    getIsMrna(\my%isMrna);
    my @cmd = ("$param{tophatPath}tophat",
               "--solexa1.3-quals",
               "--max-multihits $param{maxHits}",
               "--initial-read-mismatches $param{maxDisc}",
               "--bowtie-n",
               "--rg-id $sample --rg-sample $sample",
               "--library-type $libType",
               "--output-dir $tophatDir",
               $isMrna{$sample} ?  "--GTF $param{refSeqPath}/$param{refSeq}/refGene_$param{refSeq}.gtf" : "",
               "$param{refSeqPath}/$param{refSeq}/$param{refSeq}",
               "$readFiles{read1}",
               $param{unpaired} ? "" : "$readFiles{read2}" );
    my $cmd = join(" ", @cmd);
    status("running tophat as:\n$cmd\n\n");
    system($cmd);
}

