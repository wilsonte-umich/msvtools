#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %gffFields %bowtieFields));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $param{samType} or $param{samType} = 'bam'; 
defined $command{getSamGroups} or $command{getSamGroups} = ['multiThread', '1:00:00', 500, 0];

sub getSamGroups { #parse and print the libraries and associated read groups in a bam file
    my ($sample) = @_;
    #@RG     ID:3034_7       PL:SLX  LB:129P2_SLX_500_HC_1   PI:500  DS:129P2_Mouse_Genome   SM:129P2        CN:SC
    getReadFiles($sample, 'sam', \my%samFiles);
    my $command = "$param{samPath}samtools view -H ";
    $command .= "$samFiles{$param{samType}} ";
    $command .= "| perl $param{vampPath}/bin/sam/bin/getSamGroups.pl ";
    my $return = qx/$command/;
    $return and status("$return\n");
}

1;
