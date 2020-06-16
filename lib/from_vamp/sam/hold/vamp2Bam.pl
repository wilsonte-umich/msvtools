#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %gffFields %bowtieFields));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{vamp2Bam} or $command{vamp2Bam} = ['multiThread', '24:00:00', 10000, 0];


sub vamp2Bam { #converts vamp Oracle Frags table to SAM or BAM
    my ($sample) = @_;
    my $fragsTable = getTableName('Frags', $sample);
    
    my $fragsSQL = "SELECT trunc(PAIRID/1E8) QNAME, (PAIRID/1E8 - trunc(PAIRID/1E8))*1E8 readGroup,
                    CHROMOSOME1 RNAME, CHROMOSOME2 RNEXT,
                    POSITION1 pos, POSITION2 PNEXT,
                    LENGTH1 || 'M' CIGAR, 
                    POSITION2 - POSITION1 TLEN,

                           decode(STRAND1,1,0,1) isReverse,
                           decode(STRAND2,1,0,1) nextIsReverse

                    FROM $fragsTable
                    WHERE FRAGMENTTYPE = $types{Frags}{Normal}
                       OR (FRAGMENTTYPE > $types{Frags}{Normal} AND FRAGMENTTYPE <= $types{Frags}{DiffChrom} AND NSETSPAIR = 1 AND NSETSFRAG = 1)";
                  
                  

#Normal=>0, Deletion=>1, Insertion=>2, Inversion=>4, Duplication=>8,
#                    DiffChrom=>16,        
#                 


#                    FRAGMENTSIZE NUMBER, EVENTSIZE NUMBER, STDEVNORMAL NUMBER, ENDTOLERANCE NUMBER, 


#our %samFlags = (
#    isMultiFrag => 0x1, #  template having multiple fragments in sequencing
#    isAligned => 0x2, #  each fragment properly aligned according to the aligner
#    isUnmapped => 0x4, #  fragment unmapped
#    nextIsUnmapped => 0x8, #  next fragment in the template unmapped
#    isReverse => 0x10, #  SEQ being reverse complemented
#    nextIsReverse  => 0x20, #  SEQ of the next fragment in the template being reversed
#    isFirst => 0x40, #  the rst fragment in the template
#    isLast => 0x80, #  the last fragment in the template
#    isSecondary => 0x100, #  secondary alignment
#    failedQuality => 0x200, #  not passing quality controls
#    isDuplicate => 0x400, #  PCR or optical duplicate
#);  
#    
#our %samFields = (
#    FLAG => 1,
#    RNAME => 2,
#    POS => 3,
#    MAPQ => 4,
#    CIGAR => 5,
#    RNEXT => 6,
#    PNEXT => 7,
#    TLEN => 8,
#    SEQ => 9,
#    QUAL => 10
#);


#                          
#    getDirectories($sample, \my%dirs); #sample directory must exist and contain $sample.bam
#    -d $dirs{read1} or mkdir $dirs{read1};
#    -d $dirs{read2} or mkdir $dirs{read2};
#    getReadFiles($sample, 'fasta', \my%fastaFiles);
#    getReadFiles($sample, 'sam', \my%samFiles);
#    my $rgFile = "$dirs{sample}/$sample.rg"; #filter to read groups specific by user in rg file
#    my $command = "$param{samPath}/samtools view ";
#    -e $rgFile and $command .= "-R $rgFile "; 
#    $command .= "$samFiles{bam} ";
#    $command .= "| perl $param{vampPath}/bin/sam/bin/bam2Fasta.pl ";
#    $command .= "$fastaFiles{read1} $fastaFiles{read2} $param{baseOffset} $param{readLength}";
#    my $return = qx/$command/;
#    $return and status("$return\n");
}

1;
