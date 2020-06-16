#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#ExtractPairs_Bowtie.pl contains code specific to reading standard Bowtie output files.
#########################################################################

use vars(qw(%param %types %fields %refSeqs));

our %bowtieFields = ( #default Bowtie output format
    name => 0,
    strand => 1,
    chromosome => 2,
    position => 3,
    sequence => 4,
    qualities => 5,
    nOtherInstances => 6,
    mismatches => 7);
    
my %strands = ('+' => 1, '-' => 2);

sub getLine_bowtie{
    my ($line)= @_;
    chomp $$line;
    my @in = split(/\t/, $$line);
    my $pairID = $in[$bowtieFields{name}];    
    my $chrom = $refSeqs{$param{refSeqBase}}{$in[$bowtieFields{chromosome}]};
    my $strand = $strands{$in[$bowtieFields{strand}]};
    $param{isCircles} and $strand = ($strand % 2) + 1;    
    my $pos = $in[$bowtieFields{position}];      
    if ($pos < 0){$pos = 0}
    my ($discs) = getDiscrepancies_bowtie($in[$bowtieFields{mismatches}]);
    my $length = length($in[$bowtieFields{sequence}]);
    my $map = join(":", ($chrom, $pos, $length, $strand, $discs));
    return ($pairID, \$map);        
}

sub getDiscrepancies_bowtie{
    #Comma-separated list of mismatch descriptors. 
    #If there are no mismatches in the alignment, this field is empty. 
    #A single descriptor has the format offset:reference-base>read-base. 
    #The offset is expressed as a 0-based offset from the high-quality (5') end of the read. 
    my ($mismatches) = @_;
    $mismatches or return 0;
    my $discs = 0;
    my %encountered;
    my @mismatches = split(",", $mismatches);
    foreach my $mismatch(@mismatches){
        $mismatch =~ m/^(.*):(.)>(.)$/;
        $discs = addDiscrepancy($discs, $1 + 1, $types{Discs}{$3}, $types{Discs}{None}, \%encountered)
    }
    return $discs;
}

1;

