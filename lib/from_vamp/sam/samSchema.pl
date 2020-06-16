#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %geneticCode));

our %samFields = (
    QNAME => 0,
    FLAG => 1,
    RNAME => 2,
    POS => 3,
    MAPQ => 4,
    CIGAR => 5,
    RNEXT => 6,
    PNEXT => 7,
    TLEN => 8,
    SEQ => 9,
    QUAL => 10,
    OPT => 11
);
our %samFlags = (
    isMultiFrag => 0x1, #  template having multiple fragments in sequencing
    isAligned => 0x2, #  each fragment properly aligned according to the aligner
    isUnmapped => 0x4, #  fragment unmapped
    nextIsUnmapped => 0x8, #  next fragment in the template unmapped
    isReverse => 0x10, #  SEQ being reverse complemented
    nextIsReverse  => 0x20, #  SEQ of the next fragment in the template being reversed
    isFirst => 0x40, #  the rst fragment in the template
    isLast => 0x80, #  the last fragment in the template
    isSecondary => 0x100, #  secondary alignment
    failedQuality => 0x200, #  not passing quality controls
    isDuplicate => 0x400, #  PCR or optical duplicate
);  

sub getSamFlagBit { #return a boolean telling whether named flag bit is marked
    my ($flag, $flagName) = @_; #flag as value
    return $flag & $samFlags{$flagName};
}




sub setSamFlagBit { #set the named flag bit to the passed boolVal
    my ($flag, $flagName, $boolVal) = @_; #flag as reference
    if ($boolVal) {
        markSamFlagBit($flag, $flagName);
    } else {
        clearSamFlagBit($flag, $flagName);
    }
}
sub markSamFlagBit { #mark the named flag bit as true
    my ($flag, $flagName) = @_; #flag as reference
    $$flag = $$flag | $samFlags{$flagName}; 
}
sub clearSamFlagBit { #mark the named flag bit as false
    my ($flag, $flagName) = @_; #flag as reference
    $$flag = $$flag & ~$samFlags{$flagName}; 
}
sub invertSamFlagBit { #invert the value of the named flag bit
    my ($flag, $flagName) = @_; #flag as reference
    $$flag = $$flag ^ $samFlags{$flagName}; 
}




sub fixBamSexChrom {
    my ($rName) = @_;
    $reverseRefSeqs{$param{refSeq}}{$$rName} =~ m/X/ and $$rName = 'X' and return;
    $reverseRefSeqs{$param{refSeq}}{$$rName} =~ m/Y/ and $$rName = 'Y';
}



sub getSamChrom {
    my ($rName) = @_;
    my $chrom = $refSeqs{$param{refSeqBase}}{$rName};
    defined $chrom and return $chrom;
    $rName = "\U$rName";
    $chrom = $refSeqs{$param{refSeqBase}}{$rName};
    defined $chrom and return $chrom;
    if ($rName =~ m/^CHR(.*)/){
        $chrom = $refSeqs{$param{refSeqBase}}{$1};
        defined $chrom and return $chrom;
    } else {
        $chrom = $refSeqs{$param{refSeqBase}}{"chr$rName"};
        defined $chrom and return $chrom;
        $chrom = $refSeqs{$param{refSeqBase}}{"Chr$rName"};
        defined $chrom and return $chrom;
        $chrom = $refSeqs{$param{refSeqBase}}{"CHR$rName"};
        defined $chrom and return $chrom; 
    }
    return undef;
    #die "$rName is not a recognized chromosome name for refSeq $param{refSeqBase}";
}
sub getSamStrand {
    my ($flag, $flagName) = @_; #$flagName expected to be either 'isReverse' or 'nextIsReverse'
    my $strand = 1;
    getSamFlagBit($flag, $flagName) and $strand = 2;
    return $strand;
}
sub getSamOption {
    my ($optional, $tag_) = @_; #$optional is a reference to an array of split optional fields, each = TAG:TYPE:VALUE
    foreach my $optional(@$optional){
        my ($tag, $type, $value) = split(':', $optional);
        $tag eq $tag_ and return $value;
    }
    return undef;
}

1;

