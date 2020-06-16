#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

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
#1 The read is one of a pair
#2 The alignment is one end of a proper paired-end alignment
#4 The read has no reported alignments
#8 The read is one of a pair and has no reported alignments
#16 The alignment is to the reverse reference strand
#32 The other mate in the paired-end alignment is aligned to the reverse reference strand
#64 The read is the first (#1) mate in a pair
#128 The read is the second (#2) mate in a pair

our %samViewFlags;
#negative composite flags, used in samtools view -F

#%samF holds flag information for use in samtools view
#{+} are bits which must be present
#{-} are bits which cannot be present
our %samF;
$samF{badPair}{+} = $samFlags{failedQuality} +
                      $samFlags{isDuplicate};
$samF{badPair}{-} = 0;
$samF{mappedPair}{+} = $samFlags{isMultiFrag} +
                         $samFlags{isFirst};
$samF{mappedPair}{-} = $samF{badPair}{+} +
                         $samFlags{isUnmapped} +
                         $samFlags{nextIsUnmapped};
                       

$samViewFlags


#positive composite flags, used samtools view -f
$samViewFlags{mappedPair} = $samFlags{isMultiFrag} +
                            $samFlags{isFirst};
$samViewFlags{convergent} = $samViewFlags{mappedPair} +
                            $samFlags{isMultiFrag} +
                            $samFlags{isFirst};


$samFlags{mappedPair}{reqBits} = $samFlags{isMultiFrag} +
                                  $samFlags{isUnmapped} +
                                  $samFlags{nextIsUnmapped} +
                                  $samFlags{isFirst} +
                                  $samFlags{failedQuality} +
                                  $samFlags{isDuplicate};
$samFlags{mappedPair}{setBits} = $samFlags{isMultiFrag} +
                                  $samFlags{isFirst};
                                  
$samFlags{convergent}{reqBits} = $samFlags{mappedPair}{reqBits} +
                                     $samFlags{isReverse} +
                                     $samFlags{nextIsReverse};
$samFlags{convergent}{setBits} = $samFlags{mappedPair}{setBits} +
                                     $samFlags{nextIsReverse};
                                     
my $flag = $samFlags{isMultiFrag} + $samFlags{isFirst};
my $result = checkFlag($flag, 'mappedPair');
$result = checkFlag($flag, 'convergent');
if ($result) {print "passed\n"} else { print "failed\n" }


sub checkFlag {
    my ($flag, $flagName) = @_;
    defined $samFlags{$flagName} or die "$flagName is not a valid SAM flag\n";
    return ($flag & $samFlags{$flagName}{reqBits}) ==  $samFlags{$flagName}{setBits};
}

1;

