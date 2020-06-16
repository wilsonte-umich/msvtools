#!/usr/bin/perl
use strict;
use warnings;

my %samFlags = (
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
my %samFields = (
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
    QUAL => 10
);
my ($fa1, $fa2, $baseOffset, $readLength) = @ARGV;
open my $fa1H, ">", $fa1 or die "could not open $fa1\n";
open my $fa2H, ">", $fa2 or die "could not open $fa2\n";
my $counter = 0;
my (%pairs, %committed);
while (my $line = <STDIN>) {
    my @line = split("\t", $line);
    my $name = $line[$samFields{QNAME}];   
    unless($committed{$name}){
        my $flag = $line[$samFields{FLAG}];
        checkSamFlagBit($flag, 'failedQuality') and next;
        checkSamFlagBit($flag, 'isMultiFrag') or next;
        my $readN = 0;
        (checkSamFlagBit($flag, 'isFirst') and $readN = 1) or
        (checkSamFlagBit($flag, 'isLast') and $readN = 2) or
        next;
        my $read = $line[$samFields{SEQ}];
        checkSamFlagBit($flag, 'isReverse') and $read = reverseComplement($read); 
        $read = substr($read, $baseOffset, $readLength);
        $pairs{$name}{$readN} = $read;
        my $otherReadN = ($readN % 2) + 1;
        if ($pairs{$name}{$otherReadN}){
            $counter++;
            print $fa1H ">$counter\n$pairs{$name}{1}\n";
            print $fa2H ">$counter\n$pairs{$name}{2}\n";    
            $committed{$name} = 1;
            delete $pairs{$name};
        }  
    }
}
close $fa1H;
close $fa2H;

sub checkSamFlagBit {
    my ($flag, $flagName) = @_;
    return $flag & $samFlags{$flagName};
}

sub reverseComplement{
    my ($sequence) = @_;
    $sequence = reverse "\U$sequence";
    $sequence =~ tr/ACGT/TGCA/;
    return $sequence;
}

1;

