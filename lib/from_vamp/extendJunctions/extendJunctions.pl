#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %samFields %samFlags));

#callable parameters and commands in this script
defined $param{extReadLength} or $param{extReadLength} = 0;
defined $command{extendJunctions} or $command{extendJunctions} = ['singleThread', '48:00:00', 2000, 0];

my %qseqColumns = (machine=>0,run=>1,lane=>2,tile=>3,x=>4,y=>5,index=>6,read=>7,sequence=>8,quality=>9,filter=>10);
    
sub extendJunctions {
    my ($sample, $sourceSample, $setID) = @_;
    $param{extReadLength} or die "must specify parameter 'extReadLength'\n";
    processSetSidesEJ($sample, $setID, \my%setInfo, \my%reads);
    my @reads = keys %reads;
    scalar(@reads) or die "did not recover any reads for set\n";
#    foreach my $read(@reads){ print "$read\n" } 
    getExtReadsEJ($sourceSample, \@reads, \my%extReads);
    my @extReads = keys %extReads;
    scalar(@extReads) or die "did not recover any extended reads for set\n";
    foreach my $extRead(sort {$a cmp $b} @extReads){ 
        my $rcExtRead = reverseComplement($extRead);
        status("$extRead\t$rcExtRead\t$extReads{$extRead}\n") 
    } 
}
sub processSetSidesEJ {
    my ($sample, $setID, $setInfo, $reads) = @_;
    foreach my $side(1..2) {
        getSetInfoEJ($sample, $setID, $side, $setInfo);
        getExtendJunctionFrags($sample, $setID, $side, $setInfo, \my@frags);
        getExtendedJunctionSeqs($side, $setInfo, \@frags, $reads);
    }
    my $divider = '-' x 100;
    status("setID\t$setID\n");
    status("setType\t$$setInfo{type}\n");
    status("$divider\n");
    foreach my $side(1..2) {
        status("side\t$side\n");
        status("chrom\t$$setInfo{$side}{chrom}\n");
        status("pos\t$$setInfo{$side}{pos}\n");
        status("strand\t$$setInfo{$side}{strand}\n");
        status("$$setInfo{$side}{refExtRead}\n");
        status("$divider\n");
    }
}
sub getSetInfoEJ {
    my ($sample, $setID, $side, $setInfo) = @_;
    my $setsTable = getTableName('Sets', $sample);
    my $chromSQL = getChromSQL($side);
    my $posSQL = "decode(strand1,1,overlapstart,spanstart)";
    $side == 2 and $posSQL = "decode(strand2,2,overlapend,spanend)";
    my $setsSQL = "SELECT settype, $chromSQL chromosome, $posSQL position, strand$side strand
                   FROM $setsTable
                   WHERE setid = $setID";     
    runSQL($setsSQL, \my($type, $chrom, $pos, $strand));     
    fetchRow();
    $type or die "could not find set $setID in $setsTable\n";
    $$setInfo{type} = $type;
    $$setInfo{$side}{chrom} = $chrom;
    $$setInfo{$side}{pos} = $pos;
    $$setInfo{$side}{strand} = $strand;
}
sub getChromSQL {
    my ($side) = @_;
    my $chromSQL = "chromosome1";
    $side == 2 and $chromSQL = "decode(chromosome2,0,chromosome1,chromosome2)";
    return $chromSQL;
}
sub getExtendJunctionFrags {
    my ($sample, $setID, $side, $setInfo, $frags) = @_;
    my $fragsTable = getTableName('Frags', $sample);
    my $chromSQL = getChromSQL($side);
    my $fragSQL = "SELECT length$side length, (discrepancies$side || 'x') discs
                   FROM $fragsTable
                   WHERE fragmenttype = $$setInfo{type}
                     AND $chromSQL = $$setInfo{$side}{chrom}
                     AND position$side = $$setInfo{$side}{pos}
                     AND strand$side = $$setInfo{$side}{strand}";
    runSQL($fragSQL, \my($length, $discs));
    while (fetchRow()) { push @$frags, [$length, $discs] }
    scalar(@$frags) or die "could not find frag corresponding to set side $side\n";
}
sub getExtendedJunctionSeqs {
    my ($side, $setInfo, $frags, $reads) = @_;
    our %discTypes = (1=>'A', 2=>'C', 3=>'G', 4=>'T', 5=>'N'); 
    foreach my $frag(@$frags){
        my ($readlength, $discs) = @$frag;
        my ($extLength, $extStart) = ($param{extReadLength} - $readlength, $$setInfo{$side}{pos} + $readlength);
        $$setInfo{$side}{strand} == 2 and $extStart -= $param{extReadLength};
        my $refRead = getDNASegment_SeqFile($$setInfo{$side}{chrom}, $$setInfo{$side}{pos}, undef, $readlength);
        my $refExt =  getDNASegment_SeqFile($$setInfo{$side}{chrom}, $extStart,  undef, $extLength);
        $$setInfo{$side}{refExtRead} = $refRead.$refExt;
        $$setInfo{$side}{strand} == 2 and $$setInfo{$side}{refExtRead} = $refExt.$refRead;
        my $read = $refRead;
        $discs =~ m/(.*)x$/; #strip the 'x' appended to force Perl to use $discreps as string
        while ($1){
            #$1 = discreps remaining; $2 = relPos, $3 = discType, $4 = extra
            $1 =~ m/(.*)(..)(.)(.)$/ or $1 =~ m/(.*)(.)(.)(.)$/;
            if($3){
                my $leading = $2 - 1;
                my $newBase = $discTypes{$3};
                defined $newBase or die "dang, I encountered an indel and I'm stuck, heh heh\n";
                $read = substr($read,0,$leading).$newBase.substr($read,$leading+1);   
            }
        }
        $read =~ s/N//g;
        $$setInfo{$side}{strand} == 2 and $read = reverseComplement($read);
        $$reads{$read}++;
    }
}
sub getDNASegment_SeqFile{
    my ($chrom, $start, $end, $length) = @_; #1-referenced positions, must provide either end or length
    $start--;
    $end and $end--;
    $length or $length = $end - $start + 1;
    my $file = "$param{refSeqPath}/$param{refSeq}/$reverseRefSeqs{$param{refSeq}}{$chrom}.seq";
    open my $fileH, "<", $file or die "could not open $file: $!\n";   
    seek($fileH, $start, 0);
    read($fileH, my $sequence, $length);
    close $fileH;    
    return "\U$sequence";
}
sub getExtReadsEJ {
    my ($sample, $reads, $extReads) = @_;
    getDirectories($sample, \my%dirs);
    
    
    foreach my $readN(1..2) {
    #my $readN = 1;
    
        foreach my $tileFile (<$dirs{"qseq$readN"}/*qseq.txt>){ #read files from qseq subdir
        #foreach my $tileFile ("/home/wilsonte/data/yeastCNV/yeastMMS/1/qseq/s_8_1_1101_qseq.txt"){ #read files from qseq subdir
        
        
            if ($tileFile =~ m/.*\/(s_\d)_\d_(\d\d\d\d)_qseq.txt$/){
                open my $inH, "<", $tileFile or die "could not open $tileFile: $!\n";
                while (<$inH>){
                    chomp $_;
                    my @in = split(/\t/, $_);
                    if ($in[$qseqColumns{filter}]){ 
                        my $shortRead = getUserSequence($in[$qseqColumns{sequence}]);
                        foreach my $read(@$reads){ 
                            if($shortRead =~ m/$read/){
                                my $readID = join(":", @in[$qseqColumns{machine}..$qseqColumns{y}]);
                                $readID .= "/$readN";
                                my $extRead = $in[$qseqColumns{sequence}];
                                $$extReads{$extRead} .= "$readID,";
                                last;
                            }  
                        }
                    }
                }
                close $inH;
            }
       }
       
       
    } 
    
    
}

1;

