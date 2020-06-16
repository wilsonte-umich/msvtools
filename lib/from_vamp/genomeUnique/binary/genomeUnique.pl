#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#binSize is required but already defined by VAMP
#defined $param{upstreamBins} or $param{upstreamBins} = 5; #number of bins preceding the gene start
defined $command{genomeUnique} or $command{genomeUnique} = ['singleThread', '24:00:00', 10000, 0];

my (%segBinFileHs, $segDir);

sub genomeUnique {

    status("initializing variables...\n"); 
    $segDir = "$param{refSeqPath}/$param{refSeq}";
    -d $segDir or mkdir $segDir;
    $segDir .= '/segs';
    -d $segDir or mkdir $segDir;
    my %binVals =   (A=>'00',C=>'01',G=>'10',T=>'11');
    my %rcBinVals = (A=>'11',C=>'10',G=>'01',T=>'00');
    my ($bitsPerBase, $bitsPerByte) = (2,8);
    my ($chromByteLength, $chromBlockLength) = (1, 8);
    my ($posByteLength, $posBlockLength) = (4, 32);
    my $seqBitLength = $param{readLength} * $bitsPerBase;
    my $seqByteLength = int(($seqBitLength/$bitsPerByte) + 0.99); #in bytes
    my $seqBlockLength = $seqByteLength * $bitsPerByte; #in bits
    my $increment = $param{readLength} - 1;
    my $dupBinSize = 2 * $param{dupBinSize};    

    status("loading $param{refSeq} chromosome sequences...\n"); 
    my $lineSize = getLineSize(); 
    my $chromLines = loadChromLines();     

    status("writing genome x-mers to tmp files...\n"); #subs found in findHomozygous.pl
    foreach my $chrom(1..nChrom()){
        status("  chrom $chrom");
        my $binChrom = pack("C", $chrom); #one byte, 8-bits, for chrom
        my $pos = 0; 
        my (@bases, @rcBases);
        foreach my $zrLine(sort {$a <=> $b} keys %{$$chromLines{$chrom}}){
            my $maxLineI = length($$chromLines{$chrom}{$zrLine}) - 1;
            foreach my $i(0..$maxLineI){
                my $base = substr($$chromLines{$chrom}{$zrLine}, $i, 1);
                $base = "\U$base"; #force uppercase
                push @bases, $binVals{$base}; 
                unshift @rcBases, $rcBinVals{$base};   
                scalar(@bases) == $param{readLength} or next; 
                $pos++;
                unless(grep {!defined($_)} @bases){ #any N causes segment to be marked as bad (since it can never be marked as good) 
                    my $bitString = join('', @bases);
                    my $rcBitString = join('', @rcBases);
                    ($bitString) = sort {$a cmp $b} ($bitString, $rcBitString); #take seq as sorted value of seq and its reverse complement
                    my $binSeq = substr($bitString, 0, $dupBinSize);   
                    my $binFileH = getSegBinFileH($binSeq);  
                    print $binFileH pack("B$seqBlockLength", $bitString); #seq fills as many bytes as dictated by readLength
                    print $binFileH $binChrom; #chrom fills the next byte, 8 bits
                    print $binFileH pack("L", $pos); #pos fills the next 4 bytes, 32 bits
                }  
                shift @bases;         
                pop @rcBases;                        
            }  
        }
        status(", maxPos $pos\n");
    }
    closeSegBinFileHs();
    %$chromLines = ();
    $chromLines = undef;
    
    status("writing uniq position file\n");   
    my $uniqPosFile = "$segDir/$param{refSeq}_uniqPos.csv";
    open my $uniqPosH, ">", $uniqPosFile or die "could not open $uniqPosFile\n";
    foreach my $binSeq (keys %segBinFileHs) { #purge duplicate pairs within each binned temporary read file
        my (%good, %bad);
        my $binFile = getSegBinFile($binSeq);
        open my $binFileH, "<:raw", $binFile or die "could not open $binFile\n";
        my ($seq, $chrom, $pos);
        while (read $binFileH, $seq, $seqByteLength) {
            read $binFileH, $chrom, $chromByteLength;
            read $binFileH, $pos, $posByteLength;
            $bad{$seq} and next;            
            if($good{$seq}){
                $bad{$seq} = 1;
                delete $good{$seq};
            } else {
                $good{$seq} = [$chrom,$pos];
            }
        }            
        close $binFileH;
        #unlink $binFile;
        foreach my $seq(keys %good){
            my ($chrom, $pos) = @{$good{$seq}};
            $chrom = unpack("C", $chrom);
            $pos = unpack("L", $pos);
            print $uniqPosH "$chrom,$pos\n";          
        }        
    }    
    close $uniqPosH;
    
}


sub getSegBinFileH {
    my ($binSeq) = @_;
    unless($segBinFileHs{$binSeq}){
        my $binFile = getSegBinFile($binSeq);
        open my $binFileH, ">:raw", $binFile or die "could not open $binFile\n"; #binary file output
        $segBinFileHs{$binSeq} = $binFileH;
    }
    return $segBinFileHs{$binSeq};
}
sub getSegBinFile {
    my ($binSeq) = @_;
    return "$segDir/$param{refSeq}_$binSeq.txt";
}
sub closeSegBinFileHs {
    foreach my $binSeq (keys %segBinFileHs) {
        my $binFileH =  $segBinFileHs{$binSeq};
        close $binFileH;
    }
}

1;


