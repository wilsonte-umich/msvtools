#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#ExtractPairs_PASS.pl contains code specific to reading PASS GFF-format output files.
#########################################################################

use vars(qw(%param %types %fields %refSeqs));

our %gffFields = ( #gff format
    chromosome => 0,
    source => 1,
    feature => 2,
    start => 3,
    end => 4,
    score => 5,
    strand => 6,
    frame => 7,
    attributes => 8);
    
my %strands = ('+' => 1, '-' => 2);

sub getLine_pass{
    my ($line)= @_;
    chomp $$line;
    my @in = split(/\t/, $$line);
    unless ($in[$gffFields{attributes}] =~ m/Name=(.+?);/){return 0}
    my $pairID = $1;
    my $chrom = $refSeqs{$param{refSeqBase}}{$in[$gffFields{chromosome}]};
    my $strand = $strands{$in[$gffFields{strand}]};
    $param{isCircles} and $strand = ($strand % 2) + 1;
    my ($before, $discs) = getDiscrepancies_pass($in[$gffFields{attributes}]);
    my $pos = $in[$gffFields{start}] - $before;
    if ($pos < 0){$pos = 0}
    my $length;
    if ($param{expandEnds} and $param{readLength}){
        $length = $param{readLength};
    } else {
        $length = $in[$gffFields{end}] - $in[$gffFields{start}] + 1;
    }
    my $map = join(":", ($chrom, $pos, $length, $strand, $discs));
    return ($pairID, \$map);        
}

#genome reference position = start(=position) + relativePosition
sub getDiscrepancies_pass{
    my ($attributes) = @_;
    my ($before, $after, $discs, $offset) = (0, 0, 0, 0);
    my %encountered = ();
    #account for the mismatched ends of reads that pass doesn't report, will record these as ambiguous
    if ($param{expandEnds} and $param{readLength}
        and $attributes =~ m/P="(\d+)-(\d+?)";/){ #P="35-2";Note="M:0,G:0";Hits=1;
        if ($1 < $2) {  #calculate end offsets, i.e. how many read bases pass didn't report on
            $before = $1 - 1;
            $after = $param{readLength} - $2;
        } else {
            $before = $param{readLength} - $1;
            $after = $2 - 1;
        }  
    }   
    #collect the mismatches and gaps
    unless($param{noDisc}){        
        if ($attributes =~ m/Note="(.+?)";/) { #do nothing if -info_gff was not used during pass execution
            my ($M, $G) = split(",", $1); #M:0,G:0  or  M:2 -> 3/3 ?/T 23/24 T/G,G:1 -> Q18/G
            if ($M =~ m/-> (.+)/){ #M:0   or   M:2 -> 3/3 ?/T 23/24 T/G
                my @mismatches = split(" ", $1); #3/3 ?/T 23/24 T/G
                for (my $i=0; $i<=$#mismatches; $i+=2){ #23/24 T/G = queryRelPos/refRelPos and queryValue/refValue
                    my ($queryRelPos, $refRelPos) = split("/", $mismatches[$i]);
                    my ($queryValue) = split ("/", $mismatches[$i+1]);
                    my $discType = ($types{Discs}{$queryValue} or $types{Discs}{'?'}); #convert any unknown sequence characters to ? 
                    $discs = addDiscrepancy($discs, ($refRelPos + 1 + $before), $discType, $types{Discs}{None}, \%encountered);
                }
            }
            if ($G =~ m/-> (.+)/){ #G:0   or   G:2 -> Q18/G R40/A/T
                my @gaps = split(" ", $1); #Q18/G R40/A/T
                foreach my $gap (@gaps){ #Q18/G = Query gap, relPos/refValue
                    if ($gap =~ m/^Q(\d+)/){#18
                        $discs = addDiscrepancy($discs, ($1 + 1 + $before), $types{Discs}{Missing}, $types{Discs}{None}, \%encountered);                        
                        $offset++;
                    } elsif ($gap =~ m|^R(\d+)/(.+)/|){ #R40/A/T = Reference gap (i.e. Query insertion), relPos/queryValue/nextRefValue(i.e. after gap)
                        $discs = addDiscrepancy($discs, ($1 - 1 + 1 + $before), $types{Discs}{Extra}, $types{Discs}{$2}, \%encountered);
                        $offset--;
                    }
                }
            }
        }
        #enter padded ambiguous bases at ends of reads is using expand ends
        #adjust right-end relative positions by the number of prior Missing and Extra as needed
        if ($before > 0){
            for my $i(1..$before){    
                $discs = addDiscrepancy($discs, $i, $types{Discs}{'?'}, $types{Discs}{None}, \%encountered);
            }
        }
        if ($after > 0){
            for my $i(($param{readLength} - ($after - 1))..$param{readLength}){
                $discs = addDiscrepancy($discs, $i + $offset, $types{Discs}{'?'}, $types{Discs}{None}, \%encountered);
            }
        }
    }   
    return ($before, $discs);
}

1;

