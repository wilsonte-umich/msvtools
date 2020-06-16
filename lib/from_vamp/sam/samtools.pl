#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %samFields %samFlags));

my $maxTries = 10;

#note that samtoolsView does NOT stream data, but rather loads all recovered reads into memory
#therefore NOT suitable for large downloads, only small target regions!
#use bam2Fasta stream for large downloads
sub samtoolsView { #iteratively attempt to download target regions from a bam file (allowing for failed FTP...)
    my ($bamPath, $regions, $samResults) = @_; #regions is an arrayref of arrayrefs, samResults is an arrayref where bam line data are put
    getRegionList($regions, \my$regionList);
    my $samCommand = "$param{samPath}samtools view -h $bamPath $regionList";
    my ($nTries, $success, @samResults) = (0, 0);
    while(!$success and $nTries < $maxTries){ #iteratively attempt to download data
        my $stdout = qx/$samCommand/;
        @$samResults = split("\n", $stdout); #place retrieved data into line array
        $success = $$samResults[0] =~ m/^\@HD/; #verify download by presence of header
        $nTries++;
        $success or sleep(1); #give FTP a second to recover if failed
    }
    $success or die "failed to retrieve bam data: $samCommand\n";  
}

sub getRegionList {
    my ($regions, $regionList, $useChrPrefix) = @_;
    $$regionList = "";
    if($regions){ #extract passed target regions
        foreach my $region(@$regions){
           my ($chrom, $start, $end) = @$region;
           $useChrPrefix and $chrom = "chr$chrom";
           $$regionList .= "$chrom:$start-$end "; 
        }           
    }  
}

sub sam2FastaArray { #takes sam line array output from previous sub and dumps to fasta files
    my ($samResults, $writeType, $readFiles, $pairs, $pairID) = @_; #$samResults as above, $writeType = ">" or ">>"
    open my $fa1H, $writeType, $$readFiles{read1} or die "could not open $$readFiles{read1}\n"; #open passed target read files
    open my $fa2H, $writeType, $$readFiles{read2} or die "could not open $$readFiles{read2}\n"; 
    my %committed;
    foreach my $line(@$samResults){ #samResults is an array of lines
        $line =~ m/^\@/ and next; #ignore header lines
        my @line = split("\t", $line);    
        my $name = $line[$samFields{QNAME}];  
        unless($committed{$name}){
            my $flag = $line[$samFields{FLAG}];
            getSamFlagBit($flag, 'failedQuality') and next;
            my $read = $line[$samFields{SEQ}];
            my $readLength = length($read); #library read length, as sequenced
            getSamFlagBit($flag, 'isReverse') and $read = reverseComplement($read); 
            $read = substr($read, 0, $param{readLength}); #take only desired portion of read
            my ($rName,$pos,$rNext,$pNext) = ($line[$samFields{RNAME}], $line[$samFields{POS}], $line[$samFields{RNEXT}], $line[$samFields{PNEXT}]);
            my $nextIsReverse = getSamFlagBit($flag, 'nextIsReverse');
            $$pairs{$name}{$read} = [$rName,$pos,$rNext,$pNext,$nextIsReverse,$readLength];
            if (scalar(keys %{$$pairs{$name}}) > 1){ #print to file when both reads are known
                $$pairID++;
                my ($read1, $read2) = keys %{$$pairs{$name}};
                print $fa1H ">$$pairID\n$read1\n";
                print $fa2H ">$$pairID\n$read2\n";    
                $committed{$name} = 1; #keep track of which read pairs are completely known
                delete $$pairs{$name};  #%pairs hash will end up with pairs for which only one read was encountered
            }  
        }
    }
    close $fa1H;
    close $fa2H;         
}

sub getPartnerRegions { #collect the partner regions for pairs from above sub where only one read was encountered
    my ($pairs, $regions, $names) = @_; 
    foreach my $name(keys %$pairs){
        my ($read) = keys %{$$pairs{$name}}; #pairs hash contains only one read per pair
        my ($rName,$pos,$rNext,$pNext,$nextIsReverse,$readLength) = @{$$pairs{$name}{$read}};
        $rNext eq "=" and $rNext = $rName;
        my ($low, $high);
        if($nextIsReverse){ #take the partner portion of genome as the read surrogate
            $high = $pNext + $readLength - 1;
            $low = $high - $param{readLength} + 1;
        } else {
            ($low, $high) = ($pNext, $pNext + $param{readLength} - 1);
        }
        push @$regions, [$rNext, $low, $high]; 
        push @{$$names{"$rNext:$low-$high"}}, [$name, $nextIsReverse]; 
    }
}

sub samtoolsFastaIndex { #used indexed fasta to quickly retrieve list of genome regions from genome fasta file
    my ($regions, $useChrPrefix) = @_;
    getRegionList($regions, \my$regionList, $useChrPrefix);
    $regionList or die "no target regions provided to samtoolsFastaIndex\n";
    my $refSeqFile = "$param{refSeqPath}/$param{refSeq}/$param{refSeq}.fa";
    my $indexFile = "$refSeqFile.fai";
    my $samCommand = "$param{samPath}samtools faidx $refSeqFile";
    -e $indexFile or system($samCommand); #generate index first if not already present
    my $stdout = qx/$samCommand $regionList/;
    return $stdout; #output is a fasta format list ">chr:start-end\nsequenc\n"
}

sub partnerRegions2Fasta { #takes faidx fasta format output as surrogate read2 values for unkown partners
    my ($partners, $names, $writeType, $readFiles, $pairs, $pairID) = @_;
    open my $fa1H, $writeType, $$readFiles{read1} or die "could not open $$readFiles{read1}\n"; #open passed target read files
    open my $fa2H, $writeType, $$readFiles{read2} or die "could not open $$readFiles{read2}\n";            
    my @partners = split("\n", $partners); #split faidx output into lines
    my $region;
    foreach my $partnerLine(@partners){
        if($partnerLine =~ m/^>chr(.*)/){ #fasta name line = chr region
            $region = $1;
        } elsif ($region) {  
            foreach my $indexRead(@{$$names{$region}}){
                my ($name, $nextIsReverse) = @$indexRead;
                my ($read1) = keys %{$$pairs{$name}};
                my $read2 = "\U$partnerLine";
                $nextIsReverse and $read2 = reverseComplement($read2);
                $$pairID++;
                print $fa1H ">$$pairID\n$read1\n";
                print $fa2H ">$$pairID\n$read2\n";       
            }
        }
    } 
    close $fa1H;
    close $fa2H;     
}

1;
