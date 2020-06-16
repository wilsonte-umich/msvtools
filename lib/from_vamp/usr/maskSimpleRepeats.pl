#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));

my @repLen = (1,2,3,4);
my %fNotAllowed = (1=>0.65,2=>0.55,3=>0.5,4=>0.5);
my %repeats = (1=>['A','C','G','T']);
foreach my $repLen(2..4){
    foreach my $base(@{$repeats{1}}){
        foreach my $ext(@{$repeats{$repLen-1}}){
            push @{$repeats{$repLen}}, "$base$ext"
        }    
    }
}

sub purgeSimpleRepeats{
    my ($inFile) = @_;
    open my $inFileH, "<", $inFile;
    open my $outFileH, ">", "$inFile.tmp";
    while(my $nameline = <$inFileH>){
        $nameline =~ m/^>/ or next;
        my $readLine = <$inFileH>;    
        unless(checkSimpleRepeats($readLine)){print $outFileH "$nameline$readLine"} 
    }
    close $inFileH;
    close $outFileH;
    system("mv $inFile $inFile.backup");
    system("mv $outFileH $inFile"); 
}

sub checkSimpleRepeats{
    my ($read) = @_;
    chomp $read;
    my $readLen = length($read);
    foreach my $repLen(@repLen){    
        my $scalar = $repLen / $readLen;
        my $fNotAllowed = $fNotAllowed{$repLen};
        foreach my $repeat(@{$repeats{$repLen}}){
            scalar(my @count = $read =~ m/$repeat/g) * $scalar >= $fNotAllowed and return $repeat
        }
    }
    return undef;
}

#sub countSimpleRepeats{
#    my ($sample) = @_;
#    getFixedReadFiles($sample, \my %readFiles);
#    getFixedMapFiles($sample, \my %mapFiles);
#    open my $outH, ">", "$param{vampPath}/simpleRepeatCounts.xls";
#    foreach my $type ('fx','ufx','random'){   
#        print $outH "\n$type\n";
#        my %maxRepeats;        
#        my $maxI = 100000;
#        if($type eq 'random'){
#            foreach my $i(1..$maxI){
#                my $read;
#                foreach my $base(1..36){$read .= $repeats{1}[int(rand(4))]}  
#                getMaxRepeats($read, \%maxRepeats);                
#            }
#        } else {
#            open my $fileH, "<", $readFiles{$type}{read1}; 
#            my $i;   
#            while(my $line = <$fileH>){
#                $line =~ m/^>/ or next;
#                my $read = <$fileH>;
#                getMaxRepeats($read, \%maxRepeats);
#                $i++;
#                $i > $maxI and last;
#            }
#            close $fileH;          
#        }
#        foreach my $repLen(sort {$a <=> $b} @repLen){
#            foreach my $f (sort {$a <=> $b} keys %{$maxRepeats{$repLen}}){
#                print $outH "$repLen\t$f\t${$maxRepeats{$repLen}}{$f}\n"
#            }
#        }     
#    }
#    close $outH;
#}

#sub getMaxRepeats{
#    my ($read, $maxRepeats) = @_;
#    chomp $read;
#    my $readLen = length($read);
#    foreach my $repLen(@repLen){    
#        my $scalar = $repLen / $readLen;
#        my $maxF = 0;
#        foreach my $repeat(@{$repeats{$repLen}}){
#            my $f = scalar(my @count = $read =~ m/$repeat/g) * $scalar;
#            $f > $maxF and $maxF = $f;
#        }
#        $$maxRepeats{$repLen}{(int($maxF/0.05))*0.05}++;
#    }
#}


1;


