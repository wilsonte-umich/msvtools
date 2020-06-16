#!/usr/local/bin/perl
use strict;
use warnings;
use Math::Trig;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

my $theta50 = atan2(1, 1); #theta value for 50% B allele frequency, i.e. A = B
my %calls;

#$param{inputPath} = "/home/shared/data";

sub generateAffyRef{
    my ($refSample) = @_;
    getAffyCalls($refSample, \%calls);  
    generateRefs($refSample); 
}

sub getAffyCalls{
    my ($refSample, $calls) = @_;
    print "getting affy calls...\n";
    my ($fileH, $header) = getAffyHeader("$param{inputPath}/$refSample/quant-norm.pm-only.brlmm-p.calls.txt", "\t");
    while (my $line = <$fileH>) {
        chomp $line;
        my @line = split("\t", $line);  
        my $id = $line[$$header{probeset_id}];
        foreach my $sample (keys %$header){
            $sample eq "probeset_id" and next;
            $$calls{$id}{$sample} = $line[$$header{$sample}]; 
        }
    }    
    close $fileH;
}

sub generateRefs{
    my ($refSample) = @_;
    print "generating ref values...\n";
    my ($fileH, $header) = getAffyHeader("$param{inputPath}/$refSample/quant-norm.pm-only.brlmm-p.summary.txt", "\t");
    open my $outH, ">", "$param{inputPath}/$refSample/R_theta_refs.txt";   
    print $outH "probeset_id\tR_ref\ttheta_ref\n";    
    while (my $line = <$fileH>) {
        my ($id, $alleleA, $aVals) = parseSummaryLine($line, $header);
        $alleleA eq "A" or  die "allele A not where expected\n";
        $line = <$fileH>;
        my ($idB, $alleleB, $bVals) = parseSummaryLine($line, $header);
        $id eq $idB or die "failed id pairing in summary file\n";
        my ($rSum, $rN, $thetaSum, $thetaN);
        foreach my $sample (keys %{$header}){
        
#            $sample eq "probeset_id" and next;
            $sample =~ m/^F1cross(\d*)/ or next;
            $1 eq "12" and next; ############## only use F1cross, NOT #12 which is a test sample ######################
            
            my $a = $$aVals[$$header{$sample}];
            my $b = $$bVals[$$header{$sample}];
            my ($R, $theta) = transform($a, $b);
            $rSum += $R; #use all values to normalize R
            $rN++;
            if ( $calls{$id}{$sample} == 1) {
                $thetaSum += $theta;  #only use informative samples to normalize BAF
                $thetaN++;   
            }  
        }  
        my $rRef = $rSum / $rN;
        my $thetaRef = $theta50;
        $thetaN and $thetaRef = $thetaSum / $thetaN;
        print $outH "$id\t$rRef\t$thetaRef\n";    
    }
    close $outH;  
}

sub parseSummaryLine {
    my ($line, $header) = @_;
    chomp $line;
    my @line = split("\t", $line); 
    my $id = $line[$$header{probeset_id}];
    $id =~ m/(.*)-(.)/ or die "bad summary id on line:\n$line\n";
    return ($1, $2, \@line);
}

sub transform {
    my ($a, $b) = @_;
    my ($R, $theta);
    $R = $a + $b;
    $theta = atan2($b, $a);
    return ($R, $theta);
}

sub getAffyHeader{
    my ($file, $delimiter) = @_;
    open my $fileH, "<", $file or die "could not open $file";
    my %header;
    while (my $line = <$fileH>){
        chomp $line;
        $line =~ s/\"//g;
        if ( !($line =~ m/^\#/) ) {
            my @header = split($delimiter, $line);
            foreach my $i (0..((scalar @header) - 1)) { $header{$header[$i]} = $i }
            return ($fileH, \%header);   
        }
    }
}

1;

