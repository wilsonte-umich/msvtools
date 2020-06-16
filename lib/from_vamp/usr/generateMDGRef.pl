#!/usr/local/bin/perl
use strict;
use warnings;
use Math::Trig;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

my $theta50 = atan2(1, 1); #theta value for 50% B allele frequency, i.e. A = B
my %calls;
my %values;

my $dir = "/home/wilsonte/mdg/cel/test/celsToProcess";

sub generateMDGRef{
    my ($refSample) = @_;
    getAffyCalls($refSample, \%calls); 
    loadMDGAB();
    generateMDGRefs(); 
}

sub loadMDGAB{
    print "loading AB values...\n";
    foreach my $abFile(<$dir/*.AB.csv>){
        $abFile =~ m/$dir\/(.*)\.AB\.csv/;
        my $sample = $1;
        $sample =~ m/^F1cross(\d*)/ or next;
        open my $abFileH, "<", $abFile or die "could not open $abFile";
        my $discard = <$abFileH>;
        while (my $abLine = <$abFileH>) {
            chomp $abLine;
            $abLine =~ s/\"//g;
            my ($id, $a, $b) = split(",", $abLine); 
            ($a, $b) = (2 ** $a, 2 ** $b);  #convert to linear coordinates
            my $offset = 150;
            my ($ao, $bo) = ($a - $offset, $b - $offset); #reset the origin
            $ao < 0 and $ao = 0;
            $bo < 0 and $bo = 0;
            my ($R, $theta) = transform($ao, $bo);
            $values{$id}{$sample}{R} = $R;
            $values{$id}{$sample}{theta} = $theta;
        }
    }
}

sub generateMDGRefs{
    print "generating ref values...\n";
    open my $outH, ">", "$dir/R_theta_refs.txt";   
    print $outH "probeset_id\tR_ref\ttheta_ref\n";  
    foreach my $id (keys %values){
        my ($rSum, $rN, $thetaSum, $thetaN);
        foreach my $sample(keys %{$values{$id}}){
            $rSum += $values{$id}{$sample}{R}; #use all values to normalize R
            $rN++;        
            if ( $calls{$id}{$sample} == 1) {
                $thetaSum += $values{$id}{$sample}{theta};  #only use informative samples to normalize BAF
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

1;

