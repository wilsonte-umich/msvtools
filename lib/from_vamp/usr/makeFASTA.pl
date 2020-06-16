#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));

sub makeFASTA{


print "1\n";


    open my $fileH, ">", "$param{refSeqPath}/Ty.fa";
    my $sql = "SELECT NAME_, CHROMOSOME, START_, END_ 
               FROM SGDOTHER_SACCER2
               WHERE NAME_ LIKE '%Ty%'";
    runSQL($sql, \my($name, $chrom, $start, $end));
    while(fetchRow()){
        print $fileH ">$name\n";
        my $seq = getDNASegmentHR2($chrom, $start, $end);
        while ($seq and $seq =~ m/^(.{50})(.*)/){
            print $fileH "$1\n";
            $seq = $2;
        }
        $seq and print $fileH "$seq\n";
    }
    close $fileH;
}

1;



