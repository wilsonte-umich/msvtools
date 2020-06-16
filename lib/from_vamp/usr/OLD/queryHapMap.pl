#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames %bowtieFields));

my ($runDir, @samples, $fileH, %dataFields);

sub queryHapMap{

    my %match090;
    my %counts;
    my %matches;

    runSQL("SELECT TABLE_NAME FROM USER_TABLES WHERE TABLE_NAME like 'ARRAY_HM%'");
    my $tablesRef = fetchAllHashRef();
    foreach my $tableRef(@$tablesRef){
        my $table = $$tableRef{TABLE_NAME};

        runSQL("
SELECT q.CHROMOSOME || ':' || q.START_ || ':' || q.END_ id_,
       q.match090,
       case when q.cnvType * avg(log(2, t.RATIO)) >= 0.15 then 1 else 0 end cnvMatch
FROM
(SELECT CHROMOSOME, START_, END_, 
        case when NORMALIZEDRATIO < 0 then -1 else 1 end cnvType, 
        case when score_090 >= 5 then 1 else 0 end match090 FROM INT_7974_090) q, 
$table t
WHERE q.chromosome = t.chromosome
  AND t.POSITION >= q.START_ 
  AND t.POSITION <= q.END_
GROUP BY q.CHROMOSOME, q.START_, q.END_, q.match090, q.cnvType", \my($id, $match090, $cnvMatch));
        
        while (fetchRow()){
            $match090{$id} = $match090;
            $counts{$id}++;
            $matches{$id} += $cnvMatch;
        } 
    }
    
    my %matches090;
    my %unmatches090;
    foreach my $id(keys %match090) {
        my $freq = int( ($matches{$id} / $counts{$id} / 0.05) + 0.025) * 0.05;
        if ($match090{$id}) {
            $matches090{$freq}++;
        } else {
            $unmatches090{$freq}++;
        }
    }
    
    print "matches 090\n";
    foreach my $freq ( sort { $a <=> $b } keys %matches090) {
        print "$freq\t$matches090{$freq}\n";
    }

    print "\nunmatches 090\n";
    foreach my $freq ( sort { $a <=> $b } keys %unmatches090) {
        print "$freq\t$unmatches090{$freq}\n";
    }

}



1;



