#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames %bowtieFields));

my $dir = "/home/wilsonte/wga/cview";

my $binSize = 0.25;
my $nBins = int(20 / $binSize);

#-------------------------------------------------------------
#my @samples = ("diploid1","diploid2","diploid3","diploid4","diploid5","diploid6",
#               "sperm01","sperm02","sperm03","sperm04","sperm05","sperm06","sperm07","sperm08","sperm09");
my @samples = ("diploid1","sperm09");
my $arraySuffix = "_NIM_HG18";
#-------------------------------------------------------------

sub covar{

    foreach my $sample(@samples){
	my %bins;
        my $table = "ARRAY_$sample$arraySuffix";  
	my $lrrBin = "round(log(2,RATIO)/$binSize)*$binSize";
	my $binSql = "SELECT bin, count(*) / (SELECT count(*) FROM $table WHERE RATIO > 0) freq
			FROM (SELECT $lrrBin bin FROM $table WHERE RATIO > 0)
			GROUP BY bin
			ORDER BY bin";
	my $binTable =  "$table\_bins";
	runSQL("CREATE TABLE $binTable AS $binSql");
	runSQL("SELECT bin, freq FROM $binTable", \my($bin,$freq));
	while (fetchRow()){ $bins{$bin}{0} = $freq;  }
	foreach my $int(1..$nBins){
		my $sql = "SELECT bin, lag(freq, $int, 0) OVER (ORDER BY bin) + lead(freq, $int, 0) OVER (ORDER BY bin) intFreq
		   	FROM $binTable
		   	ORDER BY bin";
		runSQL($sql, \my($bin,$intFreq));
		while (fetchRow()){ $bins{$bin}{$int} = $intFreq }
	}

	my $lrrSql = "SELECT CHROMOSOME, POSITION, $lrrBin bin FROM $table WHERE RATIO > 0";
	my $sql = "SELECT bin, abs(bin - lead(bin, 1, 0) OVER (ORDER BY CHROMOSOME, POSITION)) / $binSize int
  			FROM ($lrrSql)
  			ORDER BY CHROMOSOME, POSITION";
	runSQL($sql,\my($bn,$int));
	my $sum;
	while(fetchRow()){
		my $expected = $bins{$bn}{$int};
		$sum += $expected;
	}

  
	print "$sample = $sum\n";



#foreach my $bin(sort {$a <=> $b} keys %bins){
#	foreach my $int (sort {$a <=> $b} keys %{$bins{$bin}}){
#		#print join("\t", $bin, $int, $bins{$bin}{$int})."\n";
#	}
#}

	dropTable($binTable);
    }

}

1;



