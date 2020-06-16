#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %gffFields %bowtieFields));

#callable parameters and commands in this script
defined $param{textToPrint} or $param{textToPrint} = '';
defined $command{printText} or $command{printText} = ['singleThread', '10:00', 1000, 0];
defined $command{printTables} or $command{printTables} = ['singleThread', '10:00', 1000, 0];

sub printText{
    my (@toPrint) = @_;
    print join(" ", @toPrint)."\n";
    $param{textToPrint} or $param{textToPrint} = 'parameter \'textToPrint\' was not specified';
    print "$param{textToPrint}\n";
    print "vampPath = $param{vampPath}\n";
    print "passPath = $param{passPath}\n";
    print "bowtiePath = $param{bowtiePath}\n";
    print "refSeqPath = $param{refSeqPath}\n";
    print "blastPath = $param{blastPath}\n";      
    print "inputPath = $param{inputPath}\n"; 
}


sub printTables{
    runSQL("SELECT table_name FROM user_tables ORDER BY table_name", \my$tableName);
    while(fetchRow()){ status("$tableName\n") }
}


1;














