#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %bowtieFields));

sub dbExport{
    my $exportPath = "$param{vampPath}/dbExport";
    my $tablesFile = "$exportPath/exportList.txt";   
    print $tablesFile, "\n";
    open my $tablesFileH, "<", $tablesFile;
    while (my $table = <$tablesFileH>) {
        $table =~ s/\n//g;
        $table =~ s/\r//g;
        print $table, "\n";
        tableExists($table) or next;
        my $exportFile = "$exportPath/$table.xls";
        open my $exportFileH, ">",  $exportFile;
        my $sth = runSQL("SELECT * FROM $table");
        while ( my @row = ($sth->fetchrow_array) ) {
            my $delimiter = "";
            foreach my $value (@row){
                (defined $value) or $value = "";
                print $exportFileH "$delimiter$value";
                $delimiter = "\t";
            }
            print $exportFileH "\n";
        }
        close $exportFileH;
    }
    close $tablesFileH;
}

1;


