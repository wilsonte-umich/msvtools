#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
defined $command{backupAllArrays} or $command{backupAllArrays} = ['singleThread', '24:00:00', 500, 0];

#backupAllArrays makes archive copies of all ARRAY and NMA files
#but specifically blocks dropping of table, i.e. is a true backup function

sub backupAllArrays {
    $param{archiveType} = 'array';
    $param{dropInput} = 0;
    my $sql = "SELECT table_name FROM user_tables WHERE table_name like 'ARRAY%'";
    runSQL($sql);
    my $arrayTables = fetchAllHashRef();  #must collect all since more SQL to be called
    foreach my $arrayTable(@$arrayTables){
        my ($tableName) = $$arrayTable{TABLE_NAME};
        if ($tableName =~ m/^ARRAY_(.*)_NIM_/ 
         or $tableName =~ m/^ARRAY_(.*)_ILL_/ 
         or $tableName =~ m/^ARRAY_(.*)_AFFY_/) {
            my $sample = $1;
            #print "$tableName\t$sample\n";
            archive $sample;
         }
    }
}

1;
