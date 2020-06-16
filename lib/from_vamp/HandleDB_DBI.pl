#!/usr/bin/perl
use strict;
use warnings;
use DBI;

#########################################################################
#HandleDB_DBI is a collection of functions for db input and output
#########################################################################

use vars(qw(%param %fields %types %partitions));
my $dbh;
my $sth;

#########################################################################
#opening and closing database connection
#------------------------------------------------------------------------
sub openOracle{
    my ($user, $leftover) = split("/", $param{dbLogin});
    my ($pw, $db) = split("@", $leftover);
    $dbh = DBI->connect("dbi:Oracle:", $param{dbLogin}, '', {AutoCommit => 1, PrintError=>1})
    or $dbh = DBI->connect("dbi:Oracle:host=localhost;port=1521;service_name=$db", $user, $pw, {AutoCommit => 1, PrintError=>1})   
    or die "openOraclefailed";
}
#------------------------------------------------------------------------
sub closeDB{
    if ($sth) {$sth->finish};
    $dbh->disconnect or die "closeDB error: $DBI::errstr";
}
#########################################################################


#########################################################################
#executing SQL and returning SELECT results
#------------------------------------------------------------------------
#establish a statement handle based on input SQL and execute it
sub runSQL{  
    my ($sql, @bindColRefs) = @_; #@bindColRefs = optional array containing references to scalars to receive the retrieved data
    $sth = $dbh->prepare($sql);
    $sth->execute;
    @bindColRefs and $sth->bind_columns(@bindColRefs);
    return $sth;  #VAMP currently never captures the returned $sth, but you may want to
}
#------------------------------------------------------------------------
#functions to return data; each assume you have already used runSQL to create a statement handle
sub fetchRow{ #return the next row with data bound to the scalars specified with @bindColRefs in runSQL
    return $sth->fetch;
}
#------------------------------------------------------------------------
sub fetchRowArray{ #return the next row as an array list
    return $sth->fetchrow_array;
}
#------------------------------------------------------------------------
sub fetchRowHashRef{ #return the next row as a hash reference
    return $sth->fetchrow_hashref;
}
#------------------------------------------------------------------------
sub fetchAllHashRef{ #returns a reference to an array that contains one hash reference per row
    return $sth->fetchall_arrayref({})
}
#########################################################################


#########################################################################
#table management
#------------------------------------------------------------------------
sub newTable{
    my ($tableType, $sample) = @_;
    #get table and type information from Schema.pl
    my $table = getTableName($tableType, $sample);
    my $fields = $fields{$tableType};
    my $partitionsRef = $partitions{$tableType};
    my ($chromPartField, $typeSubPartField, $typeSubPartsRef, $partByPosField) = @$partitionsRef;
    #drop table if it already exists
    tableExists($table) and dropTable($table); 
    #parse the partitioning information, see Schema.pl for details
    my $partitionList = "";
    if ($chromPartField){
        my $subPartitionList = "";
        if ($typeSubPartField and $typeSubPartsRef) {
            my %typeSubParts = %$typeSubPartsRef;
            $subPartitionList = " SUBPARTITION BY LIST ($typeSubPartField)
                                    SUBPARTITION TEMPLATE(" ;
            foreach my $subPartition (keys %typeSubParts){
                #apparently 'Normal' became a reserver word in Oracle
                #fix by adding '_' to all subpartition names
                $subPartitionList .= " SUBPARTITION $subPartition\_ VALUES ($typeSubParts{$subPartition}),"    
            }
            $subPartitionList = substr($subPartitionList, 0, length($subPartitionList) - 1);
            $subPartitionList .= ")";
        } elsif ($partByPosField) {
            $subPartitionList = " SUBPARTITION BY LIST ($partByPosField)
                                  SUBPARTITION TEMPLATE(
                                    SUBPARTITION bin1 VALUES (0),
                                    SUBPARTITION bin2 VALUES (1),
                                    SUBPARTITION bin3 VALUES (2),
                                    SUBPARTITION bin4 VALUES (3),
                                    SUBPARTITION bin5 VALUES (4),
                                    SUBPARTITION bin6 VALUES (5),
                                    SUBPARTITION bin7 VALUES (6),
                                    SUBPARTITION bin8 VALUES (7),
                                    SUBPARTITION bin9 VALUES (8),
                                    SUBPARTITION bin10 VALUES (9),
                                    SUBPARTITION bin11 VALUES (10),
                                    SUBPARTITION bin12 VALUES (11),
                                    SUBPARTITION bin13 VALUES (12))";
        }
        $partitionList = " PARTITION BY RANGE ($chromPartField)
                            $subPartitionList
                            (PARTITION Chr1_2 VALUES LESS THAN (3),
                             PARTITION Chr3_4 VALUES LESS THAN (5),
                             PARTITION Chr5_6 VALUES LESS THAN (7),
                             PARTITION Chr7_8 VALUES LESS THAN (9),
                             PARTITION Chr9_10 VALUES LESS THAN (11),
                             PARTITION Chr11_12 VALUES LESS THAN (13),
                             PARTITION Chr13_14 VALUES LESS THAN (15),
                             PARTITION Chr15_16 VALUES LESS THAN (17),
                             PARTITION Chr17_18 VALUES LESS THAN (19),
                             PARTITION Chr19_20 VALUES LESS THAN (21),
                             PARTITION Chr21_22 VALUES LESS THAN (23),
                             PARTITION Chr23_24 VALUES LESS THAN (MAXVALUE))"; 
    }
    #create and return the table
    runSQL("CREATE TABLE $table ($fields) $partitionList");
    return $table;
}
#------------------------------------------------------------------------
sub dropTable{
    my(@tables) = @_;
    scalar(@tables) or return;
    foreach my $table(@tables){
        tableExists($table) or next;
        runSQL("DROP TABLE $table PURGE"); 
    }
}
sub dropView{
    my($view) = @_;
    $view or return;
    viewExists($view) or return;
    runSQL("DROP VIEW $view");   
}

#------------------------------------------------------------------------
sub updateTable{  #for bulk updates of all fields in a table, MUCH faster than looping "INSERT INTO"
    my ($tableType, $table, $selectSql) = @_;
    my $tableTMP = newTable($tableType, $table); #TMP table called type_type_sample
    runSQL("INSERT INTO $tableTMP $selectSql");    
    dropTable($table);
    runSQL("RENAME $tableTMP TO $table");
}
#------------------------------------------------------------------------
sub tableExists{
    my($table) = @_;
    runSQL("SELECT 1 TABLEEXISTS FROM USER_TABLES WHERE TABLE_NAME = '\U$table'");
    return fetchRowArray();
}
sub viewExists{
    my($view) = @_;
    runSQL("SELECT 1 VIEWEXISTS FROM USER_VIEWS WHERE VIEW_NAME = '\U$view'");
    return fetchRowArray();
}
#########################################################################


#########################################################################
#bulk loading of data from files into Oracle using sqlldr
#writing to a file and bulk loading a table is MUCH faster than looping "INSERT INTO"
#------------------------------------------------------------------------
sub loadData{
    my ($dataFile, $table, $delimiter, $fieldsList, $suppressDrop) = @_;
    my $controlFile = "$dataFile.ctl";
    my $logFile = "$dataFile.log";
    my $badFile = "$dataFile.bad";
    open my $controlFileH, ">", $controlFile;
    print $controlFileH "unrecoverable load data\n";
    print $controlFileH "infile $dataFile\n";
    print $controlFileH "append\n";
    print $controlFileH "into table $table\n";
    print $controlFileH "fields terminated by \"$delimiter\"\n";
    print $controlFileH "($fieldsList)\n";    
    close $controlFileH;
    system("sqlldr $param{dbLogin} control=$controlFile log=$logFile bad=$badFile silent=all direct=true errors=10000");
    unlink $controlFile;
    $suppressDrop or unlink $dataFile;
}
sub loadData2{
    my ($dataFile, $table, $delimiter, $fieldsList, $suppressDrop) = @_;
    my $controlFile = "$dataFile.ctl";
    my $logFile = "$dataFile.log";
    my $badFile = "$dataFile.bad";
    open my $controlFileH, ">", $controlFile;
    print $controlFileH "unrecoverable load data\n";
    print $controlFileH "infile '$dataFile'\n";
    print $controlFileH "append\n";
    print $controlFileH "into table $table\n";
    print $controlFileH "fields terminated by \"$delimiter\"\n";
    print $controlFileH "($fieldsList)\n";    
    close $controlFileH;
    system("sqlldr $param{dbLogin} control=$controlFile log=$logFile bad=$badFile silent=all direct=true errors=10000");
    unlink $controlFile;
    $suppressDrop or unlink $dataFile;
}
#########################################################################


1;














