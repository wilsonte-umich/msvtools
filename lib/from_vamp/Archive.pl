#!/usr/bin/perl -w
use strict;
use warnings;

###################################################################
#Archive.pl contains scripts for dumping VAMP Oracle table to 
#flat text files, and restoring them.
#Oracle exp is NOT used so that archived files could be read
#into any other database, etc.
#File also includes script to (un)compress next-gen sequence files,
#including both raw inputs and map files.
###################################################################

#note: refSeq NOT archive since all files can just be reloaded using loadUCSC in critical failure

use vars(qw(%param %types %fields %fieldNames %refSeqs %archiveTables));

sub archive {
    my ($sample) = @_;
    -d $param{inputPath} or die "could not find directory $param{inputPath}";
    my $archiveDir = "$param{inputPath}/archive";
    -d $archiveDir or mkdir $archiveDir;
    $archiveDir .= "/$param{archiveType}";
    -d $archiveDir or mkdir $archiveDir;
    $sample or die "no sample provided\n";
    $sample = "\U$sample";
    status("archiving sample $sample...\n");
    $archiveDir .= "/$sample";
    -d $archiveDir or mkdir $archiveDir;
    my $sql;
    if ($param{archiveType} eq 'sequence'){
        $sql = "SELECT TABLE_NAME FROM USER_TABLES WHERE ";
        my $delimiter = ' ';
        #collect all tables, histogram tables, and all compound sample tables
        #(e.g. Discs) where sample is in 2nd position, i.e. the test sample
        #NOTE: '_' is a special character is Oracle!  The '#' ... escape '#'
        #escapes out the '_' to force it as a literal
        foreach my $tableType(@{$archiveTables{sequence}}){
            $sql .= " \U $delimiter TABLE_NAME LIKE '$tableType#_$sample' escape '#' ";
            $delimiter = ' OR ';
            $sql .= " \U $delimiter TABLE_NAME LIKE 'h#_$tableType#_$sample'  escape '#' ";
            $sql .= " \U $delimiter (TABLE_NAME LIKE '$tableType%' AND TABLE_NAME LIKE '%$sample') ";   
        }
    } elsif ($param{archiveType} eq 'array') {
        #collect all Array and NMA tables for single array
        #or where sample is in 2nd position, i.e. the test sample
        $sql = "\USELECT TABLE_NAME FROM USER_TABLES 
                  WHERE TABLE_NAME LIKE 'Array#_$sample#_%' escape '#'
                     OR TABLE_NAME LIKE 'NMA#_$sample' escape '#'
                     OR (TABLE_NAME LIKE 'NMA%' AND TABLE_NAME LIKE '%$sample') ";                  
    }
    $sql or die "unrecognized archiveType\n"; 
    runSQL($sql);
    my $tablesRef = fetchAllHashRef();
    scalar @$tablesRef or die "$sample has no tables in database\n";
    foreach my $tableRef(@$tablesRef){     
        my $table = $$tableRef{TABLE_NAME};
        status("  table $table");
        my $dataFile = "$archiveDir/$table.csv";
        my $headerFile = "$dataFile.header";
        open my $dataFileH, ">", $dataFile or die "could not open $dataFile";
        open my $headerFileH, ">", $headerFile or die "could not open $headerFile";
        runSQL("SELECT * FROM $table");
        my $rowCounter = 0;
        my @header;
        while (my $rowRef = fetchRowHashRef()){
            unless ($rowCounter) {
                @header = sort {$a cmp $b} keys %$rowRef;
                print $headerFileH join(",", @header)."\n";
            }
            my $delimiter = "";
            foreach my $column(@header) { 
                my $val = $$rowRef{$column};
                defined $val or $val = "";
                print $dataFileH "$delimiter$val";
                $delimiter = ",";
            }
            print $dataFileH "\n";
            $rowCounter++;
        }
        close $dataFileH;
        close $headerFileH;
        my $zipFile = "$dataFile.gz";
        -e $zipFile and unlink $zipFile;
        system ("gzip -q $dataFile");
        if (-e $zipFile) {
            $param{dropInput} and dropTable($table);
            status(", archived $rowCounter rows\n");
        }  
        if ($rowCounter == 0){
            status("    deleting empty table archive...\n");
            unlink $zipFile;
            unlink $headerFile;
        }
    }
}

my %tableTypes;
sub restore {
    my ($sample) = @_;
    $sample or die "no sample provided\n";
    $sample = "\U$sample";
    status("restoring sample $sample...\n"); 
    my $archiveDir = "$param{inputPath}/archive/$param{archiveType}/$sample";
    -d $archiveDir or die "could not find directory $archiveDir";    
    foreach my $tableType(keys %fields){$tableTypes{"\U$tableType"} = $tableType;}
    foreach my $zipFile (<$archiveDir/*.csv.gz>){  
        system ("gunzip -q $zipFile");
        $zipFile =~ m/($archiveDir\/.*\.csv)\.gz/;
        my $dataFile = $1;
        my $headerFile = "$dataFile.header";
        open my $fileH, "<", $headerFile or die "could not open $headerFile";
        my $header = <$fileH>;
        close $fileH;
        $dataFile =~ m/$archiveDir\/(.*).csv/;
        my $table = $1;
        $table =~ m/^(.*?)_(.*)/;
        my $tableType = $1;
        my $tableSample = $2;
        $tableTypes{$tableType} or next;
        my $newTable = newTable($tableTypes{$tableType}, "$tableSample");
        loadData2($dataFile, $newTable, ",", $header, 1);   
        if ($param{dropInput}){
            unlink($dataFile);
        } else {
            system ("gzip -q $dataFile");
        }
    }
}

sub compress {
    my ($sample) = @_;
    $sample or die "no sample provided\n";
    status("compressing sequence and map files for sample $sample...\n"); 
    my $sampleDir = "$param{inputPath}/$sample";
    -d $sampleDir or die "could not find directory $sampleDir";   
    compressSequenceFiles($sampleDir);
}

my @compressTypes = ('sequence.txt', 'qseq.txt', '.fa', '.fa.p', 
                     '.gff', '.bowtie', '.bowtie.not_aligned');
sub compressSequenceFiles{
    my ($dir) = @_;
    FILE: foreach my $file(<$dir/*>){
        if (-d $file){
            compressSequenceFiles($file);
        } else {
            foreach my $extension(@compressTypes){
                if ($file =~ m/.*$extension$/){
                    status("  $file\n"); 
                    my $zipFile = "$file.gz";
                    unless(-e $zipFile){ system ("gzip -q $file"); } 
                    next FILE;
                }
            }
        }
    } 
}

sub uncompress {
    my ($sample) = @_;
    $sample or die "no sample provided\n";
    status("uncompressing sequence and map files for sample $sample...\n"); 
    my $sampleDir = "$param{inputPath}/$sample";
    -d $sampleDir or die "could not find directory $sampleDir";   
    uncompressSequenceFiles($sampleDir);
}

sub uncompressSequenceFiles{
    my ($dir) = @_;
    FILE: foreach my $file(<$dir/*>){
        if (-d $file){
            uncompressSequenceFiles($file);
        } else {
            foreach my $extension(@compressTypes){
                if ($file =~ m/.*$extension.gz$/){
                    status("  $file\n"); 
                    system ("gunzip -q $file");
                    next FILE;
                }
            }
        }
    } 
}

1;

