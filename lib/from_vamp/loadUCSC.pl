#!/usr/bin/perl
use strict;
use warnings;

###########################################################
#loadUCSCData create Oracle tables readable by VAMP
#corresponding to, and downloaded from, UCSC data tables
#which are passes as the only input parameter
###########################################################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#The UCSC refGene tables have a systematic error in them!!
#If you retrieve exon DNA based on the exonStarts and exonEnds provided
#you always get one too many bases at the 5', but NOT the 3' end
#This problem is solved by adding 1 to the exonStart, using exonEnd as is.
#The UCSC values are imported into VAMP to prevent data inconsistency,
#but then all calling code (including the refGeneExon code below)
#must remember to act accordingly on gene and exon start values.
#In some instance, have also created a field CORRSTART_ = START_ + 1
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

our (%param, %types, %fields, %refSeqs, %reverseRefSeqs);

sub loadUCSCData {
    my ($ucscTable) = @_;
    my @uscsTables;
    if ($ucscTable eq 'all' or !$ucscTable){
        my @uscsFiles = <$param{vampPath}/UCSCSchema/$param{refSeqBase}/*.csv>;
        @uscsFiles or die "found no csv files in $param{vampPath}/UCSCSchema/$param{refSeqBase}";
        foreach my $ucscFile(@uscsFiles){
            $ucscFile =~ /$param{vampPath}\/UCSCSchema\/$param{refSeqBase}\/(.*).csv/;
            push @uscsTables, $1;
        }
    } else {
        @uscsTables  = ($ucscTable);
    }
    foreach my $ucscTable(@uscsTables){
        if ($ucscTable eq 'refGeneExons') { loadRefGeneExons() and return }
        status("\nimporting UCSC table $ucscTable\n");
        createUCSCTable($ucscTable, \my$outTable, \my%outCols, \my$outFields);
        getUCSCDataFile($ucscTable, \my$ucscFile);
        writeUCSCData(\$ucscFile, \$outTable, \%outCols, \my$outFile);
        loadData($outFile, $outTable, ',', $outFields);
        unlink($ucscFile);    
    }
}

sub createUCSCTable{
    my ($ucscTable, $outTableRef, $outColsRef, $outFieldsRef) = @_;
    status("  reading table schema...\n");
    my ($schemaFileH, @inCols, %inCols, @outCols, @outFields, $outCols, $delimiter);    
    open $schemaFileH, "<", "$param{vampPath}/UCSCSchema/$param{refSeqBase}/$ucscTable.csv" or 
    open $schemaFileH, "<", "$param{vampPath}/UCSCSchema/$param{refSeqBase}/$ucscTable.xls" or
    open $schemaFileH, "<", "$param{vampPath}/UCSCSchema/$param{refSeqBase}/$ucscTable.txt" or
        die "could not open schema file $param{vampPath}/UCSCSchema/$param{refSeqBase}/$ucscTable.csv, xls or txt";
    getUCSCLine($schemaFileH, \@inCols, \$delimiter) or die "schema file is empty!";    
    foreach my $i(0..((scalar @inCols) - 1)){ $inCols{$inCols[$i]} = $i}
    defined $inCols{vampField} or die "must specify column vampField in schema table containing Oracle-compatabile names";
    defined $inCols{vampType} or die "must specify column vampType in schema table containing Oracle-compatabile data types";
    my $i = -1;
    push @outCols, "VAMPID NUMBER";
    push @outFields, "VAMPID";
    while (getUCSCLine($schemaFileH, \@inCols, \$delimiter)){
        $i++;        
        $inCols[$inCols{vampField}] or next;
        $$outColsRef{$i} = $inCols[$inCols{vampField}];        
        push @outCols, "$inCols[$inCols{vampField}] $inCols[$inCols{vampType}]";
        push @outFields, $inCols[$inCols{vampField}];
    }
    close $schemaFileH;    
    $outCols = join(", ", @outCols);
    $$outFieldsRef = join(", ", @outFields);
    $$outTableRef = "\U$ucscTable\_$param{refSeqBase}";
    dropTable($$outTableRef);
    status("  creating VAMP table...\n");
    runSQL("CREATE TABLE $$outTableRef ($outCols)");   
}

sub getUCSCLine{
    my ($fileH, $colsRef, $delimiterRef) = @_;
    my $line = <$fileH>;
    $line or return 0;
    chomp $line;
    $line =~ s/\r//g;
    unless($$delimiterRef){ if($line =~ m/\t/) {$$delimiterRef = "\t"} else {$$delimiterRef = "," } }
    $$delimiterRef eq "\t" and $line =~ s/,/:/g;  #some UCSC fields are comma delimited, which screws up sdlldr
    @$colsRef = split($$delimiterRef, $line);
    return 1;
}

sub getUCSCDataFile{
    my ($ucscTable, $ucscFileRef) = @_;
    $$ucscFileRef = "$ucscTable.txt";
    my $gzFile = "$$ucscFileRef.gz";
    if(!wGetUCSCFile($gzFile)){
        gunzipUCSCFile($gzFile);    
    } else {
        my @chr_ucscFiles;
        foreach my $chrom(1..$refSeqs{$param{refSeqBase}}{nChrom}){
            my $chr = $reverseRefSeqs{$param{refSeqBase}}{$chrom};
            my $chr_ucscFile = "$chr\_$$ucscFileRef";        
            my $chr_gzFile = "$chr_ucscFile.gz";
            !wGetUCSCFile($chr_gzFile) or die "could not retrieve UCSC data file(s)";
            gunzipUCSCFile($chr_gzFile);
            push @chr_ucscFiles, $chr_ucscFile;
        }
        system("cat @chr_ucscFiles > $$ucscFileRef");
        unlink(@chr_ucscFiles);
    }
}

sub wGetUCSCFile{
    my ($gzFile) = @_;
    status("  downloading data file $gzFile...\n");
    return system("wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/$param{refSeqBase}/database/$gzFile' -O $gzFile");
}

sub gunzipUCSCFile{
    my ($gzFile) = @_;
    status("  unzipping data file $gzFile...\n");    
    return system("gunzip -f $gzFile");
}

sub writeUCSCData {
    my ($ucscFileRef, $outTableRef, $outColsRef, $outFileRef, $outSub) = @_;
    status("  loading data into VAMP table...\n");
    my (@inCols, $delimiter);    
    open my $inFileH, "<", $$ucscFileRef or die "could not open $$ucscFileRef";
    $$outFileRef = "$$outTableRef.csv";
    open my $outFileH, ">","$$outTableRef.csv"; 
    my %strands = ('+' => 1, '-' => 2);
    my $vampID = 0;
    LINE: while(getUCSCLine($inFileH, \@inCols, \$delimiter)){  
        $vampID++;  
        my @outCols;
        push @outCols, $vampID; 
        foreach my $i(0..((scalar @inCols) - 1)){ 
            if($$outColsRef{$i}){
                ($$outColsRef{$i} =~ m/CHROMOSOME/ or $$outColsRef{$i} =~ m/chromosome/) and $inCols[$i] = $refSeqs{$param{refSeqBase}}{$inCols[$i]};
                ($$outColsRef{$i} =~ m/STRAND/ or $$outColsRef{$i} =~ m/strand/) and $inCols[$i] = $strands{$inCols[$i]};
                defined $inCols[$i] or next LINE;                
                push @outCols, $inCols[$i]; 
            }  
        }  
        if ($outSub) {
            &$outSub($outFileH, \@outCols);
        } else {
            my $outLine = join(",", @outCols);  
            print $outFileH "$outLine\n";
        }
    } 
    close $inFileH;
    close $outFileH;
}

sub loadRefGeneExons {
    createUCSCTable('refGeneExons', \my$outTable, \my%outCols, \my$outFields);
    getUCSCDataFile('refGene', \my$ucscFile);
    writeUCSCData(\$ucscFile, \$outTable, \%outCols, \my$outFile, \&writeRefGeneExons);
    loadData($outFile, $outTable, ',', $outFields);
    unlink($ucscFile);  
    calculateRefGeneOffsets($outTable);   
    fixRefGeneStarts();
}

sub writeRefGeneExons {
    my ($outFileH, $outCols) = @_;
    my ($vampID,$name1, $chrom, $strand, $exonStarts, $exonEnds, $name2, $exonFrames) = @$outCols;
    my @exonStarts = split(":", $exonStarts); #remember, vamp replaced UCSC ',' with ':'
    my @exonEnds = split(":", $exonEnds);
    my @exonFrames = split(":", $exonFrames);
    my $nExons = scalar(@exonStarts);
    foreach my $i (0..$nExons-1) {
        $vampID += 0.0001;
        length($exonStarts) > 0 or next;
        #all exon starts and ends propagate exactly as specific by UCSC (i.e. these values carry the start error)
        print $outFileH "$vampID,$name1,$chrom,$strand,$exonStarts[$i],$exonEnds[$i],$name2,$exonFrames[$i]\n";
    }
}

sub calculateRefGeneOffsets {
    my ($exonTable) = @_;
    my $file = "$exonTable.csv";
    open my $fileH, ">", $file or die "could not open $file\n";
    writeRefGeneOffsets($exonTable, 1, "start_", $fileH);
    writeRefGeneOffsets($exonTable, 2, "end_ desc", $fileH);
    close $fileH;
    runSQL("ALTER TABLE $exonTable ADD (CORRSTART_ NUMBER, mRnaStart NUMBER)");
    runSQL("DELETE FROM $exonTable WHERE 1=1");
    loadData($file, $exonTable, ',', "VAMPID, NAME1, CHROMOSOME, STRAND, START_, CORRSTART_, END_, NAME2, FRAME, mRnaStart");
}

sub writeRefGeneOffsets {
    my ($exonTable, $exonStrand, $orderBy, $fileH) = @_; 
    #uses START_ + 1 to calculate the relative exon offsets to
    #correct the sysematic error in start value in UCSC refGene tables
    #place both the UCSC and the corrected values into the table
    runSQL("SELECT VAMPID, NAME1, CHROMOSOME, STRAND, START_, START_ + 1 corrSTART_, END_, NAME2, FRAME
            FROM $exonTable
            WHERE strand = $exonStrand
            ORDER BY name1, $orderBy",
            \my($vampID,$name1,$chrom,$strand,$start,$corrStart,$end,$name2,$frame));    
    my ($prevName1, $exonSum) = ("noName1"); 
    while(fetchRow()){
        $name1 eq $prevName1 or $exonSum = 0; 
        print $fileH "$vampID,$name1,$chrom,$strand,$start,$corrStart,$end,$name2,$frame,$exonSum\n";
        my $exonLength = $end - $corrStart + 1;
        $exonSum += $exonLength; 
        $prevName1 = $name1;
    }    
}

sub fixRefGeneStarts {
    #also update the refGene table itself to also have the corrected start available
    my $refGeneTable = "refGene_$param{refSeq}";
    runSQL("ALTER TABLE $refGeneTable ADD (CORRSTART_ NUMBER)");
    runSQL("UPDATE $refGeneTable SET CORRSTART_ = START_ + 1");  
}


1;

