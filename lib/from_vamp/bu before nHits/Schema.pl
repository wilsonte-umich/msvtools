#!/usr/bin/perl -w
use strict;
use warnings;

#########################################################################
#Schema.pl is a storehouse for script requires file naming conventions and db schema
#########################################################################

use vars(qw(%param %types %fields %fieldNames %partitions));

#########################################################################
#load vamp required scripts
#------------------------------------------------------------------------
-d "$param{vampPath}/bin" or die "\n$param{vampPath} is not a valid VAMP path\n";
foreach my $script(<$param{vampPath}/bin/*.pl>){ $script =~ m/ExecuteInstruction.pl/ or require $script }
#########################################################################


#########################################################################
#load user-supplied scripts
#------------------------------------------------------------------------
#load any scripts found in 'XXX/vamp/bin/usr' where XXX is either the
#home directory or vampPath, with the home directory taking precedence
#if a script of the same name is found in both locations
#Any subs in usr script must be accessed by command 'run subName parameterList'
#------------------------------------------------------------------------
my $path = "$ENV{HOME}/vamp/bin/usr";
-d $path or next;
foreach my $script(<$path/*.pl>){ require $script }
#########################################################################


#########################################################################
#subs that determine naming conventions
#------------------------------------------------------------------------
sub getDirectories{
    my ($sample, $dirsRef) = @_;
    $$dirsRef{sample} = "$param{inputPath}/$sample"; 
    $$dirsRef{read1} = "$$dirsRef{sample}/1";
    $$dirsRef{read2} = "$$dirsRef{sample}/2";
    $$dirsRef{qseq1} = "$$dirsRef{read1}/qseq";
    $$dirsRef{qseq2} = "$$dirsRef{read2}/qseq";
}
sub getReadFiles{
    my ($sample, $readType, $filesRef) = @_;
    my %dirs;
    getDirectories($sample, \%dirs);
    if ($readType eq 'solexa_fastq'){
        $$filesRef{read1} = "$dirs{read1}/$sample\_1_sequence.txt";
        $$filesRef{read2} = "$dirs{read2}/$sample\_2_sequence.txt"; 
    } elsif ($readType eq 'fastq'){
        $$filesRef{read1} = "$dirs{read1}/$sample\_1.fq";
        $$filesRef{read2} = "$dirs{read2}/$sample\_2.fq";        
    } elsif ($readType eq 'fasta') {
        $$filesRef{read1} = "$dirs{read1}/$sample\_1.fa";
        $$filesRef{read2} = "$dirs{read2}/$sample\_2.fa";
    } elsif ($readType eq 'purgedFasta') {
        $$filesRef{read1} = "$dirs{read1}/$sample\_1.fa.p";
        $$filesRef{read2} = "$dirs{read2}/$sample\_2.fa.p";        
    } else {
        die "'$readType' is not a recognized readType"
    }
}
sub getMapFiles{
    my ($sample, $filesRef) = @_;
    my %dirs;
    getDirectories($sample, \%dirs);
    my $suffix;
    if ($param{mapType} eq 'pass') {
        $suffix = 'gff';
    } elsif ($param{mapType} eq 'bowtie') {
        $suffix = 'bowtie';
    }
    $suffix or die "$param{mapType} is not a recognized mapType";
    $$filesRef{read1} = "$dirs{read1}/$sample\_1.$suffix";
    $$filesRef{read2} = "$dirs{read2}/$sample\_2.$suffix";    
}
sub getTableName{
    my ($tableType, $sample) = @_;
    return "\U$tableType\_$sample"; #force lane name to upper case for Oracle
}
#########################################################################


#########################################################################
#type enumerations
#------------------------------------------------------------------------
%types = (Pairs => {Convergent=>0, Divergent=>1, Colinear=>2, DiffChrom=>3, SingleRead=>4},
          Discs => {None=>0, A=>1, C=>2, G=>3, T=>4, '.' => 5, '?'=>5, N=>5, X=>5, Missing=>6, Extra=>7},
          Frags => {Discard => -1, Normal=>0, Deletion=>1, Insertion=>2, Inversion=>4, Duplication=>8,
                    DiffChrom=>16, ReverseNormal=>32, TooSmall=>64, Rooted=>128, SingleRead=>256}, 
          Consequences => {Silent=>0, Missense=>1, Nonsense=>2, InDel=>3},
          RevDiscs => {1=>'A', 2=>'C', 3=>'G', 4=>'T', 5=>'N'}, 
          RevDiscsFull => {0=>'-', 1=>'A', 2=>'C', 3=>'G', 4=>'T', 5=>'N', 6=>'Missing', 7=>'Extra'} );
$types{Sets} = $types{Frags};
$types{Sets}{Z} = 3;
#########################################################################


#########################################################################
#table field definitions
#------------------------------------------------------------------------
$fields{Reads} = "PAIRID NUMBER, MAPS VARCHAR2(4000)";
$fields{Pairs} = "PAIRTYPE NUMBER,
                    CHROMOSOME1 NUMBER, POSITION1 NUMBER, LENGTH1 NUMBER, STRAND1 NUMBER, DISCREPANCIES1 NUMBER,
                    CHROMOSOME2 NUMBER, POSITION2 NUMBER, LENGTH2 NUMBER, STRAND2 NUMBER, DISCREPANCIES2 NUMBER,
                    FRAGMENTSIZE NUMBER, PAIRID NUMBER, NFRAGS NUMBER";
$fields{Discs} = "DISCREPANCYID NUMBER, 
                    CHROMOSOME NUMBER, POSITION NUMBER, DISCREPANCYTYPE NUMBER, EXTRA NUMBER,
                    S1_SAMECOUNT NUMBER, S1_CONSISTENTCOUNT NUMBER, S1_READCOUNT NUMBER, 
                    S2_SAMECOUNT NUMBER, S2_CONSISTENTCOUNT NUMBER, S2_READCOUNT NUMBER,
                    SAMELOD NUMBER(*,3), CONSISTENTLOD NUMBER(*,3), CONSEQUENCE NUMBER"; 
$fields{Stats} = "STATNAME VARCHAR2(255), STATVALUE NUMBER(*,3)";
$fields{Frags} = "FRAGMENTID NUMBER, FRAGMENTTYPE NUMBER,
                    CHROMOSOME1 NUMBER, POSITION1 NUMBER, LENGTH1 NUMBER, STRAND1 NUMBER, DISCREPANCIES1 NUMBER,
                    CHROMOSOME2 NUMBER, POSITION2 NUMBER, LENGTH2 NUMBER, STRAND2 NUMBER, DISCREPANCIES2 NUMBER,
                    FRAGMENTSIZE NUMBER, PAIRID NUMBER, NFRAGS NUMBER,
                    EVENTSIZE NUMBER, STDEVNORMAL NUMBER, ENDTOLERANCE NUMBER, 
                    NSETSFRAG NUMBER, NSETSPAIR NUMBER";    
$fields{FMap} = "CHROMOSOME NUMBER, POSITION NUMBER, COVERAGE NUMBER, NORMALIZEDCOVERAGE NUMBER(*,3)";
$fields{RMap} = $fields{FMap};  
$fields{ZMap} = "CHROMOSOME NUMBER, POSITION NUMBER, COVERAGE NUMBER, MEAN NUMBER, Z NUMBER(*,1)";                         
#$fields{FRatio} = "CHROMOSOME NUMBER, POSITION NUMBER,
#                    COVERAGERATIO NUMBER(*,3), COVERAGERATIO1E4 NUMBER(*,3), COVERAGERATIO1E5 NUMBER(*,3)";                            
$fields{Sets} = "SETID NUMBER, SETTYPE NUMBER,
                    CHROMOSOME1 NUMBER, CHROMOSOME2 NUMBER,
                    SPANSTART NUMBER, SPANEND NUMBER, OVERLAPSTART NUMBER, OVERLAPEND NUMBER,
                    FRAGMENTMEAN NUMBER, EVENTMEAN NUMBER, STDEV NUMBER, 
                    MEANDELTA NUMBER(*,3), FRACTIONDELTA2 NUMBER(*,3), FRACTIONOVERLAPLEFT NUMBER(*,3), FRACTIONOVERLAPRIGHT NUMBER(*,3),
                    STRAND1 NUMBER, STRAND2 NUMBER,
                    NFRAGSTOTAL NUMBER, NFRAGSSAMPLE NUMBER, FRACBADFRAGS NUMBER(*,3),
                    MARK NUMBER, DESCRIPTION VARCHAR2(255)";                
$fields{IDs} = "SETID NUMBER, ID_ NUMBER";
$fields{CDS} = "CDSID NUMBER, NAME1 VARCHAR2(255), NAME2 VARCHAR2(255),
                CHROMOSOME NUMBER, STRAND NUMBER, EXON NUMBER,
                START_ NUMBER, CARRYOVER NUMBER, END_ NUMBER";
$fields{NCDS} = "NCDSID NUMBER, NAME1 VARCHAR2(255), NAME2 VARCHAR2(255),
                CHROMOSOME NUMBER, STRAND NUMBER, EXON NUMBER, 
                START_ NUMBER, END_ NUMBER";
$fields{h} = "SERIES NUMBER, X NUMBER, Y NUMBER"; #h = histogram, kept short since name appended to parent table name
$fields{rmsk} = "CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER, STRAND NUMBER, 
                 REPEATNAME VARCHAR(255), REPEATCLASS VARCHAR(255), REPEATFAMILY VARCHAR(255)";    
$fields{CNVs} = "CNVNAME VARCHAR2(255), CNVTYPE NUMBER, CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER";          
$fields{Array} = "CHROMOSOME NUMBER, POSITION NUMBER, RATIO NUMBER(*,5), NORMALIZEDRATIO NUMBER(*,5),
                  BFREQUENCY NUMBER(*,5), ZYGOSITY NUMBER(*,2), NORMALIZEDZYGOSITY NUMBER(*,2), CNVVALUE NUMBER(*,5),
                  GCSCORE NUMBER(*,3), X NUMBER(*,3), Y NUMBER(*,3), THETA NUMBER(*,3), R NUMBER(*,3), BF NUMBER(*,3)";          
$fields{NMA} = "NMAID NUMBER, NMATYPE VARCHAR(255), NPROBES NUMBER, NSTDDEV NUMBER(*,1), 
                CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER, CNVSIZE NUMBER, NPROBESINSET NUMBER, 
                TESTMEAN NUMBER(*,3), TESTSTDEV NUMBER(*,3), REFERENCEMEAN NUMBER(*,3), REFERENCESTDEV NUMBER(*,3), 
                NORMALIZEDRATIO NUMBER(*,3), NORMALIZEDZYGOSITY NUMBER(*,3), ZSCORE NUMBER(*,1),
                COPYNUMBER NUMBER, DESCRIPTION VARCHAR2(255) ";  
$fields{LD} = "CHROMOSOME NUMBER, REFERENCEPOSITION NUMBER, IMPUTEDPOSITION NUMBER, 
                DPRIME NUMBER(*,3), RSQUARED NUMBER(*,3), LOD NUMBER(*,3)"; 
$fields{AF} = "CHROMOSOME NUMBER, POSITION NUMBER, UNINFORMATIVEFREQ NUMBER(*,3)"; 
$fields{Events} = "EVENTID NUMBER, EVENTTYPE NUMBER,
                    CHROMOSOME NUMBER, SPANSTART NUMBER, SPANEND NUMBER, 
                    NSETS NUMBER ";  
$fields{CNCs} = "CNCID NUMBER, EVENTID NUMBER, CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER, 
                 COPYNUMBER NUMBER, NORMALIZEDCOVERAGE NUMBER(*,3) ";  
                      
$fieldNames{Sets} = "SETID, SETTYPE,
                    CHROMOSOME1, CHROMOSOME2,
                    SPANSTART, SPANEND, OVERLAPSTART, OVERLAPEND,
                    FRAGMENTMEAN, EVENTMEAN, STDEV, 
                    MEANDELTA, FRACTIONDELTA2, FRACTIONOVERLAPLEFT, FRACTIONOVERLAPRIGHT,
                    STRAND1, STRAND2,
                    NFRAGSTOTAL, NFRAGSSAMPLE, FRACBADFRAGS,
                    MARK, DESCRIPTION"; 
$fieldNames{CNVs} = "CNVNAME, CNVTYPE, CHROMOSOME, START_, END_"; 
$fieldNames{Array} = "CHROMOSOME, POSITION, RATIO, NORMALIZEDRATIO, BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, 
                      CNVVALUE, GCSCORE, X, Y, THETA, R, BF";         
$fieldNames{NMA} = "NMAID, NMATYPE, NPROBES, NSTDDEV, 
                    CHROMOSOME, START_, END_, CNVSIZE, NPROBESINSET, 
                    TESTMEAN, TESTSTDEV, REFERENCEMEAN, REFERENCESTDEV, 
                    NORMALIZEDRATIO, NORMALIZEDZYGOSITY, ZSCORE,
                    COPYNUMBER, DESCRIPTION";
$fieldNames{LD} = "CHROMOSOME, REFERENCEPOSITION, IMPUTEDPOSITION, DPRIME, RSQUARED, LOD"; 
$fieldNames{AF} = "CHROMOSOME, POSITION, UNINFORMATIVEFREQ";    
$fieldNames{Events} = "EVENTID, EVENTTYPE, CHROMOSOME, SPANSTART, SPANEND, NSETS";  
$fieldNames{CNCs} = "CNCID, EVENTID, CHROMOSOME, START_, END_, COPYNUMBER, NORMALIZEDCOVERAGE"; 
#########################################################################
  
               
#########################################################################
#table partition definitions
# = [$chromPartField, $typeSubPartField, $typeSubPartsRef]
#------------------------------------------------------------------------              
$partitions{Reads} = []; #not partitioned
$partitions{LU} = [];
$partitions{Pairs} = ["CHROMOSOME1", "PAIRTYPE", $types{Pairs}];
#must purge ? and . out of Disc Types for subpartition, since will cause Oracle name error
my %discSubParts = (None=>0, A=>1, C=>2, G=>3, T=>4, N=>5, Missing=>6, Extra=>7);
$partitions{Discs} = ["CHROMOSOME", "DISCREPANCYTYPE", \%discSubParts];
$partitions{Stats} = [];
$partitions{Frags} = ["CHROMOSOME1", "FRAGMENTTYPE", $types{Frags}];
$partitions{FMap} = ["CHROMOSOME"];
$partitions{RMap} = $partitions{FMap};
$partitions{ZMap} = $partitions{FMap};
#$partitions{FRatio} = ["CHROMOSOME"];
$partitions{Sets} = ["CHROMOSOME1", "SETTYPE", $types{Sets}];
$partitions{IDs} = [];
$partitions{CDS} = [];
$partitions{NCDS} = [];
$partitions{h} = [];
$partitions{rmsk} = ["CHROMOSOME"];
$partitions{CNVs} = ["CHROMOSOME"];
$partitions{Array} = ["CHROMOSOME"];
$partitions{NMA} = ["CHROMOSOME"];
$partitions{LD} = ["CHROMOSOME"];
$partitions{AF} = ["CHROMOSOME"];
$partitions{Events} = [];
$partitions{CNCs} = [];
#########################################################################


#########################################################################
#genetic code
#------------------------------------------------------------------------ 
our %geneticCode = (
        # * - Stop
        'TAA'=>'x','TAG'=>'x','TGA'=>'x',
        # A - Alanine
        'GCT'=>'A','GCC'=>'A','GCA'=>'A','GCG'=>'A',
        # C - Cysteine
        'TGT'=>'C','TGC'=>'C',
        # D - Aspartic Acid
        'GAT'=>'D','GAC'=>'D',
        # E - Glutamic Acid
        'GAA'=>'E','GAG'=>'E',
        # F - Phenylalanine
        'TTT'=>'F','TTC'=>'F',
        # G - Glycine
        'GGT'=>'G','GGC'=>'G','GGA'=>'G','GGG'=>'G',
        # H - Histidine
        'CAT'=>'H','CAC'=>'H',
        # I - Isoleucine
        'ATT'=>'I','ATC'=>'I','ATA'=>'I',
        # K - Lysine
        'AAA'=>'K','AAG'=>'K',
        # L - Leucine
        'CTT'=>'L','CTC'=>'L','CTA'=>'L','CTG'=>'L',
        'TTA'=>'L','TTG'=>'L',
        # M - Methionine
        'ATG'=>'M',
        # N - Asparagine
        'AAT'=>'N','AAC'=>'N',
        # P - Proline
        'CCT'=>'P','CCC'=>'P','CCA'=>'P','CCG'=>'P',
        # Q - Glutamine
        'CAA'=>'Q','CAG'=>'Q',
        # R - Arginine
        'CGT'=>'R','CGC'=>'R','CGA'=>'R','CGG'=>'R',
        'AGA'=>'R','AGG'=>'R',
        # S - Serine
        'TCT'=>'S','TCC'=>'S','TCA'=>'S','TCG'=>'S',
        'AGT'=>'S','AGC'=>'S',
        # T - Threonine
        'ACT'=>'T','ACC'=>'T','ACA'=>'T','ACG'=>'T',
        # V - Valine
        'GTT'=>'V','GTC'=>'V','GTA'=>'V','GTG'=>'V',
        # W - Tryptophan
        'TGG'=>'W',
        # Y - Tyrosine
        'TAT'=>'Y','TAC'=>'Y',
);
#########################################################################


1;

