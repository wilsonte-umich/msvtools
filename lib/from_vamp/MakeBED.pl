#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#MakeBED.pl creates a single BED file of found sets
#for uploading to the UCSC Genome Browser
#########################################################################

use vars(qw(%param %types %fields %refSeqs $uid %reverseRefSeqs));

my %brightColors;
my %darkColors;
my %shortTypes;

my $setSizeLimit = 1000000;

sub makeSetsBED{
    my ($sample) = @_;
    initializeHashes();
    my $setsTable = getTableName('Sets', $sample);
    my %dirs;
    getDirectories($sample, \%dirs);
    (-d $dirs{sample}) or mkdir($dirs{sample});
    my $bedFile = "$dirs{sample}/$sample\_SetsBED.txt";
    open my $bedFH, ">", $bedFile;
    my $track = "$sample\_Sets";        
    print $bedFH "track name=$track ";
    print $bedFH "description=$track ";
    print $bedFH "visibility=pack ";
    print $bedFH "itemRgb=\"On\" \n";
    runSQL("SELECT SETTYPE, NFRAGSTOTAL, NFRAGSSAMPLE, EVENTMEAN, CHROMOSOME1, OVERLAPSTART, OVERLAPEND
           FROM $setsTable
           WHERE (SETTYPE = $types{Sets}{Deletion}
                OR SETTYPE = $types{Sets}{Insertion}
                OR SETTYPE = $types{Sets}{Inversion}
                OR SETTYPE = $types{Sets}{Duplication})
                AND (OVERLAPEND - OVERLAPSTART + 1) <= $setSizeLimit
           ORDER BY OVERLAPSTART, OVERLAPEND",
           \my($setType, $nFragTotal, $nFragsSample, $eventMean, $chrom, $oStart, $oEnd));
    while (fetchRow()){
        my $color = $darkColors{$setType};
        ($nFragTotal = $nFragsSample) and $color = $brightColors{$setType};
        my $name = "$shortTypes{$setType}($eventMean)";
        print $bedFH join(" ", ($reverseRefSeqs{$param{refSeqBase}}{$chrom},
                                $oStart, $oEnd,
                                $name, 0, '.', 0, 0, $color))."\n";  
    }
   close $bedFH;
   if ($param{refSeq} eq 'hg19'){hg19tohg18($bedFile)}
}

sub initializeHashes{
    %brightColors = (
        $types{Sets}{Deletion} => '0,0,255', #blue
        $types{Sets}{Insertion} => '0,0,0',  #black
        $types{Sets}{Inversion} => '255,0,0',  #red
        $types{Sets}{Duplication} => '0,255,0',  #green
    );
    %darkColors = (
        $types{Sets}{Deletion} => '0,0,125', #blue
        $types{Sets}{Insertion} => '0,0,0',  #black
        $types{Sets}{Inversion} => '125,0,0',  #red
        $types{Sets}{Duplication} => '0,125,0',  #green
    );
    %shortTypes = (
        $types{Sets}{Deletion} => 'Del', #blue
        $types{Sets}{Insertion} => 'Ins',  #black
        $types{Sets}{Inversion} => 'Inv',  #red
        $types{Sets}{Duplication} => 'Dup',  #green        
    );
}

sub commify {
    my $input = shift;
    $input =~ s/(?<=\d)(?=(?:\d\d\d)+$)/,/g;  #regex to add commas to numbers
    return $input;
}

sub hg19tohg18{
    my ($bedFilehg19) = @_;
    my $bedFilehg19Stripped = "$bedFilehg19.stripped";
    open my $bF19H, "<", $bedFilehg19;    
    open my $bF19SH, ">", $bedFilehg19Stripped;
    my $line1 = <$bF19H>;
    while (<$bF19H>){print $bF19SH $_}
    close $bF19H;
    close $bF19SH;
    my $lODir = "$param{vampPath}/liftOver";
    my $bedFilehg18TMP = "$bedFilehg19.hg18.tmp";
    system("$lODir/liftOver $bedFilehg19Stripped $lODir/hg19ToHg18.over.chain $bedFilehg18TMP $bedFilehg18TMP.unmapped");
    my $bedFilehg18 = "$bedFilehg19.hg18";
    open my $bF18TH, "<", $bedFilehg18TMP;
    open my $bF18H, ">", $bedFilehg18;
    print $bF18H $line1;
    while (<$bF18TH>){print $bF18H $_}
    close $bF18H;
    close $bF18TH;
}


#1=Chromosome
#2=Start
#3=End
#4. name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
#5. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
#     shade 	  	  	  	  	  	  	  	  	 
#     score in range   	= 166 	167-277 	278-388 	389-499 	500-611 	612-722 	723-833 	834-944 	= 945
#6. strand - Defines the strand - either '+' or '-'.
#7. thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
#8. thickEnd
#9. itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
#10. blockCount - The number of blocks (exons) in the BED line.
#11. blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
#12. blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.

1;






