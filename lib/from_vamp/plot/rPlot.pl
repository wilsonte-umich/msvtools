#!/usr/bin/perl
use strict;
use warnings;

#wholly generic functions for creating R plots of VAMP data

use vars(qw(%param %command));

#call first, sets up common plot defaults
sub intializeRPlot { 
    my ($dataFile, $jpgExtension, $main,
        $xlab, $xmin, $xmax, $xlog,
        $ylab, $ymin, $ymax, $ylog) = @_;
        
    #adjust inputs
    if($jpgExtension){ $jpgExtension = ".$jpgExtension" } else { $jpgExtension = "" }  
    $main or $main = "";   
    my $log = "";
    $xlog and $log = "x";
    $ylog and $log = $log."y"; 
     
    #set up default plot parameters
    my @r;
    push @r, "dataFile  <- '$dataFile'",
             "data      <- read.table(dataFile,header=TRUE,sep=',')",
             "jpgFile   <- paste(dataFile,'$jpgExtension.jpg',sep='')",
             "main      <- '$main'",
             "width     <- 7",
             "height    <- 7",
             "units     <- 'in'",
             "pointsize <- 12",
             "res       <- 300",
             "pch       <- 20",
             "cex       <- 0.20",
             "xlab      <- '$xlab'",
             "xmin      <- $xmin",
             "xmax      <- $xmax",
             "ylab      <- '$ylab'",
             "ymin      <- $ymin",
             "ymax      <- $ymax",
             "log       <- '$log'",
             "bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)",
             "plot(1,1,pch='',log=log,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=xlab,ylab=ylab,main=main)";       
    return \@r;
}

#a series of intermediate calls to configure the rules and plotted data series
sub addRDiagonalRule { 
    my ($r, $x1, $y1, $x2, $y2, $col, $lty) = @_;
    defined $x1 or $x1 = 'xmin'; 
    defined $x2 or $x2 = 'xmax'; 
    defined $y1 or $y1 = 'ymin'; 
    defined $y2 or $y2 = 'ymax'; 
    $col or $col = 'black';
    $lty or $lty = 1;
    push @$r, "lines(c($x1,$x2),c($y1,$y2),col='$col',lty=$lty)";
}
sub addRHorizonatalRule {
    my ($r, $y, $col, $lty) = @_; 
    $col or $col = 'black';
    $lty or $lty = 1;
    push @$r, "lines(c(xmin,xmax),c($y,$y),col='$col',lty=$lty)";
}
sub addRVerticalRule {
    my ($r, $x, $col, $lty) = @_; 
    $col or $col = 'black';
    $lty or $lty = 1;
    push @$r, "lines(c($x,$x),c(ymin,ymax),col='$col',lty=$lty)";
}
sub addRLines {
    my ($r, $x, $y, $col, $lty) = @_;
    $col or $col = 'black';
    $lty or $lty = 1;
    push @$r, "lines(data\$$x,data\$$y,col='$col',lty=$lty)";
}
sub addRPoints {
    my ($r, $x, $y, $col, $pch, $cex) = @_; #$pch must carry quotes if needed
    $col or $col = 'black';
    $pch or $pch = 'pch';
    $cex or $cex = 'cex';
    push @$r, "points(data\$$x,data\$$y,pch=$pch,cex=$cex,col='$col')";
}
sub addRLegend {
    my ($r, $samples, $colors) = @_; 
    $samples = "'".join("','", @$samples)."'";
    $samples = "c($samples)";
    $colors = "'".join("','", @$colors)."'";
    $colors = "c($colors)";
    push @$r, "legend(xmin,ymax,$samples,fill=$colors)";
}

#call last, actually creates R script and jpg plot
sub createRPlot { 
    my ($r, $dataFile, $rFile) = @_;
    $rFile or $rFile = "$dataFile.R";
    open my $rFileH, ">", $rFile or die "could not open $rFile for writing: $!\n";
    print $rFileH join("\n", @$r)."\n";
    close $rFileH;
    system("Rscript $rFile");
}

#sub to output sql data to csv file for R handling and also storage and sharing
sub printPlotData {
    my ($sql, $fileRoot, @columns) = @_;
    runSQL($sql);
    my $file = "$param{inputPath}/$fileRoot.csv";
    open my $fileH, ">", $file or die "could not open $file for writing: $!\n";
    print $fileH join(",", @columns)."\n";
    while (my @vals = fetchRowArray()){
        fixOracleNulls(\@vals);
        print $fileH join(",", @vals)."\n";
    }
    close $fileH;
    return $file;
}
#sub to make sure nulls values print correctly for R handling
sub fixOracleNulls {
    my ($vals) = @_;
    foreach my $i(0..(scalar(@$vals)-1)){
        defined $$vals[$i] or $$vals[$i] = "";
    }                
}

1;

