#!/usr/bin/perl
use strict;
use warnings;

#@RG     ID:3034_7       PL:SLX  LB:129P2_SLX_500_HC_1   PI:500  DS:129P2_Mouse_Genome   SM:129P2        CN:SC
my (%readGroups, %fragSizes, %libraries);
while (my $line = <STDIN>) {
    $line =~ m/^\@RG/ or next;
    my @line = split("\t", $line);
    my $library = getSamTagValue(\@line, 'LB');
    $library or $library = 'NA';
    my $readGroup = getSamTagValue(\@line, 'ID');
    push @{$readGroups{$library}}, $readGroup;
    $fragSizes{$library} = getSamTagValue(\@line, 'PI');
    $fragSizes{$library} or $fragSizes{$library} = 'NA';
    $libraries{$readGroup} = $library;
}
foreach my $library(sort {$a cmp $b} keys %readGroups){
    print "library $library (nominal fragment size = $fragSizes{$library})\n";
    print "  associated read groups:\n";
    foreach my $readGroup(sort {$a cmp $b} @{$readGroups{$library}}){
        print "    $readGroup\n";
    }
}
sub getSamTagValue {
    my ($line, $tag) = @_;
    foreach my $tagField(@$line){ 
        $tagField =~ m/$tag:(.*)/ or next;
        return $1;   
    }
    return undef;
}

1;

