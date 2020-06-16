#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#Causes VAMP to require all scripts in all folders specified by parameter requireFolders.
#At least one script in folder must contain definitions for any new callable
#parameters and commands as follows:
#  defined $param{paramName} or $param{paramName} = defaultValue;
#  defined $command{commandName} or $command{commandName} = ['multiThread', '1:00:00', 500, 0];
#########################################################################

use vars(qw(%param));

sub requireFolders {
    $param{require} or return;
    my @requires = split(",", $param{require});
    REQUIRE: foreach my $require(@requires){
        foreach my $sourceFolder("$param{vampPath}/bin/", ""){
            my $requireFolder = "$sourceFolder$require";
            requireFolders_($requireFolder) and next REQUIRE;
        }
    }
}

sub requireFolders_ {
    my ($requireFolder) = @_;
    -d $requireFolder or return undef;
    foreach my $script(<$requireFolder/*.pl>){ require $script }   
    return 1;
}

1;
