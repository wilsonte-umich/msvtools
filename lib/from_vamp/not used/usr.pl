#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#usr.pl loads any scripts found in 'XXX/vamp/bin/usr' where XXX is either the
#home directory or vampPath, with the home directory taking precedence
#if a script of the same name is found in both locations
#########################################################################

use vars(qw(%param));

my %usrScripts;
foreach my $path("$ENV{HOME}/vamp","$param{vampPath}"){
    $path .= "/bin/usr";
    -d $path or next;
    foreach my $script(<$path/*.pl>){ 
        $usrScripts{$script} or $usrScripts{$script} = require $script;
    }
}

#require "$param{vampPath}/bin/usr/sim.pl";
#require "$param{vampPath}/bin/usr/AssembleEvents.pl";

1;

