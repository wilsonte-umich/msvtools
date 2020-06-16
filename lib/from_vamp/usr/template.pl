#!/usr/bin/perl
use strict;
use warnings;

####################################################
#template.pl provides a simple template for you to use
#to create new subs to be run under the VAMP envelope.
#The variables provided by 'use vars', provide access
#to most VAMP parameters, etc.  If this script is included
#in Directory <vampPath>/bin/usr, it will be loaded whenever
#VAMP is run.  Thus, any subs in this script also have access
#to any of the other scripts, subs or data structures 
#provided by VAMP (especially database access).
#Subs contained in scripts under <vampPath>/bin/usr
#are called using the special VAMP command 'run' as follows:
#   run command parameter1 [parameter2 ...]
#For example, the sub in this template could be run as:
#   run printMyText Hello_World
#which would print the text 'Hello_World' to stdout.
####################################################

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

sub printMyText{
    my ($toPrint) = @_;
    print "$toPrint\n";
}

1;

