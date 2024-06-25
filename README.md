
--------------------------------------------------------------------------------
OVERVIEW
--------------------------------------------------------------------------------

msvtools explores microarray data for genomic structural variants, specifically
CNVs. msvtools implements the analysis logic of VAMP as described here:

  Arlt MF, Ozdemir AC, Birkeland SR, Lyons RH Jr, Glover TW, Wilson TE.
  Comparison of constitutional and replication stress-induced genome structural
  variation by SNP array and mate-pair sequencing. Genetics. 2011 187:675-83.

--------------------------------------------------------------------------------
USAGE PHASES
-------------------------------------------------------------------------------- 

msvtools is a suite of scripts that manage data flow through a series of
open-source third-party programs and applications. 

Command line tools run under a Linux or compatible environment.
  
--------------------------------------------------------------------------------
REQUIREMENTS
--------------------------------------------------------------------------------

### Analysis Phase

System prerequisites for running msvtools in the analysis phase are:

Perl:      http://www.perl.org/

R:         http://www.r-project.org/ (plus package mixtools for the `train` command)
           

in addition to other standard Linux system utilities. The installation 
utility will apprise you if you are missing a prerequisite. Note that msvtools 
cannot be run on Windows in the analysis phase.

Folder 'msvtools/q-scripts' carries q master scripts 
for running msvtools analysis using the 
[q-pipeline-manager](https://github.com/wilsonte-umich/q-pipeline-manager).

### Visualization Phase

System prerequisites for running msvtools in the visualization phase are:

Perl:      http://www.perl.org/

R:         http://www.r-project.org/

in addition to specific perl packages listed in INSTALL. 

--------------------------------------------------------------------------------
INSTALLATION
--------------------------------------------------------------------------------

Read the INSTALL file to learn how to use the 'configure.pl' installation script.

--------------------------------------------------------------------------------
ANALYSIS COMMAND SUMMARY
--------------------------------------------------------------------------------

msvtools encompasses a series of commands called as follows:

```sh
    msvtools <command> [options]
```

Use:

```sh
    msvtools -h/--help
    msvtools <command> -h/--help
```

for more information on individual commands and their options.

--------------------------------------------------------------------------------
INPUT DATA
--------------------------------------------------------------------------------




--------------------------------------------------------------------------------
USAGE EXAMPLES
--------------------------------------------------------------------------------



