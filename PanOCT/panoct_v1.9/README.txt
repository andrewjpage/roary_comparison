Copy (C) 2011-2012  The J. Craig Venter Institute (JCVI).  All rights reserved

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

PATCHES
-------

Suggested updates should be directed to Derrick Fouts (dfouts@jcvi.org) for consideration.

INTRODUCTION
------------

PanOCT, Pan-genome Ortholog Clustering Tool, a heuristic computer program, was created as a tool for pan-genomic analysis of closely related prokaryotic species or strains.  For more information please visit the PanOCT website at http://panoct.sourceforge.net/

PanOCT was written by: 

Derrick E. Fouts, Ph.D.
Assistant Profesor
Genomic Medicine

and 

Granger Sutton, Ph.D.
Professor
Informatics

The J. Craig Venter Institute (JCVI)
9704 Medical Center Drive
Rockville, MD  20850
(301) 795-7874
dfouts@jcvi.org

SYSTEM REQUIREMENTS
-------------------

The programs should run on all Unix platforms.  It has been tested on CentOS Linux and Mac OS X 10.6 operating systems.  It has not been tested on all Unix platforms.

Memory usage guide (based on 4 Mbp bacterial genomes, estimations not guaranteed):
1 GB RAM	1-5 genomes
2 GB RAM	1-8 genomes
4 GB RAM	1-12 genomes
8 GB RAM	1-18 genomes
14 GB RAM	1-25 genomes

SOFTWARE REQUIREMENTS/DEPENDENCIES
----------------------------------

PanOCT.pl requires the following programs or packages for full functionality:

PERL version 5.10 or later (http://www.perl.org)
NCBI BLASTALL version 2.2.10 or later [Linux] (ftp://ftp.ncbi.nih.gov/blast/executables/release/) 
WUBLAST 2.0 (now AB-BLAST at http://blast.advbiocomp.com/)
Getopt::Std [should be installed with PERL]
Data::Dumper PERL module (available from http://search.cpan.org/~smueller/Data-Dumper-2.131/Dumper.pm).  Used for debugging.

INCLUDED IN DISTRIBUTION
------------------------

-PERL SCRIPTS AND MODULES- 

~/panoct/bin/PanOCT.pl:  The main PERL script for finding orthologous proteins.

INSTALLATION
------------

First, place the distribution tarball to your home directory (~/)

Second, uncompress the distribution tarball by typing:

% tar -xvzf panoct_v1.6.tar.gz

REQUIRED INPUT FILES
--------------------

1) WU-BLAST or NCBI (-m 8 or 9 option) btab input file
2) A text file containing unique genome identifiers, one identifier per line, to determine which genome is to be treated as the reference genome in the output files and which genomes to include in the analysis.  The genome identifier can be associated with specific proteins in two ways:  a) by placing the identifier after the protein identifier (e.g. NT08AB0001-GENOME_IDENTIFIER) or b) in the gene attribute file.
3) The genome attribute file is a tab-delimited file containing the following data:  contig id, protein identifier (e.g. locus), 5Õ coordinate, 3Õ coordinate, annotation, and genome identifier. 
4) The protein fasta file used in the all-versus-all BLASTP searches.  The protein fasta file is used by PanOCT to calculate the length of each protein, which is necessary in order to compute the BSR. 

INVOCATION
----------

Usage: PanOCT.pl <options>
Example: PanOCT.pl -t example_blast.txt -f example_tags.txt -g example.gene_att -P example.pep -i 20 -F 1.33 -N Y -M Y -H Y -G Y
Version: ver1.9
 Switch: -h for help
 Option:
     -b: base directory path [DEFAULT = PWD]
     -p: path to btab file (DEFAULT = base directory)
     -t: name of btab (wublast-style or ncbi -m8 or -m9) input file [REQUIRED]
     -f: file containing unique genome identifier tags [REQUIRED]
     -g: gene attribute file (asmbl_id<tab>protein_identifier<tab>end5<tab>end3<tab>annotation<tab>genome_tag)
     -P: name of concatinated .pep file [REQUIRED to calc protein lengths]
     -Q: path to .pep file [DEFAULT = base directory)
     -i: aa % identity cut-off [DEFAULT = 35.0]
     -E: E-value [DEFAULT = 0.00001]
     -L: Minimum % match length [DEFAULT = 1]
     -V: Want to create an ortholog matchtable?  (y)es or (n)o [DEFAULT = YES]
     -N: Want to create a normalized BLAST score file?  (y)es or (n)o [DEFAULT = NO]
     -M: Want to create microarray-like data for normalized BLAST scores?  (y)es or (n)o [DEFAULT = NO]
     -H: Want to create a table of hits (y)es or (n)o [DEFAULT = NO]
     -G: Want to create score histograms (y)es or (n)o [DEFAULT = NO]
     -F: Deprecate shorter protein fragments when protein is split due to frameshift or other reason
         Takes an argument between 1.0 and 2.0 as a length ratio test - recommended value is 1.33
     -a: Number of amino acids at the beginning or end of a match that can be missing and still be
         considerd a full length match - must be between 0 and 100 - default is 20
     -s: Number of blast matches needed to confirm a protein fragment/frameshift - default 1
     -d: takes no argument, if present overrides default and does not output deprecated fragments
     -D: DEBUG MODE (DEFAULT = 0) 
 Output:  All stored within a subdirectory of the current working directory ($PWD)
          1) panoct_matchtable.txt:  a tab-delimited file containing PanOCT clusters, one cluster per line, one protein identifier per column.  The first column is the reference 
                                     genome and all subsequent columns are the remaining genomes in the order specified in the genome identifier "tags" file (specified with option f).
                                     e.g. NT08AB0001	NT16AB0001	NT17AB0001	NT20ABA0020
          2) panoct_matchtable_id.txt:  a tab-delimited file similar to panoct_matchtable.txt, but also containing the percent identity of each target protein in parentheses.
                                        e.g. NT08AB0001      NT16AB0001 (99.78%)     NT17AB0001 (99.26%)     NT20ABA0020 (99.78%)
          3) panoct_BSR.txt:  a tab-delimited file containing reference protein identifier, annotation, and BSR of each ortholog from target genomes where 1 = perfect match and 0 = no match
                              e.g. NT08AB0001      chromosomal replication initiator protein DnaA  1       0.863285556780595       1
          4) panoct_micro.txt:  a tab-delimited file containing reference protein identifier, annotation, and BSR rescaled such that 1 = perfect match, 100 = no match.
                                this was geared toward visualizing ortholog data using microarray software tools.
                                e.g. NT08AB0001	chromosomal replication initiator protein DnaA	1.00	14.53	1.00
          5) panoct_hits.txt:  a tab-delimited file containing reference protein identifier followed by the protein identifier and annotation of each ortholog from target genomes.
                               e.g. NT08AB0001 [chromosomal replication initiator protein DnaA]     NT16AB0001 [chromosomal replication initiator protein DnaA]     NT17AB0001 [chromosomal replication initiator protein DnaA]     NT20ABA0020 [chromosomal replication initiator protein DnaA]

          6) panoct_id.txt:  a tab-delimited file containing reference protein identifier,  annotation and percent identities to orthologs in each target genome.
                             e.g. NT08AB0001	chromosomal replication initiator protein DnaA	99.78	99.26	99.78
          7) panoct_frameshifts.txt:  a tab-delimited file containing proteins that are likely split due to frame-shifts.  It is organized by genome, assembly/contig and protein identifer 
                                      e.g. >genome ntab08
                                           >asmbl_id 1
                                           NT08AB3019	NT08AB3018	NT08AB3020
          8) panoct_report.txt:  a file containing runtime parameters used (e-value, %id, match length, and blast file used)
 
          
EXAMPLE DATA
------------

Sample data can be found in ~/panoct_v1.6/example_dir

