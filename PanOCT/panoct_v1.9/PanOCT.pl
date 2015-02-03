#!/usr/local/bin/perl -w
#Copy (C) 2011-2012  The J. Craig Venter Institute (JCVI).  All rights reserved
#Written by Derrick E. Fouts, Ph.D. and Granger Sutton, Ph.D.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Revision notes
# Match table and other revisions 01/07/2004
# Added NCBI -m9 input option 12/08/2010
# Moved genome identification from feat/locus_name-tag to the gene_att file 12/08/2010
my $prog = $0;
$prog =~ s/.*\///;

use strict;
use Getopt::Std;
use Data::Dumper;
getopts ('dhDb:p:t:f:i:I:E:L:M:G:H:V:P:Q:N:g:A:F:a:s:');#GGS added F, a, s, d, options F option for frameshift and a option for amount of missing amino acids to still be full length and s option for an evidence threshold for frameshifts and d option for whether to output deprecated fragments as singleton clusters

our ($opt_d, $opt_h,$opt_D,$opt_b,$opt_p,$opt_t,$opt_f,$opt_i,$opt_I,$opt_E,$opt_L,$opt_M,$opt_G,$opt_H,$opt_V,$opt_P,$opt_Q,$opt_N,$opt_g,$opt_A,$opt_F,$opt_a,$opt_s);#GGS added $opt_F, $opt_a, $opt_s, and $opt_d

my ($basedir,$btabpath,$pep_path,$att_file,$pep_file,$btabfile,$tagfile,$percentid,$fspercentid,$evalue,$min_hit_length,$microarray,$histogramfile,$hitsfile,$vennfile,$normalizefile,$dupsfile,$DEBUG,$frameshiftfile,$frameshift_length_ratio, $max_missing_aa, $frames_link_thresh, $output_fragments);#GGS added frameshift variables
#GGS $frameshiftfile is a boolean for detecting and outputting frameshifted or truncated protein fragments true if -F
#GGS -F $frameshift_length_ratio is a value between 1 and 2 (default 1.33) used to determine that we are detecting adjacent protein fragments rather than adjacent paralogs
#GGS -a $max_missing_aa is how many amino acids at either end of a blast match can be missing and still be considered full length
#GGS -s $frames_link_thresh is a threshold for number of blast matches indicating a frameshift to believe it
#GGS -d $output_fragments is a boolean for whether to output deprecated protein fragments as singleton clusters
#GGS -G yes outputs a histograms file
$microarray = 0; # set to zero to initialize
$hitsfile = 0;
$histogramfile = 0;
$normalizefile = 0;
$dupsfile = 0;
$vennfile = 1;

## use boolean logic:  TRUE = 1, FALSE = 0

my $version = "ver1.9";
#GGS process frameshift option
if ($opt_F) {
    $frameshiftfile = 1;
    if (($opt_F < 2) && ($opt_F > 1)) {#GGS ratio must be strictly between 1 and 2
	$frameshift_length_ratio = $opt_F;
    } else {
	&option_help;
    }
} else {
    $frameshiftfile = 0;
}
#GGS set max_missing_aa
if ($opt_a) {
    if (($opt_a > 0) && ($opt_a <= 100)) {#GGS must be between 0 and 100
	$max_missing_aa = $opt_a;
    } else {
	&option_help;
    }
} else {
    $max_missing_aa = 20;#default
}
#GGS set frameshift evidence threshold
if ($opt_s) {
    if ($opt_s > 0) {#GGS must be > 0
	$frames_link_thresh = $opt_s;
    } else {
	&option_help;
    }
} else {
    $frames_link_thresh = 1;#default
}
#GGS set protein fragment output option
if ($opt_d) {
    $output_fragments = 0;#do not output deprecated fragments as singleton clusters
} else {
    $output_fragments = 1;#default is to output deprecated fragments as singleton clusters
}
if ($opt_h) { &option_help; } # quit with help menu
if ($opt_b) {$basedir = $opt_b;} else { $basedir = $ENV{'PWD'}; } # if no value for option b (base or working directory) set it to current directory
if ($opt_p) {$btabpath = $opt_p;} else { $btabpath = $basedir; } # if no value given for path to btab file, make default = base directory
if ($opt_Q) {$pep_path = $opt_Q;} else { $pep_path = $basedir; } # if no value given for path to btab file, make default = base directory
if (($opt_g) && (-s "$opt_g")) {$att_file = $opt_g;} else { print STDERR "Error with -g\n"; &option_help; } # if no value for option g (name of gene_att file), quit with help menu
if (($opt_P) && (-s "$pep_path/$opt_P")) {$pep_file = $opt_P;} else { print STDERR "Error with -P $pep_path/$opt_P\n"; &option_help; } # if no value for option P (pepfile), quit with help menu
if (($opt_t) && (-s "$btabpath/$opt_t")) {$btabfile = $opt_t;} else { print STDERR "Error with -t $btabpath/$opt_t\n"; &option_help; } # if no value for option t (name of btab file), quit with help menu
if (($opt_f) && (-s "$opt_f")) {$tagfile = $opt_f;} else { print STDERR "Error with -f"; &option_help; }  # must supply the file containing the list of unique genome names to analyze
if ($opt_i) {$percentid = $opt_i;} else { $percentid = 35.0; } # the minimum cutoff to use blast scores for possible matches
if ($opt_I) {$fspercentid = $opt_I;} else { $fspercentid = 50.0; } # the minimum cutoff to use blast scores for frameshift detection
if ($opt_E) {$evalue = $opt_E;} else { $evalue = 0.00001; } # if no E-value cut-off given, make default 0.00001
if ($opt_L) {$min_hit_length = $opt_L;} else { $min_hit_length = 1; }
if ($opt_M) {
    if (($opt_M eq "Y") or ($opt_M eq "y") or ($opt_M eq "Yes") or ($opt_M eq "yes") or ($opt_M eq "YES"))  {$microarray = 1;} 
    elsif (($opt_M eq "N") or ($opt_M eq "n") or ($opt_M eq "No") or ($opt_M eq "no") or ($opt_M eq "NO")) {$microarray = 0;}
    else { $microarray = 0; }  # Want to create microarray-like data of normalized BLAST scores  0 = NO, 1 = YES [DEFAULT = NO]
}
if ($opt_G) {
    if (($opt_G eq "Y") or ($opt_G eq "y") or ($opt_G eq "Yes") or ($opt_G eq "yes") or ($opt_G eq "YES")) {$histogramfile = 1;} 
    elsif (($opt_G eq "N") or ($opt_G eq "n") or ($opt_G eq "No") or ($opt_G eq "no") or ($opt_G eq "NO")) {$histogramfile = 0;}
    else { $histogramfile = 0; }  # Want to create histograms  0 = NO, 1 = YES [DEFAULT = NO]
}
if ($opt_H) {
    if (($opt_H eq "Y") or ($opt_H eq "y") or ($opt_H eq "Yes") or ($opt_H eq "yes") or ($opt_H eq "YES")) {$hitsfile = 1;} 
    elsif (($opt_H eq "N") or ($opt_H eq "n") or ($opt_H eq "No") or ($opt_H eq "no") or ($opt_H eq "NO")) {$hitsfile = 0;}
    else { $hitsfile = 0; }  # Want to create table of hits  0 = NO, 1 = YES [DEFAULT = NO]
}
if ($opt_V) {
    if (($opt_V eq "Y") or ($opt_V eq "y") or ($opt_V eq "Yes") or ($opt_V eq "yes") or ($opt_V eq "YES")) {$vennfile = 1;}
    elsif (($opt_V eq "N") or ($opt_V eq "n") or ($opt_V eq "No") or ($opt_V eq "no") or ($opt_V eq "NO")) {$vennfile = 0;}
    else { $vennfile = 1; } # Want to create match table file?  0 = NO, 1 = YES [DEFAULT = YES]
}
if ($opt_N) {
    if (($opt_N eq "Y") or ($opt_N eq "y") or ($opt_N eq "Yes") or ($opt_N eq "yes") or ($opt_N eq "YES")) {$normalizefile = 1;}
    elsif (($opt_N eq "N") or ($opt_N eq "n") or ($opt_N eq "No") or ($opt_N eq "no") or ($opt_N eq "NO")) {$normalizefile = 0;}
    else { $normalizefile = 0; } # Want to create normalized BLAST score file?  0 = NO, 1 = YES [DEFAULT = NO]
}
if  ( $opt_A ) { 
    if(($opt_A eq "Y") or ($opt_A eq "y") or ($opt_A eq "Yes") or ($opt_A eq "yes") or ($opt_A eq "YES")) {$dupsfile = 1;}
    elsif (($opt_A eq "N") or ($opt_A eq "n") or ($opt_A eq "No") or ($opt_A eq "no") or ($opt_A eq "NO")) {$dupsfile = 0;}
    else { $dupsfile = 0; } # Want to create a list of paralogs?  0 = NO, 1 = YES [DEFAULT = NO]
}
if ($opt_D) {$DEBUG = 1;} else { $DEBUG = 0; } # Debug mode.

my %query_evalue_cutoff = ();# Key = query id, Value = query specific evalue cutoff for low scoring queries
my %paralog_cutoff = ();    # Key1 = genome tag, Key2 = genome tag, Value = [0,1] normalized blast score cutoff for paralogs between these two genomes
my %genome_hash = ();       # Key1 = genome tag, Key2 = asmbl_id, Value = array of protein ids sorted by end5 and then end3
my %feat_hash = ();         # Key = feat_name Key2 = struct members with their values
my %relationship_hash = (); # Key1 = Query ($a[0]), Key2 = Subject (hit), Key3 = struct members with their values
my %Qbytaghash = ();        # Key1 = Query ($a[0]), Key2 = Subject (hit), Key3 = tag, Key4 = same as relationship_hash Key3 (actually is the same pointer)
my %orf_counter = ();       # Key = genome tag, Value = orf counts
my %Tagbyfeatnamehash = (); # Key1 = genome tag, Key2 = feat_name (query)
my %TagByPointer = ();      # Key = feat_name-tag, value = array position location
my @tag_array = ();         # array of genome tags
my %Paralogs = ();         # Key1 = Query tag, Key2 = Query orf, Value = subject orf (self genome hit)
my %FeatnameLookupTag_hash = ();# Key = feat_name, Value = tag.  Needed for input lacking tagged feat_names
my %LookupTagIndex = ();    # Key = genome tag, Value = index in tag_array of associated genome tag
my %AssemblyLookup_hash = ();# Generated to store the asmbl_id (value) of each feat_name (key) so we can lookup the correct array during synteny searches
my %clusters = ();       # Key = feat_name, Value = array of (hashes of feat_name and genome_tag) that are in a cluster with the feat_name key
my $minSynMatches = 2;    # Hard-code the minimum number of matches to be consider syntenous
my $min_synteny_threshold = 0;  # minimum synteny score to trust
my $genome_number;          # number of genomes defined in tagfile
my $outprefix = "$basedir/panoct";

#%relationship_hash stores information for the blast matches
###key1 = query
###key2 = subject
###value -> id = percent identity
###         eval = e-value
###         score = BLAST bits score
###         BSR = BLAST Score Ratio # added 12/13/2010 so that the true top matches were chosen, not the Blast score
###         best = 1 if best blast match in the subject genome, otherwise 0
###         bibest = 1 if best bidirectional blast match in the subject genome, otherwise 0
###         min_query = smaller match coordinate on the query protein
###         max_query = larger match coordinate on the query protein
###         min_sub = smaller match coordinate on the subject protein
###         max_sub = larger match coordinate on the subject protein
###         clique_top = number in the clique at the top of blast matches
###         clique_all = number in the clique for all blast matches


sub get_tags {  # obtain list of genomes to compare
   
    my %temp_hash = ();
    $genome_number = 0;     # total number of genomes to be processed

    open (INFILE, "<$basedir/$tagfile");
    while (<INFILE>)  {
	chomp;
	if (length($_) > 0) { #((exists $genome_hash{$b[1]})
            if ($temp_hash{$_})  {
               die ("ERROR:  You have more than one occurance of $_ in $basedir/$tagfile!\n");
            }
            else  {
		$temp_hash{$_} = 1; # used to be genome_hash
		push (@tag_array, $_); # populate the tag_array in the order of the tagfile (1st tag is the reference tag)
		$genome_number++;
		if ($genome_number == 1) {
		    print STDERR "$_ [Reference]\n";
		}
		else {
		    print STDERR "  $_\n";
		}
	    }
	}  
    }
    close(INFILE);

    my $index = 0;
    foreach my $tag (@tag_array) {
	$LookupTagIndex{$tag} = $index++
    }
}

sub get_protein_info { #formerly get_fasta_headers.  Added calculation of protein length for ncbi data 12/08/10

  my @line = ();
  my $id;
  my $title = "";
  my $sequence = "";
  my $length = "";

  unless (open (PEPFILE, "<$pep_path/$pep_file") )  {
    die ("can't open file $pep_path/$pep_file.\n");
  }
  my ($save_input_separator) = $/;
  $/="\n>";
  while (<PEPFILE>) {
    ($title,$sequence) = /^>?\s*(.*)\n([^>]+)>?/; # split the header line and sequence (very cool)
    @line = split(/\s+/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @fasta
    $id = $line[0]; # unique orf identifier is in column 0, com_name is in rest
    $id =~ s/>//;
    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet letters
    $length = length($sequence);
    $feat_hash{$id}->{'header'} = join(' ', @line[1..$#line]); # put the identifier into the hash as the "key" and the header as value "header" (joining all columns after first space/tab)
    $feat_hash{$id}->{'length'}= $length; # only store the length, not the sequence
    #print STDERR "$id ($feat_hash{$id}->{'header'}) = $feat_hash{$id}->{'length'}\n";
    $title = ""; # clear the title for the next round.
    $sequence = ""; #clear out the sequence for the next round.
  }
  $/ = $save_input_separator; # restore the input separator
  close (PEPFILE);
  return;
}


sub get_gene_att {

    my $pos = "";
    my $tag = "";
    my $end5 = "";
    my $end3 = "";
    my $asmbl_id = "";
    my $feat_name = "";
    my $anno = "";

    unless (open (ATTFILE, "<$basedir/$att_file") )  {
	die ("ERROR: can not open file $basedir/$att_file.\n");
    }
    while (<ATTFILE>) {
	my @att_line = ();
	chomp;
	@att_line = split(/\t/, $_);  # split the scalar $line on tab
	$asmbl_id = $att_line[0];
	$feat_name = $att_line[1];
	$end5 = $att_line[2];
	$end3 = $att_line[3];
	$anno = $att_line[4];
	#tag if present is in $a[5]
	if ($feat_name =~ /-/) { # if tag provided within the feat_name or locus_name, then use it
	    my @split_name = split(/-/, $feat_name); # if tag is appended like this feat_name-tag
	    $tag = $split_name[$#split_name];
	    $feat_name = join("-", @split_name[0..($#split_name - 1)]);
	} elsif (defined($att_line[5])) {
	    $tag = $att_line[5];
	}
	else { # quit on error - no tags provided for each protein
	    print STDERR "ERROR: user must provide genome tags on feat/locus_name or in column 6 of the gene_att file.\n";
	    exit(1);
	}
	if (!defined $LookupTagIndex{$tag}) {
	    print STDERR "ERROR: $tag is a genome tag in the attribute file but not in the genome tag file!\n";
	    exit(1);
	}
	$feat_hash{$feat_name}->{'full'} = 0; # initialize number of full length blast matches to 0
	$FeatnameLookupTag_hash{$feat_name} = $tag;
	push (@{ $genome_hash{$tag}{$asmbl_id} }, $feat_name);
	$AssemblyLookup_hash{$feat_name} = $asmbl_id; # new 01/12/2011
	$feat_hash{$feat_name}->{'end5'} = $end5;
	$feat_hash{$feat_name}->{'end3'} = $end3;
	if ($end5 < $end3) {
	    $feat_hash{$feat_name}->{'orient'} = 1;
	} else {
	    $feat_hash{$feat_name}->{'orient'} = -1;
	}
	$feat_hash{$feat_name}->{'fragment'} = 0;
	if (!defined $feat_hash{$feat_name}->{'header'}) {
	    print STDERR "ERROR: $feat_name in gene attribute file ($att_file) was not found in protein file ($pep_file)!\n";
	    exit(1);
	}
	if ($anno) { # check if there is any annotation from gene_att file
	    $feat_hash{$feat_name}->{'header'} = $anno;
	}
	print STDERR "$feat_name $feat_hash{$feat_name}->{'header'} $feat_hash{$feat_name}->{'end5'} $feat_hash{$feat_name}->{'end3'} $feat_hash{$feat_name}->{'length'} $feat_hash{$feat_name}->{'orient'} $tag $asmbl_id\n" if ($DEBUG);
    }
    close (ATTFILE);

    my $failed = 0;
    foreach my $feat_id (keys %feat_hash) {
	if (!defined $feat_hash{$feat_id}->{'end5'}) {
	    print STDERR "ERROR: $feat_id appears in the protein file ($pep_file) but not in the gene attribute file ($att_file)!\n";
	    $failed = 1;
	}
    }
    if ($failed) {
	exit(1);
    }

    foreach $tag (keys %genome_hash) {
	foreach $asmbl_id (keys %{ $genome_hash{$tag} }) {
	    my $sort_by_end5_end3 = sub {
		my $end5a = $feat_hash{$a}->{'end5'};
		my $end5b = $feat_hash{$b}->{'end5'};
		my $end3a = $feat_hash{$a}->{'end3'};
		my $end3b = $feat_hash{$b}->{'end3'};
		
		if ($end5a <=> $end5b) {
		    return ($end5a <=> $end5b);
		} elsif ($end3a <=> $end3b) {
		    return ($end3a <=> $end3b);
		} else {
		    print STDERR "SORRY $a and $b have the same end5 ($end5a) and end3 ($end3a), please correct the gene_att file!\n";
		    die;
		}
	    };
	    
	    @{ $genome_hash{$tag}{$asmbl_id} } = sort $sort_by_end5_end3 ( @{ $genome_hash{$tag}{$asmbl_id} } );
	    $pos = 0; # reset the array position pointer to zero for every asmbl_id (eg. chrom, plasmid)
	    foreach $feat_name ( @{ $genome_hash{$tag}{$asmbl_id} } ) {
		$TagByPointer{$feat_name} = $pos;
		$pos++;
	    }
	}
    }
    return;
}

sub select_data_from_btab { # get tab-delimited BLAST results in either WUBLAST or NCBI -m8/-m9 formats
## new code - test for btab format (WUbtab or ncbi m8/m9

    my @btab_line = (); # array variable to store split btab lines
    my $btab = ""; # 1 = WU, 0 = NCBI
    my $qmatch_length = ""; # stores the query length
    my $smatch_length = ""; # stores the subject length
    my $qid = ""; # query id
    my $sid = ""; # subject id (from database)
    my $qtag; # query genome tag
    my $stag; # subject genome tag
    my $qbegin = ""; # start query
    my $qend = ""; # end query
    my $sbegin = ""; # start subject
    my $send = ""; # end subject
    my $pid = ""; # percent identity
    my $evlu = ""; # e-value
    my $score = ""; # BLAST bit score
    my $qlength = ""; # size of query protein sequence
    my $slength = ""; # size of subject (database match) protein sequence
    open (INFILE, "<$btabpath/$btabfile") || die ("ERROR: can't open file $btabfile: $!\n");

    # this while look is interrogate the file to determine which tabular output style it is
    while (<INFILE>)  {
	chomp;
	if (!/^#/) { # don't look at REMs
	    @btab_line = split(/\t/);
	    last; # we have out info, break the loop
	}
    }
    close(INFILE);

    if ($#btab_line >= "19") { # adjusted because Badgers btab.pl script only goes to e-value or perl column 19
	print STDERR "Detected WUBLAST-style btab file ...\n";
	$btab = 1; # WU
    }
    elsif ($#btab_line == "11") {
	print STDERR "NCBI BLAST (-m 8 or -m 9) option btab file detected ...\n";
	$btab = 0; # NCBI
    }
    else {
	die ("ERROR:  BLAST data must be either WUBLAST btab or NCBI BLAST -m8 or -m9 formats.\n");
    }
    ### process BLAST results ###
    open (INFILE, "<$btabpath/$btabfile");
    while (<INFILE>)  {
	chomp;
	@btab_line = split(/\t/);
        # same variables for both btab styles
	$qbegin = $btab_line[6];
	$qend = $btab_line[7];
	$sbegin = $btab_line[8];
	$send = $btab_line[9];

	if ($btab) { # WU

            # adjusted because Badgers WU btab.pl script only goes to e-value or perl column 19
	    # ========================================================
	    # btab output for WUBLAST output
	    # column number Description (for Perl), add 1 for Unix
	    # 0       Query Sequence Name
	    # 1       Date of the Analysis
	    # 2       Query Sequence Length
	    # 3       Search Method  --  Blast family application name
	    # 4       Database Name
	    # 5       Subject Sequence Name  --  Database entry name
	    # 6       Start of alignment on query (5' nucleotide match in query)
	    # 7       End of alignment on query (3' nucleotide match in query)
	    # 8       Start of alignment on subject (5' nucleotide match in db hit)
	    # 9       End of alignment on subject (3' nucleotide match in db hit)
	    # 10      % Identity 
	    # 11      % Similarity 
	    # 12      Score (bits)
	    # 13      File Offset for Beginning of Alignment
	    # 14      File Offset for End of Alignment
	    # 15      Description (annotatioon)
	    # 16      Frame  --  1 through 6, or NULL
	    # 17      Query Strand  --  Plus, Minus or NULL
	    # 18      DB sequence length
	    # 19      Expect -- expected value
	    # 20      P-Value  --  Poisson ratio
	    # ========================================================

	    # propigate variables for WUBLAST
	    $qid = $btab_line[0];
	    $sid = $btab_line[5];
	    $pid = $btab_line[10];
	    $evlu = $btab_line[19];
	    $score = $btab_line[12];
	}
	else { # NCBI
	    # ========================================================
	    # btab output from NCBI blastn (-m 8) option:
	    # column number Description (for Perl), add 1 for Unix
	    # 0      Query_id
	    # 1	     subject_id (Hit from db)
	    # 2	     % Identity
	    # 3	     length of alignment
	    # 4	     number or mismatches
	    # 5	     number of gaps
	    # 6	     start of alignment on query (5' nucleotide match in query)
	    # 7	     end of alignment on query (3' nucleotide match in query)
	    # 8	     start of alignment on subject (5' nucleotide match in db hit)
	    # 9	     end of alignment on subject (3' nucleotide match in db hit)
	    # 10     e-value
	    # 11     score (bits)
	    # ========================================================

	    # propigate variables for NCBI m8/m9
	    if ($btab_line[0] =~ /^#/) { next;} # skip the lines beginning with #
	    $qid = $btab_line[0];
	    $sid = $btab_line[1];
	    $pid = $btab_line[2];
	    $evlu = $btab_line[10];
	    $score = $btab_line[11];
	}
	### Generic processing ###
	if ($qid =~ /-/) { # if a tag in the feat_name, parse the feat_name and the tag depending on how many dashes there are
	    my @tmpSplitQ = split(/-/, $qid);
	    $qtag = $tmpSplitQ[$#tmpSplitQ];
	    $qid = join("-", @tmpSplitQ[0..($#tmpSplitQ - 1)]);
	}
	else { # no tag provided with feat_name, so get from FeatnameLookupTag_hash
	    $qtag = $FeatnameLookupTag_hash{$qid};
	}
	if ($sid =~ /-/) {
	    my @tmpSplitS = split(/-/, $sid);
	    $stag = $tmpSplitS[$#tmpSplitS];
	    $sid = join("-", @tmpSplitS[0..($#tmpSplitS - 1)]);
	}
	else { # no tag provided with feat_name, so get from FeatnameLookupTag_hash
	    $stag = $FeatnameLookupTag_hash{$sid};
	}
        if (!defined($qtag) || !defined $LookupTagIndex{$qtag}) {print STDERR "WARNING!!! $qtag is a genome tag in the btab file but not in the genome tag file skipping this entry!\n" if ($DEBUG); next;}
        if (!defined($stag) || !defined $LookupTagIndex{$stag}) {print STDERR "WARNING!!! $stag is a genome tag in the btab file but not in the genome tag file skipping this entry!\n" if ($DEBUG); next;}
        if (!defined($qid) || !defined $feat_hash{$qid}) {print STDERR "WARNING!!! $qid is a feature identifier in the btab file but not in the gene attribute file skipping this entry!\n" if ($DEBUG); next;}
        if (!defined($sid) || !defined $feat_hash{$sid}) {print STDERR "WARNING!!! $sid is a feature identifier in the btab file but not in the gene attribute file skipping this entry!\n" if ($DEBUG); next;}
	if ($qtag ne $FeatnameLookupTag_hash{$qid}) {
	    print STDERR "ERROR: the genome tag ($qtag) in the btab file is not compatible with the genome tag in the gene attribute file ($FeatnameLookupTag_hash{$qid}) for feature identifier ($qid)!\n";
	    exit(1);
	}
	if ($stag ne $FeatnameLookupTag_hash{$sid}) {
	    print STDERR "ERROR: the genome tag ($stag) in the btab file is not compatible with the genome tag in the gene attribute file ($FeatnameLookupTag_hash{$sid}) for feature identifier ($sid)!\n";
	    exit(1);
	}
	if ((defined $genome_hash{$qtag}) && (defined $genome_hash{$stag})) {
	    $qlength = $feat_hash{$qid}->{'length'};
	    $slength = $feat_hash{$sid}->{'length'};
	    $qmatch_length = abs($qbegin - $qend)/$qlength*100; #GGS shouldn't this be (($qend - $qbegin) + 1)?
	    $smatch_length = abs($sbegin - $send)/$slength*100;  # new!
	    if (!$feat_hash{$qid}->{'bit'})  {
		$feat_hash{$qid}->{'bit'} = 1; # used to determine number of orfs processed
		$orf_counter{$qtag}->{'raw'}++;  # increment the total orf counter (we will get the total # of orfs this way)
	    }
	    if (!defined $relationship_hash{$qid}) { # this assumes that the best Blast match for a query appears first in the file
		$query_evalue_cutoff{$qid} = $evlu * 10.0;
	    }
	    if ($pid < $percentid)  {next;}
	    if (($evlu > $evalue) && ($evlu > $query_evalue_cutoff{$qid}))  {next;}
	    if ($qmatch_length < $min_hit_length)  {next;} # should we have this test only for the larger of these two?
	    if ($smatch_length < $min_hit_length)  {next;} # new!
	    if (!defined $relationship_hash{$qid}{$sid})  { # this assumes that the first Blast match in the file for these two proteins is the only important one
		$Qbytaghash{$qid}{$stag}->{$sid} = $relationship_hash{$qid}{$sid} = {}; #have Qbytaghash and relationship_hash reference the same hash
		$relationship_hash{$qid}{$sid}->{'id'} = $pid;
		$relationship_hash{$qid}{$sid}->{'eval'} = $evlu;
		$relationship_hash{$qid}{$sid}->{'score'} = $score;
		$relationship_hash{$qid}{$sid}->{'min_query'} = $qbegin;
		$relationship_hash{$qid}{$sid}->{'max_query'} = $qend;
		$relationship_hash{$qid}{$sid}->{'min_sub'} = $sbegin;
		$relationship_hash{$qid}{$sid}->{'max_sub'} = $send;
		$relationship_hash{$qid}{$sid}->{'best'} = 0;
		$relationship_hash{$qid}{$sid}->{'bibest'} = 0;
		$relationship_hash{$qid}{$sid}->{'synbest'} = 0;
		$relationship_hash{$qid}{$sid}->{'synbibest'} = 0;
		print STDERR "Query: $qid X Subject $sid = $relationship_hash{$qid}{$sid}->{'id'}\n" if ($DEBUG);
		if ((($qbegin <= $max_missing_aa + 1) && ($qend + $max_missing_aa >= $feat_hash{$qid}->{'length'})) && (($sbegin <= $max_missing_aa + 1) && ($send + $max_missing_aa >= $feat_hash{$sid}->{'length'})) && ($pid >= $fspercentid)) {
		    $feat_hash{$qid}->{'full'}++; #not a protein fragment if almost a full length match
		}
	    }
	    if (!defined $Tagbyfeatnamehash{$qtag}{$qid})  {
		$Tagbyfeatnamehash{$qtag}{$qid} = 1;
		$orf_counter{$qtag}->{'used'}++;
	    }
	    ### record paralogs ###
	    if (($qtag eq $stag) && ($qid ne $sid))  {  # if hit in same genome, but not same orf, count as a paralog
		if (!defined $Paralogs{$qtag}{$sid}->{$qid})  { # prevent bidirectional entries
		    if (defined $Paralogs{$qtag}{$sid})  { # consulidate lists (ie ORF1 matches ORF2 then ORF112 matches ORF1, write ORF112 with the ORF1 list)
			$Paralogs{$qtag}{$sid}->{$qid} = 1;
		    }
		    else  {
			$Paralogs{$qtag}{$qid}->{$sid} = 1;  # add to dups hash
		    }
		}  
		next;
	    } 
	}
    }
    close (INFILE);
}

sub calc_BSR  { # subroutine to generate the average (both directions) Blast Score Ratio (BSR) and populate Qbytaghash and relationship_hash with BSR

    my $qid = ""; # query feat_name
    my $sid = ""; # subject feat_name
    my $qtag = ""; # query genome tag
    my $stag = ""; # subject tag

    foreach $qtag (@tag_array)  {  # start looping through by order in tag file (reference is first)
	foreach $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
	    foreach $stag (@tag_array)  {
		if (defined $Qbytaghash{$qid}{$stag}) { # if query protein matches anything in subject genome, lets drill through each match
		    foreach $sid (keys %{ $Qbytaghash{$qid}{$stag} })  { # need to check value of $sid here (see comment in previous subroutine)
			# calculate the Blast Score Ratio (BSR) of each relationship
			if (!defined $relationship_hash{$sid}{$qid}) {
			    $Qbytaghash{$sid}{$qtag}->{$qid} = $relationship_hash{$sid}{$qid} = {}; #have Qbytaghash and relationship_hash reference the same hash
			    $relationship_hash{$sid}{$qid}->{'score'} = 0; #set score to zero if no match in reverse direction
			    #setting these values to make sure they are not undefined when symmetric entry exists
			    $relationship_hash{$sid}{$qid}->{'id'} = $relationship_hash{$qid}{$sid}->{'id'};
			    $relationship_hash{$sid}{$qid}->{'eval'} = $relationship_hash{$qid}{$sid}->{'eval'};
			    $relationship_hash{$sid}{$qid}->{'min_query'} = $relationship_hash{$qid}{$sid}->{'min_sub'};
			    $relationship_hash{$sid}{$qid}->{'max_query'} = $relationship_hash{$qid}{$sid}->{'max_sub'};
			    $relationship_hash{$sid}{$qid}->{'min_sub'} = $relationship_hash{$qid}{$sid}->{'min_query'};
			    $relationship_hash{$sid}{$qid}->{'max_sub'} = $relationship_hash{$qid}{$sid}->{'max_querry'};
			    $relationship_hash{$sid}{$qid}->{'best'} = 0;
			    $relationship_hash{$sid}{$qid}->{'bibest'} = 0;
			    $relationship_hash{$sid}{$qid}->{'synbest'} = 0;
			    $relationship_hash{$sid}{$qid}->{'synbibest'} = 0;
			}
			if (!defined $relationship_hash{$qid}{$qid}) {
			    print STDERR "WARNING!!! no Blast self hit for qid $qid setting self hit to highest score\n";
			    $relationship_hash{$qid}{$qid}->{'score'} = 1; # just in case there are no scores
			    foreach my $tmpid (sort {$relationship_hash{$qid}{$b}->{'score'} <=> $relationship_hash{$qid}{$a}->{'score'}} keys %{ $relationship_hash{$qid} }) {
				$relationship_hash{$qid}{$qid}->{'score'} = $relationship_hash{$qid}{$tmpid}->{'score'};
				last; #quit after storing best score
			    }
			}
			if (!defined $relationship_hash{$sid}{$sid}) {
			    print STDERR "WARNING!!! no Blast self hit for sid $sid setting self hit to highest score\n";
			    $relationship_hash{$sid}{$sid}->{'score'} = 1; # just in case there are no scores
			    foreach my $tmpid (sort {$relationship_hash{$sid}{$b}->{'score'} <=> $relationship_hash{$sid}{$a}->{'score'}} keys %{ $relationship_hash{$sid} }) {
				$relationship_hash{$sid}{$sid}->{'score'} = $relationship_hash{$sid}{$tmpid}->{'score'};
				last; #quit after storing best score
			    }
			}
			$Qbytaghash{$qid}{$stag}->{$sid}->{'BSR'} = $Qbytaghash{$sid}{$qtag}{$qid}->{'BSR'} = ($relationship_hash{$qid}{$sid}->{'score'}/$relationship_hash{$qid}{$qid}->{'score'} + $relationship_hash{$sid}{$qid}->{'score'}/$relationship_hash{$sid}{$sid}->{'score'})/2 ; # calc BSR GGS need to keep $Qbytaghash symmetrical for frameshifts
		    }
		}
	    }
	}
    }
			
}

sub calc_bibest  { # subroutine to find bidrectional best blast matches and populate Qbytaghash and relationship_hash with them

    my $qid = ""; # query feat_name
    my $sid = ""; # subject feat_name
    my $qtag = ""; # query genome tag
    my $stag = ""; # subject tag

    #mark best blast hit per genome per qid
    foreach $qtag (@tag_array)  {  # start looping through by order in tag file (reference is first)
	foreach $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
	    if (!defined $Qbytaghash{$qid}) {
		next;# need to do this to prevent keys %{ $Qbytaghash{$qid}{$stag} } causing $Qbytaghash{$qid}{$stag} to be defined later
	    }
	    foreach $stag (keys %{ $Qbytaghash{$qid} })  {
		if ($qtag eq $stag) {
		    next; #skip matches within the same genome
		}
		foreach $sid (sort {$Qbytaghash{$qid}{$stag}{$b}->{'BSR'} <=> $Qbytaghash{$qid}{$stag}{$a}->{'BSR'}} keys %{ $Qbytaghash{$qid}{$stag} })  {
		    $relationship_hash{$qid}{$sid}->{'best'} = 1;
		    last; #quit after marking best blast match
		}
	    }
	}
    }

    #mark bidirectionally best blast hit per genome per qid if it exists
    foreach $qid (keys %feat_hash)  { # go through all featnames
	foreach $sid (keys %{ $relationship_hash{$qid} } )  {
		if ($relationship_hash{$qid}{$sid}->{'best'} && $relationship_hash{$sid}{$qid}->{'best'}) {
		    $relationship_hash{$qid}{$sid}->{'bibest'}  = 1;
		}
	}
    }
			
}

sub calc_score_histograms  { # subroutine to calculate and output score histograms and score cutoffs

    my $qid = ""; # query feat_name
    my $sid = ""; # subject feat_name
    my $qtag = ""; # query genome tag
    my $stag = ""; # subject tag
    my @tmp_tag_array = @tag_array;
    my %self_hist = ();

    open (OUTHISTOGRAM, ">$outprefix" . "_histograms.txt") if ($histogramfile == 1);

    foreach $qtag (@tag_array)  {  # start looping through genomes by order in tag file
	$self_hist{$qtag} = [];
	foreach my $index (0..100) { # initialize self histogram
	    $self_hist{$qtag}[$index] = 0;
	}
	my $num_qids = 0;
	foreach $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches to paralogs in the same genome
	    $num_qids++;
	    my $match_found = 0;
	    foreach $sid (sort {$Qbytaghash{$qid}{$qtag}{$b}->{'score'} <=> $Qbytaghash{$qid}{$qtag}{$a}->{'score'}} keys %{ $Qbytaghash{$qid}{$qtag} })  {
		if ($qid eq $sid) { # skip self match
		    next;
		}
		$self_hist{$qtag}[int (100 * ($relationship_hash{$qid}{$sid}->{'score'} / $relationship_hash{$qid}{$qid}->{'score'}))]++;
		$match_found = 1;
		last; #quit after finding best nonself match
	    }
	    if (!$match_found) {
		$self_hist{$qtag}[0]++;
	    }
	}
	foreach my $index (0..100) { # normalize self histogram
	    $self_hist{$qtag}[$index] /= $num_qids;
	}
	if ($histogramfile) {
	    print OUTHISTOGRAM "$qtag:$qtag";
	    foreach my $index (0..100) { # print self histogram
		print OUTHISTOGRAM "\t$self_hist{$qtag}[$index]";
	    }
	    print OUTHISTOGRAM "\n";
	}
    }

    foreach $qtag (@tag_array)  {  # start looping through genomes by order in tag file
	shift(@tmp_tag_array);
	foreach $stag (@tmp_tag_array) { # loop through all other genomes
	    my @qtag_best_hist;
	    my @qtag_second_hist;
	    my @stag_best_hist;
	    my @stag_second_hist;
	    foreach my $index (0..100) { # initialize histograms
		$qtag_best_hist[$index] = 0;
		$qtag_second_hist[$index] = 0;
		$stag_best_hist[$index] = 0;
		$stag_second_hist[$index] = 0;
	    }
	    my $num_qtag_matches = 0;
	    my $num_qids = 0;
	    foreach $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
		$num_qids++;
		if (!defined $Qbytaghash{$qid}) {
		    $qtag_second_hist[0]++;
		    next;# need to do this to prevent keys %{ $Qbytaghash{$qid}{$stag} } causing $Qbytaghash{$qid}{$stag} to be defined later
		}
		if (!defined $Qbytaghash{$qid}{$stag}) {
		    $qtag_second_hist[0]++;
		    next;
		}
		$num_qtag_matches++;
		my $match_found = 0;
		my $best_found = 0;
		foreach $sid (sort {$Qbytaghash{$qid}{$stag}{$b}->{'score'} <=> $Qbytaghash{$qid}{$stag}{$a}->{'score'}} keys %{ $Qbytaghash{$qid}{$stag} })  {
		    if (!$best_found) {
			$qtag_best_hist[int (100 * ($relationship_hash{$qid}{$sid}->{'score'} / $relationship_hash{$qid}{$qid}->{'score'}))]++;
			$best_found = 1;
			next;
		    }
		    $qtag_second_hist[int (100 * ($relationship_hash{$qid}{$sid}->{'score'} / $relationship_hash{$qid}{$qid}->{'score'}))]++;
		    $match_found = 1;
		    last; #quit after finding second best match
		}
		if (!$match_found) {
		    $qtag_second_hist[0]++;
		}
	    }
	    foreach my $index (0..100) { # normalize histograms
		$qtag_best_hist[$index] /= $num_qtag_matches;
		$qtag_second_hist[$index] /= $num_qids;
	    }

	    my $num_stag_matches = 0;
	    my $num_sids = 0;
	    foreach $sid ( keys %{ $Tagbyfeatnamehash{$stag} } )  { # go through featnames of query genome to look for matches in other genomes
		$num_sids++;
		if (!defined $Qbytaghash{$sid}) {
		    $stag_second_hist[0]++;
		    next;# need to do this to prevent keys %{ $Qbytaghash{$qid}{$stag} } causing $Qbytaghash{$qid}{$stag} to be defined later
		}
		if (!defined $Qbytaghash{$sid}{$qtag}) {
		    $stag_second_hist[0]++;
		    next;
		}
		$num_stag_matches++;
		my $match_found = 0;
		my $best_found = 0;
		foreach $qid (sort {$Qbytaghash{$sid}{$qtag}{$b}->{'score'} <=> $Qbytaghash{$sid}{$qtag}{$a}->{'score'}} keys %{ $Qbytaghash{$sid}{$qtag} })  {
		    if (!$best_found) {
			$stag_best_hist[int (100 * ($relationship_hash{$sid}{$qid}->{'score'} / $relationship_hash{$sid}{$sid}->{'score'}))]++;
			$best_found = 1;
			next;
		    }
		    $stag_second_hist[int (100 * ($relationship_hash{$sid}{$qid}->{'score'} / $relationship_hash{$sid}{$sid}->{'score'}))]++;
		    $match_found = 1;
		    last; #quit after finding second best match
		}
		if (!$match_found) {
		    $stag_second_hist[0]++;
		}
	    }
	    foreach my $index (0..100) { # normalize histograms
		$stag_best_hist[$index] /= $num_stag_matches;
		$stag_second_hist[$index] /= $num_sids;
	    }
	    if ($histogramfile) {
		print OUTHISTOGRAM "$qtag:$stag:best";
		foreach my $index (0..100) { # print qtag best histogram
		    print OUTHISTOGRAM "\t$qtag_best_hist[$index]";
		}
		print OUTHISTOGRAM "\n";
		print OUTHISTOGRAM "$qtag:$stag:second";
		foreach my $index (0..100) { # print qtag second histogram
		    print OUTHISTOGRAM "\t$qtag_second_hist[$index]";
		}
		print OUTHISTOGRAM "\n";
		print OUTHISTOGRAM "$stag:$qtag:best";
		foreach my $index (0..100) { # print stag best histogram
		    print OUTHISTOGRAM "\t$stag_best_hist[$index]";
		}
		print OUTHISTOGRAM "\n";
		print OUTHISTOGRAM "$stag:$qtag:second";
		foreach my $index (0..100) { # print stag second histogram
		    print OUTHISTOGRAM "\t$stag_second_hist[$index]";
		}
		print OUTHISTOGRAM "\n";
	    }

	    my $max_value = 0;
	    my $max_cutoff = 1.0;
	    my $cur_value = 0;
	    for (my $index = 100; $index > 0; $index--) {
		$cur_value += $qtag_best_hist[$index];
		$cur_value += $stag_best_hist[$index];
		$cur_value -= $qtag_second_hist[$index];
		$cur_value -= $stag_second_hist[$index];
		$cur_value -= $self_hist{$qtag}[$index];
		$cur_value -= $self_hist{$stag}[$index];
		if ($cur_value > $max_value) {
		    $max_value = $cur_value;
		    $max_cutoff = $index / 100.0;
		}
	    }
	    $paralog_cutoff{$qtag}{$stag} = $paralog_cutoff{$stag}{$qtag} = $max_cutoff;
	}
    }

    print STDERR Dumper(\%paralog_cutoff) if ($DEBUG);

    close (OUTHISTOGRAM) if ($histogramfile == 1);

			
}

sub calc_synbibest  { # subroutine to find bidrectional best synteny scores and populate Qbytaghash and relationship_hash with them

    my $qid = ""; # query feat_name
    my $sid = ""; # subject feat_name
    my $qtag = ""; # query genome tag
    my $stag = ""; # subject tag

    #mark best synteny score per genome per qid
    foreach $qtag (@tag_array)  {  # start looping through by order in tag file (reference is first)
	foreach $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
	    if (!defined $Qbytaghash{$qid}) {
		next;# need to do this to prevent keys %{ $Qbytaghash{$qid}{$stag} } causing $Qbytaghash{$qid}{$stag} to be defined later
	    }
	    foreach $stag (keys %{ $Qbytaghash{$qid} })  {
		if ($qtag eq $stag) {
		    next; #skip matches within the same genome
		}
		foreach $sid (sort {$Qbytaghash{$qid}{$stag}{$b}->{'synteny'} <=> $Qbytaghash{$qid}{$stag}{$a}->{'synteny'}} keys %{ $Qbytaghash{$qid}{$stag} })  {
		    $relationship_hash{$qid}{$sid}->{'synbest'} = 1;
		    last; #quit after marking best synteny score
		}
	    }
	}
    }

    #mark bidirectionally best synteny scores per genome per qid if it exists
    foreach $qid (keys %feat_hash)  { # go through all featnames
	foreach $sid (keys %{ $relationship_hash{$qid} } )  {
		if ($relationship_hash{$qid}{$sid}->{'synbest'} && $relationship_hash{$sid}{$qid}->{'synbest'}) {
		    $relationship_hash{$qid}{$sid}->{'synbibest'}  = 1;
		}
	}
    }
			
}

sub calc_cliques  { # subroutine to find bidrectional best cliques and populate Qbytaghash and relationship_hash with them

    foreach my $qid (keys %feat_hash)  { # go through all featnames
	my %clique_ids = ();
	my %clique_tags = ();
	my $clique_already_tried = 0;
	my $clique_all = 1;
	my $clique_failed = 0;
	my $num_clique_matches;
	my $qtag = $FeatnameLookupTag_hash{$qid};
	my $clique_top = 1;
	foreach my $sid (sort {$relationship_hash{$qid}{$b}->{'BSR'} <=> $relationship_hash{$qid}{$a}->{'BSR'}} keys %{ $relationship_hash{$qid} } )  {
	    print STDERR "Clique $qid $sid ($relationship_hash{$qid}{$sid}->{'BSR'}) :$relationship_hash{$qid}{$sid}->{'bibest'}\n" if ($DEBUG);
	    if ($qid eq $sid) {
		next; #skip matches to same protein
	    }
	    my $stag = $FeatnameLookupTag_hash{$sid};
	    if ($relationship_hash{$qid}{$sid}->{'bibest'}) {
		if (defined $clique_tags{$stag}){
		    print STDERR "WARNING!!! more than one bibest for qid $qid $qtag: stag $stag: $sid\n";
		}
		if (defined $relationship_hash{$qid}{$sid}->{'clique_all'}) {
		    $clique_already_tried = 1;
		    print STDERR "already tried :$relationship_hash{$qid}{$sid}->{'clique_all'}\n" if ($DEBUG);
		    last;
		}
		$clique_ids{$sid} = 1;
		$clique_tags{$stag} = 1;
	    } else {
		$clique_all = 0; #all matches are not part of the clique
		last; #end of possible clique
	    }
	}
	if ($clique_already_tried) {
	    next; #tried this clique already from some other starting point
	}
	$num_clique_matches = keys %clique_ids;
	if ($num_clique_matches < 2) {
	    $clique_failed = 1;
	    $clique_ids{$qid} = 1; #add the query into the clique
	} else {
	    my @clique_ids_to_check = keys %clique_ids;
	    $clique_ids{$qid} = 1; #add the query into the clique
	    foreach my $cqid (@clique_ids_to_check) {
		my $num_cqid_clique_matches = 0;
		my $clique_cqid_failed = 0;
		my @sort_array = sort {$relationship_hash{$cqid}{$b}->{'BSR'} <=> $relationship_hash{$cqid}{$a}->{'BSR'}} keys %{ $relationship_hash{$cqid} };
		my $num_cqid_matches = @sort_array;
		foreach my $csid (@sort_array)  {
		    if ($cqid eq $csid) {
			$num_cqid_matches--;
			next; #skip matches to same protein
		    }
		    if ($relationship_hash{$cqid}{$csid}->{'bibest'}) {
			if (!defined $clique_ids{$csid}){
			    $clique_cqid_failed = 1; #should not have an additional bibest not in clique
			    last;
			}
			$num_cqid_clique_matches++;
			if (defined $relationship_hash{$cqid}{$csid}->{'clique_all'}) {
			    print STDERR "WARNING!!! this clique already marked as tried for qid $qid $qtag: cqid $cqid: csid $csid\n" if ($DEBUG);
			    last;
			}
		    } else {
			last; #end of possible clique
		    }
		}
		if (($clique_cqid_failed) || ($num_cqid_clique_matches != $num_clique_matches)) {
		    $clique_failed = 1;
		    last;
		}
		if ($num_cqid_matches != $num_clique_matches) {
		    $clique_all = 0;
		}
	    }
	}
	if ($clique_failed) {
	    $clique_all = 0;
	    $clique_top = 0;
	} else {
	    $clique_top = $num_clique_matches + 1; # include the query in the clique count
	    if ($clique_all) {
		$clique_all = $num_clique_matches + 1; # include the query in the clique count
	    }
	}
	if ($DEBUG){
	    my @tmp_clique_ids = (keys %clique_ids);
	    print STDERR "$clique_all : $clique_top @tmp_clique_ids\n";
	}
	foreach my $cqid (keys %clique_ids) {
	    foreach my $csid (keys %clique_ids) {
		if ($cqid eq $csid) {
		    next; #skip matches to same protein
		}
		if (defined $relationship_hash{$cqid}{$csid}) {
		    $relationship_hash{$cqid}{$csid}->{'clique_all'} = $clique_all;
		    $relationship_hash{$cqid}{$csid}->{'clique_top'} = $clique_top;
		}
	    }
	}
    }
			
    #mark any unmarked clique_top and clique_all as 0
    foreach my $qid (keys %feat_hash)  { # go through all featnames
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  {
	    if (!defined $relationship_hash{$qid}{$sid}->{'clique_top'}) {
		$relationship_hash{$qid}{$sid}->{'clique_all'} = 0;
		$relationship_hash{$qid}{$sid}->{'clique_top'} = 0;
	    }
	}
    }
			
}

sub calc_synteny  { # subroutine to calculate synteny scores and populate Qbytaghash and relationship_hash with them

   #call synteny score for each relationship
    foreach my $qid (keys %feat_hash)  { # go through all featnames
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  {
	    if ($qid eq $sid) {
		next; #skip matches to same protein
	    }
	    my $qtag = $FeatnameLookupTag_hash{$qid};
	    my $stag = $FeatnameLookupTag_hash{$sid};
	    if ($qtag eq $stag) {
		next; #skip matches to same genome
	    }
	    print STDERR "$qid $sid $relationship_hash{$qid}{$sid}->{'id'} $relationship_hash{$qid}{$sid}->{'eval'} $relationship_hash{$qid}{$sid}->{'score'} $relationship_hash{$qid}{$sid}->{'best'} $relationship_hash{$qid}{$sid}->{'bibest'} $relationship_hash{$qid}{$sid}->{'clique_top'} $relationship_hash{$qid}{$sid}->{'clique_all'}\n" if($DEBUG);
	    $relationship_hash{$qid}{$sid}->{'synteny'} = &synteny($qtag, $qid, $stag, $sid);
	}
    }
			
}

sub merge_clusters {#&merge_clusters(qid, sid) is used to test if two clusters can be merged and if so to merge them
    my ($qid, $sid) = @_;
    my $merged = [];
    my $failed_merge = 0;
    
    #assumes the cluster arrays are sorted first by genome tag and then by protein id
    my $index = 0; #need to keep track of where we are in the subject array
    foreach my $cur_qid_tag ( @{ $clusters{$qid} } ) {
	while (($index <= $#{ $clusters{$sid} }) && (($clusters{$sid}->[$index])->{'tag'} < $cur_qid_tag->{'tag'})) {
	    push(@{ $merged }, $clusters{$sid}->[$index]);
	    print STDERR Dumper($merged) if ($DEBUG);
	    $index++;
	}
	if ($index <= $#{ $clusters{$sid} }) {
	    if (($clusters{$sid}->[$index])->{'tag'} == $cur_qid_tag->{'tag'}) {
		if (($clusters{$sid}->[$index])->{'id'} eq $cur_qid_tag->{'id'}) {
		    $index++; #protein is shared between clusters so don't put into merged twice
		} else {
		    $failed_merge = 1; #cannot have two different proteins from the same genome
		    last;
		}
	    }
	}
	push(@{ $merged }, $cur_qid_tag);
	print STDERR Dumper($merged) if ($DEBUG);
    }
    if (!$failed_merge) {
	while ($index <= $#{ $clusters{$sid} }) {
	    push(@{ $merged }, $clusters{$sid}->[$index]);
	    print STDERR Dumper($merged) if ($DEBUG);
	    $index++;
	}
	foreach my $cur_id_tag ( @{ $clusters{$qid} }, @{ $clusters{$sid} } ) {
	    $clusters{$cur_id_tag->{'id'}} = $merged;
	}
    } else {
	print STDERR "failed merge\n" if ($DEBUG);
    }
    
}

sub calc_clusters { # greedily compute clusters by starting with largest relationship score to merge existing clusters (starting with every protein as a singleton) and constringing a cluster to only contain one protein from each genome

    # start by creating a sorted array of relationships from the relationship_hash
    my @relationship_array = ();
    foreach my $qid (keys %feat_hash)  { # go through all featnames
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  { # go through all relationships
	    if ($qid eq $sid) {
		next; #skip matches to same protein
	    }
	    my $qtag = $FeatnameLookupTag_hash{$qid};
	    my $stag = $FeatnameLookupTag_hash{$sid};
	    if ($qtag eq $stag) {
		next; #skip matches to same genome
	    }
	    if (!$relationship_hash{$qid}{$sid}->{'synbibest'}) {
		next; # ignore matches which are not bidirectionally best synteny matches
	    }
	    my $score = $relationship_hash{$qid}{$sid}->{'synteny'};
	    if (defined $relationship_hash{$sid}{$qid}) {
		if ($qid gt $sid) {
		    next; # only enter qid,sid pairs once not symmetrically since they indicate the same cluster join
		}
		if ($relationship_hash{$sid}{$qid}->{'synteny'} > $score) {# use the larger synteny score
		    $score = $relationship_hash{$sid}{$qid}->{'synteny'};
		}
	    }
	    my $hash_ref = {};
	    $hash_ref->{'query'} = $qid;
	    $hash_ref->{'subject'} = $sid;
	    $hash_ref->{'score'} = $score;
	    push (@relationship_array, $hash_ref);
	}
    }
    @relationship_array = sort {$b->{'score'} <=> $a->{'score'}} ( @relationship_array );
    print STDERR Dumper(\@relationship_array) if ($DEBUG);

    foreach my $qid (keys %feat_hash)  { # go through all featnames and initialize clusters as singletons
	$clusters{$qid} = [{ 'id' => $qid, 'tag' => $LookupTagIndex{$FeatnameLookupTag_hash{$qid}} }];
    }
    print STDERR Dumper(\%clusters) if ($DEBUG);
    foreach my $match ( @relationship_array ) {
	if ($match->{'score'} < $min_synteny_threshold) {
	    last; # stop using matches when the scores get too low
	}
	if ($clusters{$match->{'query'}} == $clusters{$match->{'subject'}}) {
	    next; #already in the same cluster
	}
	print STDERR "merge ($match->{'query'}, $match->{'subject'}): $match->{'score'}\n" if ($DEBUG);
	&merge_clusters($match->{'query'}, $match->{'subject'});
    }
    print STDERR "Done clustering\n" if ($DEBUG);
    print STDERR Dumper(\%clusters) if ($DEBUG);
			
}

sub write_files  {

    my ($choice, $i, $column, $query_featname, $subject_tag, $subject_featname) = @_;
    my $ratio = "";
    my $x = "";

    if ($choice == "1")  {
	if ($i == 1)  {
	    print OUTVENN "$query_featname" if ($vennfile);
	    print MICRO "$query_featname\t$feat_hash{$query_featname}->{'header'}" if ($microarray);
	    print PAULSEN "$query_featname \[$feat_hash{$query_featname}->{'header'}\]" if ($hitsfile);
	    print READ "$query_featname\t$feat_hash{$query_featname}->{'header'}" if ($normalizefile);
	    print OUTVENNID "$query_featname" if ($vennfile);
	    print OUTID "$query_featname\t$feat_hash{$query_featname}->{'header'}" if ($vennfile);
	    
	}
	else  {
	    print OUTVENN "\t$query_featname" if ($vennfile);
	    print OUTVENNID "\t$query_featname" if ($vennfile);
	}
    } elsif ($choice == "2")  {
	my $score;
	my $percent_id;
	if (!defined $relationship_hash{$query_featname}{$subject_featname}->{'id'}) {
	    print STDERR "undefined relationship_hash for $query_featname and $subject_featname id\n" if ($DEBUG);
	    $score = 0;
	    $percent_id = 0;
	} else {
	    $score = $relationship_hash{$query_featname}{$subject_featname}->{'score'};
	    $percent_id = $relationship_hash{$query_featname}{$subject_featname}->{'id'};
	}
        if ($column != 1)  {
	    print OUTVENN "\t" if ($vennfile);
	    print OUTVENNID "\t" if ($vennfile);
	    if ($i == 1) { # only print the reference as query, not the leftovers from the other genomes
		printf MICRO "\t" if ($microarray);
		print PAULSEN "\t" if ($hitsfile);
		print READ "\t" if ($normalizefile);
		printf OUTID "\t" if ($vennfile);
	    }
	}
	print OUTVENN "$subject_featname" if ($vennfile);
	printf OUTVENNID "$subject_featname (%5.2f", $percent_id if ($vennfile);
	print OUTVENNID "%)" if ($vennfile);
	if ($i == 1) { # only print the reference as query, not the leftovers from the other genomes
            #$Qbytaghash{$qid}{$stag}->{$sid}->{'BSR'}
	    $x = $score/$relationship_hash{$query_featname}{$query_featname}->{'score'};  # determine normalized BLAST score
	    ## perhaps this can be simplified by getting data from Qbytaghash?
	    $ratio = -99*$x+100; # convert to scale 1 = perfect match, 100 = no match ( arthematic supplied by Emmanuel Mongodin )
	    printf MICRO "%2.2f", $ratio if ($microarray);
	    print PAULSEN "$subject_featname \[$feat_hash{$subject_featname}->{'header'}\]" if ($hitsfile);
	    print READ "$x" if ($normalizefile);
	    printf OUTID "%5.2f", $percent_id if ($vennfile);
	}

	delete $Tagbyfeatnamehash{$subject_tag}{$subject_featname};  # remove the orf from the matched genome, so not printed again
    } elsif ($choice == "3")  {
	if ($column != 1)  {
	    print OUTVENN "\t" if ($vennfile); 
	    print OUTVENNID "\t" if ($vennfile);
	    if ($i == 1) { # only print the reference as query, not the leftovers from the other genomes
                printf MICRO "\t" if ($microarray);
                print PAULSEN "\t" if ($hitsfile);
                print READ "\t" if ($normalizefile);
                print OUTID "\t" if ($vennfile);
	    }
	}
	print OUTVENN "----------" if ($vennfile);  # if no match print the default delimiter ----------
	if ($i == 1) { # only print the reference as query, not the leftovers from the other genomes
	    $x = 0;
	    $ratio = -99*$x+100;
	    printf MICRO "%2.2f", $ratio if ($microarray);
	    print READ "$x" if ($normalizefile);
	}
    } elsif ($choice == "4")  {
	print OUTVENN "\n" if ($vennfile);
	print OUTVENNID "\n" if ($vennfile);
	if ($i == 1)  {
	    print MICRO "\n" if ($microarray);
	    print PAULSEN "\n" if ($hitsfile);
	    print READ "\n" if ($normalizefile);
	    print OUTID "\n" if ($vennfile);
	}
    }
}

sub synteny {#&synteny($genome_tag, $query_featname, $subject_tag, $subject_featname)
# need to make genome_hash two dimensional to account for multiple genomes and multiple molecules
    my ($query_tag, $query_featname, $subject_tag, $subject_featname) = @_;
    my $QueryArrayIndex = $TagByPointer{$query_featname}; # points to the position in the genome_hash array
    my $SubjectArrayIndex = $TagByPointer{$subject_featname}; # points to the position in the genome_hash array
    my $query_orient = $feat_hash{$query_featname}->{'orient'};
    my $subject_orient = $feat_hash{$subject_featname}->{'orient'};
    my $orient = $query_orient * $subject_orient;
    my $i = "";
    my $j = "";
    my $OuterLoopEnd = "";
    my $OuterLoopBeg = "";
    my $InnerLoopEnd = "";
    my $InnerLoopBeg = "";
    my $query = "";
    my $subject = "";
    my $Qasmbl_id = $AssemblyLookup_hash{$query_featname};
    my $Sasmbl_id = $AssemblyLookup_hash{$subject_featname};    
    my $query_array_last = $#{ $genome_hash{$query_tag}{$Qasmbl_id} };
    my $subject_array_last = $#{ $genome_hash{$subject_tag}{$Sasmbl_id} };
    my $synteny_range = 5; #this is the number of proteins on either side of the query and subject to compute synteny
    my $num_synteny = 0;
    my $individual_score;
    my $max_individual_score;
    my $total_score = 0;
    my %max_subject_hash = ();
    my $max_subject;
    my $max_j;
    my $abs_rel_pos_query;
    my $abs_rel_pos_subject;
    my $dist_from_center;
    my $query_5pend; #+1 if on the 5' end side of the subject ORF -1 3'end
    my $subject_5pend; #+1 if on the 5' end side of the subject ORF -1 3'end


    print STDERR "query $query_featname ($QueryArrayIndex)     subject $subject_featname ($SubjectArrayIndex)\n" if ($DEBUG);

    $OuterLoopEnd = $QueryArrayIndex + $synteny_range;
    if ($OuterLoopEnd > $query_array_last) {# check to make sure we don't extend past the last array element 
	$OuterLoopEnd = $query_array_last; # if so, make OuterLoopEnd equal to the size of the array
    }
    $OuterLoopBeg = $QueryArrayIndex - $synteny_range;
    if ($OuterLoopBeg < 0)  { # don't let this run off the beginning of the array
	$OuterLoopBeg = 0; # adjust to the first position of the array
    }
    for ($i = $OuterLoopBeg; $i <= $OuterLoopEnd; $i++) { # iterate over synteny range of query genome
	print STDERR "       qpos = $i (${ $genome_hash{$query_tag}{$Qasmbl_id} }[$i])\n" if ($DEBUG);
	$abs_rel_pos_query = abs ($i - $QueryArrayIndex);
	$query = ${ $genome_hash{$query_tag}{$Qasmbl_id} }[$i];
	if (($QueryArrayIndex - $i) * $query_orient > 0) {
	    $query_5pend = 1;
	} else {
	    $query_5pend = -1;
	}
	$InnerLoopEnd = $SubjectArrayIndex + $synteny_range;
	if ($InnerLoopEnd > $subject_array_last) {# check to make sure we don't extend past the last array element 
	    $InnerLoopEnd = $subject_array_last; # if so, make InnerLoopEnd equal to the size of the array
	}
	$InnerLoopBeg = $SubjectArrayIndex - $synteny_range;
	if ($InnerLoopBeg < 0)  { # don't let this run off the beginning of the array
	    $InnerLoopBeg = 0; # adjust to the first position of the array
	}
	$max_individual_score = 0;
	for ($j = $InnerLoopBeg; $j <= $InnerLoopEnd; $j++) { # iterate over synteny range of subject genome
	    $abs_rel_pos_subject = abs ($j - $SubjectArrayIndex);
	    if ($abs_rel_pos_query > $abs_rel_pos_subject) {
		$dist_from_center = $abs_rel_pos_query;
	    } else {
		$dist_from_center = $abs_rel_pos_subject;
	    }
	    $subject = ${ $genome_hash{$subject_tag}{$Sasmbl_id} }[$j];
	    if (($SubjectArrayIndex - $j) * $subject_orient > 0) {
		$subject_5pend = 1;
	    } else {
		$subject_5pend = -1;
	    }
	    if ((defined $relationship_hash{$query}{$subject}) && (defined $relationship_hash{$subject}{$query}) && (defined $relationship_hash{$query}{$subject}->{'BSR'}) && (defined $relationship_hash{$subject}{$query}->{'BSR'})) {#need to check for relationship_hash existence first to avoid creating it when we check for BSR existence
		$individual_score = 1;
		if ($relationship_hash{$query}{$subject}->{'best'}) {
		    $individual_score += 2;
		    if ($relationship_hash{$query}{$subject}->{'bibest'}) {
			$individual_score += 5;
			if ($relationship_hash{$query}{$subject}->{'clique_top'}) {
			    $individual_score += (5 * $relationship_hash{$query}{$subject}->{'clique_top'}) / $genome_number;
			    if ($relationship_hash{$query}{$subject}->{'clique_all'}) {
				$individual_score += (5 * $relationship_hash{$query}{$subject}->{'clique_all'}) / $genome_number;
			    }
			}
		    }
		}
		if (($orient * ($feat_hash{$query}->{'orient'} * $feat_hash{$subject}->{'orient'})) > 0) {
		    $individual_score += 1;
		}
		if ($query_5pend * $subject_5pend > 0) {
		    $individual_score += 1;
		}
		$individual_score += (($synteny_range + 1) - abs ($abs_rel_pos_query - $abs_rel_pos_subject));
		$individual_score *= exp ((log 3) *(($synteny_range + 1) - $dist_from_center));
		if ($individual_score > $max_individual_score) {
		    $max_individual_score = $individual_score;
		    $max_subject = $subject;
		    $max_j = $j;
		}
		print STDERR "         spos = $j (${ $genome_hash{$subject_tag}{$Sasmbl_id} }[$j]) [$relationship_hash{$query}{$subject}->{'id'}] <=> [$relationship_hash{$subject}{$query}->{'id'}] $individual_score $relationship_hash{$query}{$subject}->{'best'} $relationship_hash{$query}{$subject}->{'bibest'} $relationship_hash{$query}{$subject}->{'clique_top'} $relationship_hash{$query}{$subject}->{'clique_all'} (orient: $orient $feat_hash{$query}->{'orient'} $feat_hash{$subject}->{'orient'}) (5p: $query_5pend $subject_5pend)\n" if ($DEBUG);
	    }
	}
	$total_score += $max_individual_score;
	if ($max_individual_score > 0) {
	    if (!defined $max_subject_hash{$max_subject}) {
		$max_subject_hash{$max_subject} = $max_individual_score;
		if (($i != $QueryArrayIndex) && ($max_j != $SubjectArrayIndex)) { # don't count synteny for the match being evaluated
		    $num_synteny++;
		}
	    } else {
		if ($max_subject_hash{$max_subject} >= $max_individual_score) {
		    $total_score -= $max_individual_score;
		} else {
		    $total_score -= $max_subject_hash{$max_subject};
		    $max_subject_hash{$max_subject} = $max_individual_score;
		}
	    }
	}
    }
    print STDERR "     num_synteny $num_synteny\n" if ($DEBUG);
    print STDERR "     total_score $total_score\n" if ($DEBUG);
    return($total_score);

}

#GGS new subroutine to detect frameshifts/protein fragments and record them
#need to consider only doing this based on bidirectionally best blast matches
#would be better to take orientation and order of matches against the query into account
sub frameshifts {

    my %max_frames = ();        #value is counter for best/longest fragment of a frameshifted set, Key = feat_name
    my %con_comp_frames = ();   #value is counter for sets of frameshifted fragments, Keys are feat_name, feat_name

    my $qid = ""; # query feat_name
    my $sid = ""; # subject feat_name
    my $qtag = ""; # query genome tag
    my $stag = ""; # subject genome tag

    foreach $qtag (@tag_array)  {  # start looping through genomes by order in tag file
	foreach $qid (keys %{ $Tagbyfeatnamehash{$qtag} } ) { # go through featnames of each genome to look for matches in other genomes
	    foreach $stag (keys %{ $Qbytaghash{$qid} } ) {# if query protein matches anything in subject genome, lets drill through each match
		print STDERR "frameshifts: $qtag $qid $stag\n" if ($DEBUG);

		my $sort_asmbl_id_genome_order = sub {#sort routine used in frameshifts to sort proteins from a given genome by asmbl_id and then by position in the assembly/contig
		    my $asmbl_id_a = $AssemblyLookup_hash{$a};
		    my $asmbl_id_b = $AssemblyLookup_hash{$b};
		    if ($asmbl_id_a < $asmbl_id_b){ #sort ascending asmbl_id
			return -1;
		    } elsif ($asmbl_id_a > $asmbl_id_b){
			return 1;
		    } else { #sort ascending genome order
			$TagByPointer{$a} <=> $TagByPointer{$b};
		    }
		};

		if ($DEBUG) {
		    my @tmp_sids = keys %{ $Qbytaghash{$qid}{$stag} };
		    print STDERR "@tmp_sids\n";
		}

		my @qid_array_of_sids = sort $sort_asmbl_id_genome_order keys %{ $Qbytaghash{$qid}{$stag} };
		my $cur_index = 0;
		my $skip_index = 0;

		foreach my $feat_id (@qid_array_of_sids) {#remove sids which do not have full blast matches - they were created to preserve syntenny?
		    if (($relationship_hash{$qid}{$feat_id}->{"id"} < $fspercentid) || (!defined $relationship_hash{$qid}{$feat_id}->{"min_query"}) || (($relationship_hash{$qid}{$feat_id}->{"min_query"} == 0) && ($relationship_hash{$qid}{$feat_id}->{"max_query"} == 0))) { # remove low %ID matches here as well
			$skip_index++;
		    } else {
			$qid_array_of_sids[$cur_index] = $qid_array_of_sids[$skip_index];
			$cur_index++;
			$skip_index++;
		    }
		}
		$cur_index--;
		$#qid_array_of_sids = $cur_index;
		$cur_index = 0;

		if (@qid_array_of_sids < 2) {
		    next;
		}# only need to look for frameshifts with two or more matches
		my $prev_asmbl_id = -1;
		my $cur_asmbl_id;
		my $query_hits_span;
		my $min_query;
		my $max_query;
		my $cur_min_query;
		my $cur_max_query;
		my $cur_query_hit_len;
		my $sum_query_hit_len;
		my $start_index = 0;
		my $stop_index = 0;
		my $max_index;
		my $max_query_hit_len;
		my $failed_to_extend;
		my $prev_pos;
		my $cur_pos;

		print STDERR ">= 2 $#qid_array_of_sids\n" if ($DEBUG);

		if (@qid_array_of_sids == 2) {#look for fragments on the ends of two contigs (could be the same contig for circular contigs)
		    my $sid0 = $qid_array_of_sids[0];
		    my $sid1 = $qid_array_of_sids[1];
		    my $asmbl_id0 = $AssemblyLookup_hash{$sid0};
		    my $asmbl_id1 = $AssemblyLookup_hash{$sid1};
		    my $pos0 = $TagByPointer{$sid0};
		    my $pos1 = $TagByPointer{$sid1};
		    my $last_pos0 = $#{ $genome_hash{$stag}{$asmbl_id0} };
		    my $last_pos1 = $#{ $genome_hash{$stag}{$asmbl_id1} };
		    if ((($pos0 == 0) || ($pos0 == $last_pos0)) && (($pos1 == 0) || ($pos1 == $last_pos1))) {
			#matches are on ends of contigs so treat as possible fragments
			$min_query = $relationship_hash{$qid}{$sid0}->{"min_query"};
			$max_query = $relationship_hash{$qid}{$sid0}->{"max_query"};
			$query_hits_span = ($max_query - $min_query) + 1;
			$sum_query_hit_len = $query_hits_span;
			$max_query_hit_len = $query_hits_span;
			$max_index = 0;
			
			$cur_min_query = $relationship_hash{$qid}{$sid1}->{"min_query"};
			if ($cur_min_query < $min_query) {
			    $min_query = $cur_min_query;
			}
			$cur_max_query = $relationship_hash{$qid}{$sid1}->{"max_query"};
			if ($cur_max_query > $max_query) {
			    $max_query = $cur_max_query;
			}
			$query_hits_span = ($max_query - $min_query) + 1;
			$cur_query_hit_len = ($cur_max_query - $cur_min_query) + 1;
			$sum_query_hit_len += $cur_query_hit_len;
			if ($sum_query_hit_len < ($frameshift_length_ratio * $query_hits_span)){
			    if ($cur_query_hit_len > $max_query_hit_len){
				$max_query_hit_len = $cur_query_hit_len;
				$max_index = 1;
			    }
			    #record frameshift/fragments
			    if (!defined $con_comp_frames{$sid0}{$sid1}){
				$con_comp_frames{$sid0}{$sid1} = 1;
				$con_comp_frames{$sid1}{$sid0} = 1;
			    } else {
				$con_comp_frames{$sid0}{$sid1}++;
				$con_comp_frames{$sid1}{$sid0}++;
			    }
			    if (!defined $max_frames{$qid_array_of_sids[$max_index]}){
				$max_frames{$qid_array_of_sids[$max_index]} = 1;
			    } else {
				$max_frames{$qid_array_of_sids[$max_index]}++;
			    }
			}
			next; #already tested for fragments/frameshifts if on the ends of contigs
		    }
		}

		foreach $sid (@qid_array_of_sids) {
		    $cur_asmbl_id = $AssemblyLookup_hash{$sid};
		    $cur_pos = $TagByPointer{$sid};

		    print STDERR "$cur_asmbl_id $sid $cur_pos\n" if ($DEBUG);

		    if ($cur_asmbl_id == $prev_asmbl_id){
			$cur_min_query = $relationship_hash{$qid}{$sid}->{"min_query"};
			if ($cur_min_query < $min_query) {
			    $min_query = $cur_min_query;
			}
			$cur_max_query = $relationship_hash{$qid}{$sid}->{"max_query"};
			if ($cur_max_query > $max_query) {
			    $max_query = $cur_max_query;
			}
			if (($cur_min_query <= $max_missing_aa + 1) && ($cur_max_query + $max_missing_aa >= $feat_hash{$qid}->{'length'})) {
			    $failed_to_extend = 1; #not a protein fragment if almost a full length match
			} elsif ((($prev_pos + 1) != $cur_pos) && (($prev_pos - 1) != $cur_pos)) {
			    $failed_to_extend = 1; #not protein fragments/frameshift if not adjacent on the contig
			} else {
			    $query_hits_span = ($max_query - $min_query) + 1;
			    $cur_query_hit_len = ($cur_max_query - $cur_min_query) + 1;
			    $sum_query_hit_len += $cur_query_hit_len;
			    if ($sum_query_hit_len < ($frameshift_length_ratio * $query_hits_span)){
				$failed_to_extend = 0;
				if ($cur_query_hit_len > $max_query_hit_len){
				    $max_query_hit_len = $cur_query_hit_len;
				    $max_index = $cur_index;
				}
				$stop_index = $cur_index;
			    } else {
				$failed_to_extend = 1;
			    }
			}
			print STDERR "failed_to_exend = $failed_to_extend\n" if ($DEBUG);
		    }
		    


		    if (($cur_asmbl_id != $prev_asmbl_id) || $failed_to_extend || ($cur_index == $#qid_array_of_sids)){
			if ($stop_index > $start_index) {#record the frameshifted featnames
			    for (; $start_index < $stop_index; $start_index++){
				if (!defined $con_comp_frames{$qid_array_of_sids[$start_index]}{$qid_array_of_sids[$start_index + 1]}){
				    $con_comp_frames{$qid_array_of_sids[$start_index]}{$qid_array_of_sids[$start_index + 1]} = 1;
				    $con_comp_frames{$qid_array_of_sids[$start_index + 1]}{$qid_array_of_sids[$start_index]} = 1;
				} else {
				    $con_comp_frames{$qid_array_of_sids[$start_index]}{$qid_array_of_sids[$start_index + 1]}++;
				    $con_comp_frames{$qid_array_of_sids[$start_index + 1]}{$qid_array_of_sids[$start_index]}++;
				}
			    }
			    if (!defined $max_frames{$qid_array_of_sids[$max_index]}){
				$max_frames{$qid_array_of_sids[$max_index]} = 1;
			    } else {
				$max_frames{$qid_array_of_sids[$max_index]}++;
			    }
			}
			if ($cur_index == $#qid_array_of_sids){
			    last;
			}
			$min_query = $relationship_hash{$qid}{$sid}->{"min_query"};
			$max_query = $relationship_hash{$qid}{$sid}->{"max_query"};
			if (($min_query <= $max_missing_aa + 1) && ($max_query + $max_missing_aa >= $feat_hash{$qid}->{'length'})) {
			    $cur_asmbl_id = -1; #not a protein fragment if almost a full length match
			} else {
			    $query_hits_span = ($max_query - $min_query) + 1;
			    $sum_query_hit_len = $query_hits_span;
			    $max_query_hit_len = $query_hits_span;
			    $max_index = $cur_index;
			    $start_index = $cur_index;
			    $stop_index = $cur_index;
			}
		    }
		    $prev_asmbl_id = $cur_asmbl_id;
		    $prev_pos = $cur_pos;
		    $cur_index++;
		}

	    }
	}
    }

    open (OUTFRAMESHIFT, ">$outprefix" . "_frameshifts.txt");
    my $prev_asmbl_id = -1;
    my $cur_asmbl_id;
    my $prev_tag = "";
    my $cur_tag;
    my $max_feat_id;
    my $sort_fragments = sub {#sort routine used in frameshifts to sort proteins genome tag, then asmbl_id and then by best/longest fragment
	my $taga = $FeatnameLookupTag_hash{$a};
	my $tagb = $FeatnameLookupTag_hash{$b};
	my $asmbl_id_a = $AssemblyLookup_hash{$a};
	my $asmbl_id_b = $AssemblyLookup_hash{$b};
	if ($taga lt $tagb){ #sort ascending genome tag
	    return -1;
	} elsif ($taga gt $tagb){
	    return 1;
	} elsif ($asmbl_id_a < $asmbl_id_b){ #sort ascending asmbl_id
	    return -1;
	} elsif ($asmbl_id_a > $asmbl_id_b){
	    return 1;
	} else { #sort descending max fragment
	    return ($max_frames{$b} <=> $max_frames{$a});
	}
    };

    foreach $max_feat_id (sort $sort_fragments keys %max_frames) { #sort protein fragments by genome tag, then asmbl_id and then by best/longest fragment

	print STDERR "max_feat_id = $max_feat_id :full $feat_hash{$max_feat_id}->{'full'}\n" if ($DEBUG);

	if (!defined $max_frames{$max_feat_id}) {
	    next; #already included this in another set of fragments/frameshifts
	}
	$cur_tag = $FeatnameLookupTag_hash{$max_feat_id};
	$cur_asmbl_id = $AssemblyLookup_hash{$max_feat_id};
	my @linked_array = sort { $con_comp_frames{$max_feat_id}{$b} <=> $con_comp_frames{$max_feat_id}{$a} } keys %{ $con_comp_frames{$max_feat_id} };
	if ((@linked_array < 1) || (@linked_array > 2)) {
	    die ("Unexpected number of links for a protein fragment $max_feat_id - max index is $#linked_array\n");
	}
	if (($con_comp_frames{$max_feat_id}{$linked_array[0]} < $frames_link_thresh) || ($con_comp_frames{$max_feat_id}{$linked_array[0]} < $feat_hash{$max_feat_id}->{'full'}) || ($con_comp_frames{$max_feat_id}{$linked_array[0]} < $feat_hash{$linked_array[0]}->{'full'})) {
	    next; #max link less than threshold
	}
	print STDERR "frameshift weight: $con_comp_frames{$max_feat_id}{$linked_array[0]}\n" if ($DEBUG);
	print STDERR "linked_array[0] = $linked_array[0] :full $feat_hash{$linked_array[0]}->{'full'}\n" if ($DEBUG);
	delete $con_comp_frames{$max_feat_id}{$linked_array[0]};
	delete $con_comp_frames{$linked_array[0]}{$max_feat_id};
	if (@linked_array == 2) {
	    if (($con_comp_frames{$max_feat_id}{$linked_array[1]} < $frames_link_thresh) || ($con_comp_frames{$max_feat_id}{$linked_array[1]} < $feat_hash{$max_feat_id}->{'full'}) || ($con_comp_frames{$max_feat_id}{$linked_array[1]} < $feat_hash{$linked_array[1]}->{'full'})) {
		pop @linked_array;
	    } else {
		print STDERR "frameshift weight: $con_comp_frames{$max_feat_id}{$linked_array[1]}\n" if ($DEBUG);
		print STDERR "linked_array[1] = $linked_array[1] :full $feat_hash{$linked_array[1]}->{'full'}\n" if ($DEBUG);
		delete $con_comp_frames{$max_feat_id}{$linked_array[1]};
		delete $con_comp_frames{$linked_array[1]}{$max_feat_id};
	    }
	}
	if ($cur_tag ne $prev_tag){
	    print OUTFRAMESHIFT ">genome $cur_tag\n>asmbl_id $cur_asmbl_id\n";
	} elsif ($cur_asmbl_id != $prev_asmbl_id) {
	    print OUTFRAMESHIFT ">asmbl_id $cur_asmbl_id\n";
	}
	print OUTFRAMESHIFT $max_feat_id;
	my $linked_feat_id;
	while (defined ($linked_feat_id = shift(@linked_array))) {
	    print OUTFRAMESHIFT "\t$linked_feat_id";
	    delete $max_frames{$linked_feat_id};
	    #need to delete all references to this $linked_feat_id since it is a protein fragment we don't want to cluster
	    #need to rebuild %genome_hash and %TagByPointer after deleting all of the protein fragments
	    #do we want to update orf_counter to subtract out the frameshifts?
	    #$orf_counter{$FeatnameLookupTag_hash{$linked_feat_id}}->{'raw'}--;  # decrement the total orf counter when removing a fragment
	    #$orf_counter{$FeatnameLookupTag_hash{$linked_feat_id}}->{'used'}--;  # decrement the used orf counter when removing a fragment
	    #not clear what to do with %Paralogs - not touching for now

	    #delete matches between the deprecated protein fragments and other proteins
	    foreach my $relate_id (keys %{ $relationship_hash{$linked_feat_id} }) {#remove symmetric relationships
		delete $relationship_hash{$relate_id}{$linked_feat_id};
	    }
	    delete $relationship_hash{$linked_feat_id};

	    my $qtag = $FeatnameLookupTag_hash{$linked_feat_id};
	    foreach my $stag (keys %{ $Qbytaghash{$linked_feat_id} } ) { #delete references symmetrically
		    foreach my $sid (keys %{ $Qbytaghash{$linked_feat_id}{$stag} })  {
			delete $Qbytaghash{$sid}{$qtag}{$linked_feat_id};
		    }
	    }
	    delete $Qbytaghash{$linked_feat_id};

	    if ($output_fragments) {#do not delete the deprecated protein fragment so it will be output as a singleton cluster
		$feat_hash{$linked_feat_id}->{"fragment"} = 1;
	    } else {
		delete $feat_hash{$linked_feat_id};
		delete $Tagbyfeatnamehash{$qtag}{$linked_feat_id};
		delete $FeatnameLookupTag_hash{$linked_feat_id};
		delete $AssemblyLookup_hash{$linked_feat_id};
	    }

	    my @new_linked_array = keys %{ $con_comp_frames{$linked_feat_id} };
	    if (@new_linked_array > 1) {
		die ("Unexpected number of links for a protein fragment $linked_feat_id\n");
	    }
	    if (@new_linked_array < 1) {
		next;
	    }
	    my $new_linked_feat_id = shift(@new_linked_array);
	    if (($con_comp_frames{$linked_feat_id}{$new_linked_feat_id} < $frames_link_thresh) || ($con_comp_frames{$linked_feat_id}{$new_linked_feat_id} < $feat_hash{$linked_feat_id}->{'full'}) || ($con_comp_frames{$linked_feat_id}{$new_linked_feat_id} < $feat_hash{$new_linked_feat_id}->{'full'})) {
		next;
	    }
	    push (@linked_array, $new_linked_feat_id);
	    print STDERR "frameshift weight: $con_comp_frames{$linked_feat_id}{$new_linked_feat_id}\n" if ($DEBUG);
	    print STDERR "new_linked_feat_id = $new_linked_feat_id :full $feat_hash{$new_linked_feat_id}->{'full'}\n" if ($DEBUG);
	    delete $con_comp_frames{$linked_feat_id}{$new_linked_feat_id};
	    delete $con_comp_frames{$new_linked_feat_id}{$linked_feat_id};
	}
	print OUTFRAMESHIFT "\n";
	    
	$prev_tag = $cur_tag;
	$prev_asmbl_id = $cur_asmbl_id;
    }
    close (OUTFRAMESHIFT);


    #rebuild %genome_hash and %TagByPointer after deleting all of the protein fragments
    %TagByPointer = ();
    foreach my $tag (keys %genome_hash) {
	foreach my $asmbl_id (keys %{ $genome_hash{$tag} }) {
	    my $cur_index = 0;
	    my $skip_index = 0;
	    foreach my $feat_id (@{ $genome_hash{$tag}{$asmbl_id} }) {
		if ((!defined $feat_hash{$feat_id}) || ($feat_hash{$feat_id}->{"fragment"})) {
		    $skip_index++;
		} else {
		    $genome_hash{$tag}{$asmbl_id}->[$cur_index] = $genome_hash{$tag}{$asmbl_id}->[$skip_index];
		    $TagByPointer{$feat_id} = $cur_index;
		    $cur_index++;
		    $skip_index++;
		}
	    }
	    $cur_index--;
	    $#{ $genome_hash{$tag}{$asmbl_id} } = $cur_index;
	}
    }

}
			    
sub make_tables {

############ Generate a match table of bidirectional best hits ###############
    my $i = 0;
    
    open (OUTVENN, ">$outprefix" . "_matchtable.txt") if ($vennfile);
    open (OUTVENNID, ">$outprefix" . "_matchtable_id.txt") if ($vennfile);
    open (OUTID, ">$outprefix" . "_id.txt") if ($vennfile);
    open (MICRO, ">$outprefix" . "_micro.txt") if ($microarray == 1);
    open (PAULSEN, ">$outprefix" . "_hits.txt") if ($hitsfile == 1);
    open (READ, ">$outprefix" . "_BSR.txt") if ($normalizefile == 1);
    foreach my $genome_tag (@tag_array)  {  # start looping through by order in tag file (reference is first)
	my $genome_tag_index = $LookupTagIndex{$genome_tag};
	$i++; # increment the times through counter
	for my $query_featname ( sort keys %{ $Tagbyfeatnamehash{$genome_tag} } )  { # go through featnames of reference first to look for matches in other genomes
	    my $column = 0; # reset column counter with each query orf
	    my $subject_tag_index;
	    my $subject_featname = "";
	    for ($subject_tag_index = 0; $subject_tag_index < $genome_tag_index; $subject_tag_index++) {
		&write_files("3", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname);
	    }
	    # We are at $genome_tag now so we are looking at self comparisons, generate output and move on
	    &write_files("1", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname);#$subject_featname has not been defined yet here!!
	    $subject_tag_index++; #go past $genome_tag
	    # Now iterate through sorted (by subject_tag_index) cluster array to print out members of cluster
	    foreach my $cluster_member ( @{ $clusters{$query_featname} } ) {
		my $cluster_member_tag_index = $cluster_member->{'tag'};
		if ($cluster_member_tag_index < $subject_tag_index) {
		    if ($cluster_member->{'id'} ne $query_featname) {
			print STDERR "WARNING!!! cluster for $query_featname should have been output sooner!\n";
		    }
		    next; #skip $genome_tag member we already output for
		}
		while ($cluster_member_tag_index > $subject_tag_index) {
		    #output blanks for this genome
		    &write_files("3", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname);
		    $subject_tag_index++; #go to next genome
		}
		$subject_featname = $cluster_member->{'id'};
		&write_files("2", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname);
		$subject_tag_index++; #go to next genome
	    }
	    while ($subject_tag_index <= $#tag_array) {
		&write_files("3", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname);
		$subject_tag_index++; #go to next genome
	    }
	    &write_files("4", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname); #add newlines to files
	}
    }
    close (OUTVENN) if ($vennfile);
    close (MICRO) if ($microarray == 1);
    close (PAULSEN) if ($hitsfile == 1);
    close (READ) if ($normalizefile == 1);
    close (OUTVENNID) if ($vennfile);
    close (OUTID) if ($vennfile);
}

sub gather_dups  {

my $genome_tag = "";
my $orf = "";
my $dup = "";
my @temp = ();

open (OUTDUP, ">$outprefix" . "_paralogs.txt");
foreach $genome_tag (@tag_array)  {  # start looping through by order in tag file (reference is first)
  for $orf ( sort keys %{ $Paralogs{$genome_tag} } )  { # go through featnames
    print OUTDUP "$genome_tag: $orf";
    @temp = keys %{ $Paralogs{$genome_tag}{$orf} };
    if ($#temp > 1) {  #GGS should this be @temp instead of $#temp?
      for $dup ( sort keys %{ $Paralogs{$genome_tag}{$orf} } )  {
        print OUTDUP "\t$dup";
      }
    }
    else  {
      print OUTDUP "\t$temp[0]";
    }
    print OUTDUP "\n";
  }
  print OUTDUP "\n";
}
close (OUTDUP);
}

sub print_report  {

    my $key = "";
    open (OUT, ">$outprefix" . "_report.txt");
    print OUT " e-value cut-off: <= $evalue\n";
    print OUT "percent identity: >= $percentid\n";
    print OUT " length of match: >= $min_hit_length\n";
    print OUT "input .btab file: $btabpath/$btabfile\n\n";
    foreach $key (@tag_array) {
	print OUT "Raw $key = $orf_counter{$key}->{'raw'}\n";
	print OUT "Used $key = $orf_counter{$key}->{'used'}\n";
    }
    close (OUT);
}

sub option_help {

   system("clear");
   print STDERR <<_EOB_;
$prog  - Pan-genome Ortholog Clustering Tool, a heuristic computer program for pan-genomic analysis of closely related prokaryotic species or strains.
 Copy (C) 2011-2012  The J. Craig Venter Institute (JCVI).  All rights reserved

 License:  This program is free software: you can redistribute it and/or modify
           it under the terms of the GNU General Public License as published by
           the Free Software Foundation, either version 3 of the License, or
           (at your option) any later version.

           This program is distributed in the hope that it will be useful,
           but WITHOUT ANY WARRANTY; without even the implied warranty of
           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
           GNU General Public License for more details.

           You should have received a copy of the GNU General Public License
           along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Citation: Derrick E. Fouts, Lauren Brinkac, Erin Beck, Jason Inman, and Granger Sutton (2011) "PanOCT: Automated Clustering of Orthologs using 
           Conserved Gene Neighborhood for Pan-Genomic Analysis of Bacterial Strains and Closely Related Species" Nucleic Acids Res submitted.

Usage: $prog <options>
Example: PanOCT.pl -t example_blast.txt -f example_tags.txt -g example.gene_att -P example.pep -i 20 -F 1.33 -N Y -M Y -H Y -G Y
Version: $version
 Switch: -h for help\
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
 Output: All stored within a subdirectory of the current working directory (PWD)
          1) panoct_report.txt:  a file containing runtime parameters used (e-value, %id, match length, and blast file used)
          2) panoct_matchtable.txt:  a tab-delimited file containing PanOCT clusters, one cluster per line.  The first column is the reference genome and all subsequent columns are the remaining genomes
                                     in the order specified in the genome identifier "tags" file (specified with option f).
                                     e.g. NT08AB0001	NT16AB0001	NT17AB0001	NT20ABA0020
          3) panoct_matchtable_id.txt:  a tab-delimited file similar to panoct_matchtable.txt, but also containing the percent identity of each target protein in parentheses.
                                        e.g. NT08AB0001      NT16AB0001 (99.78%)     NT17AB0001 (99.26%)     NT20ABA0020 (99.78%) 
          4) panoct_id.txt:  a tab-delimited file containing reference protein annotation and percent identities to orthologs in each target genome.
                             e.g. NT08AB0001	chromosomal replication initiator protein DnaA	99.78	99.26	99.78
          5) panoct_frameshifts.txt:  a tab-delimited file containing proteins that are likely split due to frame-shifts.  It is organized by genome and assembly/contig 
                                      e.g. >genome ntab08
                                           >asmbl_id 1
                                           NT08AB3019	NT08AB3018	NT08AB3020

 Authors: Derrick E. Fouts, Ph.D. and Granger Sutton, Ph.D.
 Date: December 21, 2004; last revised November 02, 2011
_EOB_
    exit;
}

########################################  M A I N  #####################################################
print STDERR "fetching tags to search\n";
&get_tags;
print STDERR "Gathering protein sequence information from $pep_file\n";
&get_protein_info;
print STDERR "gathering gene attributes\n";
&get_gene_att;
print STDERR "getting data from .btab file\n";
&select_data_from_btab;
print STDERR "Calculating the Blast Score Ratio (BSR) ...\n";
&calc_BSR;
if ($frameshiftfile) {
    print STDERR "Determining frameshifts / protein fragments ...\n";
    &frameshifts;
}
print STDERR "Calculating reciprocal best hits (RBHs) ...\n";
&calc_bibest;
print STDERR "Calculating histograms ...\n";
&calc_score_histograms;
print STDERR "Calculating cliques of RBHs ...\n";
&calc_cliques;
print STDERR "Calculating synteny ...\n";
&calc_synteny;
print STDERR "Calculating clusters ...\n";
&calc_synbibest;
print STDERR "Calculating conserved gene neighborhood (CGN) scores ...\n";
&calc_clusters;
print STDERR "building table for microarray software\n" if ($microarray);
print STDERR "building table of hits\n" if ($hitsfile);
print STDERR "generating match table!\n" if ($vennfile);
&make_tables;
if ($dupsfile)  {
  print STDERR "Finding paralogs...\n";
  &gather_dups;
}
&print_report;
exit(0);
