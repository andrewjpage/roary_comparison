#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

### programs from BLAST
my $formatdb="/usr/bin/formatdb";
my $blastall="/usr/bin/blastall";

### programs from mcl
my $mcl="/software/pathogen/external/apps/usr/bin/mcl";

### programs from mafft
my $mafft="/software/pubseq/bin/mafft";

### programs from PHYLIP
my $seqboot="/software/pathogen/external/apps/usr/bin/seqboot";
my $neighbor="/software/pathogen/external/apps/usr/bin/neighbor";
my $consense="/software/pathogen/external/apps/usr/bin/consense";
my $dnaml="/software/pathogen/external/apps/usr/bin/dnaml";
my $dnadist="/software/pathogen/external/apps/usr/bin/dnadist";
my $dnapars="/software/pathogen/external/apps/usr/bin/dnapars";

my $count_tree=0;

my $sampleSize=8000; # when calculate the pan-genome size, we will sample $sampleSize combinations
                     # if the total combination number is larger than $sampleSize for specific genomes
                     # Surely, the number of $sampleSize is, the larger, the better. 
                     # However, the larger the $sampleSize is, the more time would be consumed.
                     # we suggest the range: 5000 ~ 20,000

#####################################################################
#  DOn't modify the following code, unless you know their functions
#####################################################################

my %opt=qw();
GetOptions(\%opt,"strains:s","input:s","output:s","cluster!","pangenome!","variation!","evolution!","function!","method:s","thread:i","score:f","evalue:f","coverage:f","local:f","global:f","identity:f","bootstrap:i","help|h!");

my @usage=qq(
====== Pan-Genome Analysis Pipeline (PGAP) ======
                   Version 1.1

Usage:   perl PGAP.pl [Options]

Options: 
  --strains    String    Input strains nicknames, and join them with '+', for example: A+B+C
  --input      String    Input data directory 
  --output     String    Result output directory
 
  --cluster              Run homologous gene clustering
  --pangenome            Run pan-genome analysis
  --variation            Run homologous clusters variation analysis
  --evolution            Run evolution analysis
  --function             Run Function analysis

  --method     String    GF for GeneFamily method,  and MP for MultiParanoid method
                           for GF: fast, but not very accurate
                               evalue, score, indentity, coverage are employed
                           for MP: slow, but more accurate
                               score, coverage, local, global are employed
  --thread     Int       Number of processors to use in blastall. [default:1]
  --score      Int       Minimum score in blastall. [default:40]
  --evalue     Decimal   Maximal E-value in blastall. [default:1e-10]
  --coverage   Decimal   Minimum alignment coverage for two homologous proteins. [default:0.5]
  --local      Decimal   Minimum local alignment overlap in MP method. [default:0.25]
  --global     Decimal   Minimum global alignment overlap in MP method. [default:0.5]
  --identity   Decimal   Minimum alignment indentity for two homologous proteins. [default:0.5]
  --bootstrap  Int       Bootstrap times for phylogenetics tree. [default:1]

  --h or help            Display this message
);


############# specified variable #############
my $inputDIR;
my $outputDIR;
my $run_cluster;
my $run_pangenome;
my $run_variation;
my $run_evolution;
my $run_function;
my $method="";
my $thread;
my $score;
my $identity;
my $evalue;
my $coverage;
my $global;
my $local;
my $bootstrap;


my %pep;
my %nuc;
my $spnum;
my @clusters;
my $Cluster;
my @SpecieCombination;
my @spID;
my %genenum;
my %aaAln;
my %ntAln;
my %cog;
my %description;
#my %aa4tree;  ### AA sequence for Phylogenetic Tree
my %nt4tree;  ### nucleotide sequence for Phylogenetic Tree
my @SNPPosition; ### SNP position
my $dieMessage="You did not run PGAP.pl in the program directory\n";
my $section;

######### common temporary variable #############
my $i;
my $j;
my $line;
my %tmpHash;
my @tmp;
my $tmp;
my $key;
my @row;
my $inparacount;
my $ClusterID;
my $orth;
my @content;
my $clusterName;
my @xdata;
my @ydata;
my @fit;

my $fit_A;
my $fit_A_interval;
my $fit_B;
my $fit_C;
my $fit_C_interval;
my $fit_Rsquare;



#### check option

my $opt_error=0;

if ((scalar(keys %opt) ==0) or (exists($opt{"help"}))) 
	{
		print join("\n",@usage)."\n";
		exit;
	}


###################### public info
### strains name
my @species;
if (exists($opt{"strains"})) 
	{
		@species=split(/\+/,$opt{"strains"});
		$spnum=scalar(@species);
	}else
		{
			print "Please assign strains nick name!\n";
			exit;
		}

### input data directory

if (exists($opt{"input"})) 
	{
		$inputDIR=$opt{"input"};
		if ($inputDIR!~/\/$/) 
			{
				$inputDIR=$inputDIR."/";
			}
	}else
		{
			print "Please assign input data directory!\n\n";
			exit;
		}
### output data directory

if (exists($opt{"output"})) 
	{
		$outputDIR=$opt{"output"};
		if ($outputDIR!~/\/$/) 
			{
				$outputDIR=$outputDIR."/";
			}
	}else
		{
			print "Please assign result output directory!\n\n";
			exit;
		}

###################### section info

if (exists($opt{"cluster"})) 
	{
		$run_cluster=1;
	}else
		{
			$run_cluster=0;
		}

if (exists($opt{"pangenome"})) 
	{
		$run_pangenome=1;
	}else
		{
			$run_pangenome=0;
		}

if (exists($opt{"variation"})) 
	{
		$run_variation=1;
	}else
		{
			$run_variation=0;
		}

if (exists($opt{"evolution"})) 
	{
		$run_evolution=1;
	}else
		{
			$run_evolution=0;
		}

if (exists($opt{"function"})) 
	{
		$run_function=1;
	}else
		{
			$run_function=0;
		}

if ($run_cluster) 
	{
		### method
		if (exists($opt{"method"})) 
			{
				$method=uc($opt{"method"});
				if ($method!~/^GF$/ and $method!~/^MP$/) 
					{
						print "Unknown method: ".$opt{"method"}."\n";
						exit;
					}
			}else
				{
					print "Please assign the cluster method!\n\n";
					exit;
				}

		##thread
				if (exists($opt{"thread"})) 
					{
						$thread=$opt{"thread"};
						if ($thread==0) 
							{
								print "please assign an applicable thread value.\n";
								exit;
							}
					}else
						{
							$thread=1;
						}

		##score
				if (exists($opt{"score"})) 
					{
						$score=$opt{"score"};
						if ($score<=0) 
							{
								print "please assign an applicable score value.\n";
								exit;
							}
					}else
						{
							$score=40;
						}

		if ($method eq "GF")
			{
		###identity
				if (exists($opt{"identity"})) 
					{
						$identity=$opt{"identity"};
						if ($identity>1 or $identity<=0) 
							{
								print "identity should be 0 ~ 1 \n";
								exit;
							}
					}else
						{
							$identity=0.5;
						}

		###evalue
				if (exists($opt{"evalue"})) 
					{
						$evalue=$opt{"evalue"};
					}else
						{
							$evalue=1e-10;
						}

		###coverage
				if (exists($opt{"coverage"})) 
					{
						$coverage=$opt{"coverage"};
						if ($coverage>1 or $coverage<=0) 
							{
								print "coverage should be 0 ~ 1 \n";
								exit;
							}
					}else
						{
							$coverage=0.5;
						}
			}


		if ($method eq "MP")
			{
		###global
				if (exists($opt{"global"})) 
					{
						$global=$opt{"global"};
						if ($global>1 or $global<=0) 
							{
								print "global coverage should be 0 ~ 1 \n";
								exit;
							}
					}else
						{
							$global=0.5;
						}
		###local
				if (exists($opt{"local"})) 
					{
						$local=$opt{"local"};
						if ($local<=0) 
							{
								print "local coverage should be 0 ~ [global coverage value] \n";
								exit;
							}
						if ($local>$global) 
							{
								print "local coverage should be less than global coverage!\n";
								exit;
							}
					}else
						{
							$local=0.25;
						}
			}
	}

if ($run_evolution) 
	{
		if (exists($opt{"bootstrap"})) 
			{
				$bootstrap=$opt{"bootstrap"};
				if ($bootstrap<=0) 
					{
						print "please assign an applicable bootstrap value.\n";
					}
			}else
				{
					$bootstrap=1;
				}
	}

print "Program begin at ".localtime()."\n";
print "The following are the parameters for current process:\n";
				print "Strains:                   ".join(",",@species)."\n";
				print "Input directory:           $inputDIR\n";
				print "Output directory:          $outputDIR\n";
if ($run_cluster) 
	{
				print "Cluster analysis:          yes\n";
				print "  Method:                  $method\n";
				print "  Thread:                  $thread\n";
		if ($method eq "GF") 
			{
				print "  E-value:                 $evalue\n";
				print "  Identity:                $identity\n";
				print "  Coverage:                $coverage\n";
				print "  Score:                   $score\n";
			}
		if ($method eq "MP") 
			{
				print "  Local:                   $local\n";
				print "  Global:                  $global\n";
				print "  Score:                   $score\n";
			}
	}else
		{
				print "Cluster analysis:          no\n";
		}
if ($run_pangenome) 
	{
				print "Pan-genome analysis:       yes\n";
	}else
		{
				print "Pan-genome analysis:       no\n";
		}
if ($run_variation) 
	{
				print "Variation analysis:        yes\n";
	}else
		{
				print "Variation analysis:        no\n";
		}
if ($run_evolution) 
	{
				print "Evolution analysis:        yes\n";
				print "  Bootstrap:               $bootstrap\n";
	}else
		{
				print "Evolution analysis:        no\n";
		}
if ($run_function) 
	{
				print "Function analysis:         yes\n";
	}else
		{
				print "Function analysis:         no\n";
		}

$section=$run_cluster.$run_pangenome.$run_variation.$run_evolution.$run_function;

###############################################
#	section 0) check input file and program
###############################################
if (!(-e $outputDIR)) 
	{
		system("mkdir $outputDIR");
	}
system("chmod +rw $outputDIR");

if (!(-w $outputDIR)) 
	{
		print "There is no WRITE permission in $outputDIR\n";
		exit;
	}
@tmp=qw();
&CheckInputFile(\@species,$inputDIR,$section,$method,\@tmp);
&CheckExtraProgram($section,$method,\@tmp);

if (scalar(@tmp)>0) 
	{
		open(R,">".$outputDIR."0.error.message");
		print R join("",@tmp)."\n";
		close(R);
		print "error!\nlog are saved in ${outputDIR}0.error.message\n";
		exit;
	}


############################################
#	section 1) cluster analysis
############################################

if ($run_cluster) 
	{
	print "\n\n############################################\n";
	print "#	section 1) cluster analysis\n";
	print "############################################\n\n\n";

#### cluster gene and return result to the array @clusters

		if ($method eq "MP") 
			{
				print "Begin cluster gene with MP method ...\n";
				&MP();
			}else
				{
					print "Begin cluster gene with GF method ...\n";
					&GF();
				}

#### output normal cluster format

	&FormatClusterOutPut(\@species,"${outputDIR}1.Orthologs_Cluster.txt",\@clusters);

#### Retrieve cluster

	&RetrieveClusterFromFile("${outputDIR}1.Orthologs_Cluster.txt",\@clusters);

##### gene distribution in each strains
		%tmpHash=();
		&GeneDistribution(\@clusters,\%tmpHash);

		open(R,">${outputDIR}1.Gene_Distribution_By_Conservation.txt");
		print R "SharedBy_Strains\t".join("\t",@species)."\n";

		for ($i=$spnum;$i>0;$i--) 
			{
				print R $i;
				for ($j=0;$j<$spnum;$j++) 
					{
						if (exists($tmpHash{$i."|".$j})) 
							{
								print R "\t".$tmpHash{$i."|".$j};
							}else
								{
									print R "\t0";
								}
					}
				print R "\n";
			}

		close(R);


	}else
		{
			print "Homologous gene clustering is skipped!\n";
		}


if ($run_pangenome) 
	{
	print "\n\n############################################\n";
	print "#	section 2) Pan-genome analysis\n";
	print "############################################\n\n\n";

#### Retrieve cluster
	&RetrieveClusterFromFile("${outputDIR}1.Orthologs_Cluster.txt",\@clusters);
	chomp(@clusters);

#### convert file into 0-1 matrix
	for ($line=0;$line<@clusters;$line++) 
		{
			@row=split(/\t/,$clusters[$line]);
			splice(@row,0,1);
			for ($i=0;$i<@row;$i++) 
				{
					if ($row[$i] eq "-") 
						{
							$row[$i]=0;
						}else
							{
								$row[$i]=1;
							}
				}
			$clusters[$line]=join("\t",@row);
		}

#### fetch gene number of each strains
for ($i=0;$i<$spnum;$i++) 
	{
		open(F,"$inputDIR$species[$i].pep");
		@tmp=<F>;
		close(F);
		@tmp=grep(/^>/,@tmp);
		$genenum{$species[$i]}=scalar(@tmp);
	}

#### pan genome size and core genome size
print "Deducing pan genome size and core genome size for each composition...\n\n";

open(PAN,">${outputDIR}2.PanGenome.Data.txt");
print PAN "ClusterConservation\tTotalGeneNumber\tPanGenome\tCoreGenome\n";

for ($i=1;$i<=scalar(@species);$i++) 
	{
		#@SpecieCombination=&Combination(\@species,$i);
		#@SpecieCombination=&Combination($spnum,$i);
		if (&ChkCombinationValue($spnum,$i) !=0) ### transfer the array reference to the subroutine
		{
			&Combination($spnum,$i,\@SpecieCombination); ## if the combination number is less than sampleSize, then fecth all, else sample
		}else
		{
			&SampleCombination($spnum,$i,\@SpecieCombination);
		}
		
		foreach $key (@SpecieCombination) 
			{
		##### count total gene number in current combination
					$tmp=0;
					@spID=split(/\t/,$key); #### speices id  in current combination
					foreach  (@spID) 
						{
							$tmp=$tmp+$genenum{$species[$_]};
						}
		##### scan pangenome and coregenome
					@tmp=split(/\t/,&PanGenomeNumber(\@spID));
					print PAN "$i\t$tmp\t".join("\t",@tmp)."\n";

			}
	}

close(PAN);

#### data fit

#### for model A

if ($spnum<3) 
	{
		print "There are $spnum strains. For pan-genome function fitting, at least 3 strains data are required.\n";
	}else
		{
		open(R,">${outputDIR}2.PanGenome.Profile.txt");
	##### genome number & pan-genome size
			@xdata=qw();
			@ydata=qw();
			&ReadData2Array("${outputDIR}2.PanGenome.Data.txt",\@xdata,0,\@ydata,2);
			&SumData(\@xdata,\@ydata,"mean");
			($fit_Rsquare, $fit_A, $fit_A_interval, $fit_B, $fit_C, $fit_C_interval)=&fit_model_A(\@xdata,\@ydata);
			print R "The relation bewteen genome number and pan-genome size\n\n";
			print R "Function model: y=A*x**B +C \n";
			print R "\ty denotes pan-genome size, x denotes genome number, and A, B, C are fitting parameters.\n\n";
			print R "Fitting result:\n";
			print R "\ty = $fit_A *x**$fit_B + $fit_C\n";
			print R "\tR-square = $fit_Rsquare\n";
			print R "\tA 95% confidence interval: ($fit_A - $fit_A_interval , $fit_A + $fit_A_interval)\n";
			print R "\tC 95% confidence interval: ($fit_C - $fit_C_interval , $fit_C + $fit_C_interval)\n\n\n\n\n";

	##### total gene number & pan-genome size
			#@xdata=qw();
			#@ydata=qw();
			#&ReadData2Array("${outputDIR}2.PanGenome.Data.txt",\@xdata,1,\@ydata,2);
			#&SumDataByMedian(\@xdata,\@ydata);
			#($fit_Rsquare, $fit_A, $fit_B, $fit_C)=&fit_model_A(\@xdata,\@ydata);
			#print R "The relation bewteen total gene number and pan-genome size\n\n";
			#print R "$fit_Rsquare, $fit_A, $fit_B, $fit_C\n";
			#print R "\ty = $fit_A *x**$fit_B + $fit_C   R-square = $fit_Rsquare\n";
			#print R "\tx: total gene number\n";
			#print R "\ty: pan-genome size\n\n\n\n\n";

	##### genome number & core genome
			@xdata=qw();
			@ydata=qw();
			&ReadData2Array("${outputDIR}2.PanGenome.Data.txt",\@xdata,0,\@ydata,3);
			&SumData(\@xdata,\@ydata,"mean");
			($fit_Rsquare, $fit_A, $fit_A_interval, $fit_B, $fit_C, $fit_C_interval)=&fit_model_B(\@xdata,\@ydata);
			print R "The relation bewteen genome number and core genome size\n\n";
			print R "Function model: y=A*exp(B*x) +C \n";
			print R "\ty denotes pan-genome size, x denotes genome number, and A, B, C are fitting parameters.\n\n";
			print R "Fitting result:\n";
			print R "\ty = $fit_A *exp($fit_B * x) + $fit_C   R-square = $fit_Rsquare\n";
			print R "\tR-square = $fit_Rsquare\n";
			print R "\tA 95% confidence interval: ($fit_A - $fit_A_interval , $fit_A + $fit_A_interval)\n";
			print R "\tC 95% confidence interval: ($fit_C - $fit_C_interval , $fit_C + $fit_C_interval)\n\n\n\n\n";
		close(R);
		}

	}


	############################################
	#	section 3) CDS variation analysis
	############################################

if ($run_variation) 
	{
		print "\n\n############################################\n";
		print "#	section 3) CDS variation analysis\n";
		print "############################################\n\n\n";

	#### Retrieve cluster
		&RetrieveClusterFromFile("${outputDIR}1.Orthologs_Cluster.txt",\@clusters);
		chomp(@clusters);

	 ## protein
		system("rm -rf *.pep");
		&PrepareFasta(\@species,$inputDIR,".pep"); ###prepare pep file
		system("cat *.pep > All.faa && rm -rf *.pep && mv  All.faa All.pep");
		&ReadSequenceInToHash("All.pep",\%pep);
	 ## nucleic
		system("rm -rf *.nuc");
		&PrepareFasta(\@species,$inputDIR,".nuc"); ###prepare nuc file
		system("cat *.nuc > All.ffn && rm -rf *.nuc && mv  All.ffn All.nuc");
		&ReadSequenceInToHash("All.nuc",\%nuc);

	## scanning SNP
		%nt4tree=();
		for ($i=0;$i<$spnum;$i++) 
			{
				$nt4tree{"S".$i}="";
			}

		open(VAR,">${outputDIR}3.CDS.variation.txt");
		print VAR "ClusterID\tStrains_Number\tGene_Number\tPosition\taaType\tntType\tntProfile\tVariation type\n";

		open(VA,">${outputDIR}3.CDS.variation.analysis.txt");
		print VA "ClusterID\tInDel Base\tNonsynonymous mutation\tSynonymous mutation\n";

		for ($line=0;$line<@clusters;$line++) 
			{
				@row=split(/\t|\,/,$clusters[$line]);
				$ClusterID=$row[0];
				splice(@row,0,1);
				@row=grep(/^S/,@row);
				if (scalar(@row) >=2) 
					{
						open(PEP,">$ClusterID.pep");
						open(NUC,">$ClusterID.nuc");
						foreach $key (@row) 
							{
								print PEP ">$key\n$pep{$key}\n";
								print NUC ">$key\n$nuc{$key}\n";
							}
						close(PEP);
						close(NUC);
						system("$mafft --quiet $ClusterID.pep > $ClusterID.pal");
						#system("perl ./pal2nal.pl $ClusterID.pal $ClusterID.nuc -output fasta > $ClusterID.nal");
						$tmp=&pal2nal("$ClusterID.pal","$ClusterID.nuc","$ClusterID.nal");
						if ($tmp == 0) 
							{
								system("rm -rf $ClusterID.*");
								next;
							}

						@tmp=&DetectSNP();
						if (scalar(@tmp)>0) 
							{
								print VA $ClusterID."\t".&VarAnalysis(\@tmp)."\n";
								print VAR join("",@tmp);
						### core orthologs
								@row=split(/\t/,$clusters[$line]);
								splice(@row,0,1);
								if ((&CountGeneInCluster(join("\t",@row)) ==$spnum) and (&CountSpeicesInCluster(join("\t",@row)) == $spnum) )
									{
										$count_tree++;
										%tmpHash=();
										foreach  (@row) 
											{
												$tmpHash{$_}="";
											}
										&RemoveHeadGap("$ClusterID.nal",\%tmpHash);
										&ExtractSNP4tree(\%tmpHash,\%nt4tree);
									}
							}
						system("rm -rf $ClusterID.*");
					}
			}
		close(VAR);
		close(VA);

		open(R,">${outputDIR}3.CDS.variation.for.evolution.txt");
		foreach $key (keys %nt4tree) 
			{
				$_=$key;
				s/s//gi;
				print R ">$species[$_]\n$nt4tree{$key}\n";
			}
		close(R);
		print $count_tree."\n\n";

###
		system("rm All.nuc All.pep");
	}else
		{
			print "CDS variation is skipped.\n";
		}

	############################################
	#	section 4) CDS variation analysis
	############################################

if ($run_evolution) 
	{
#### Retrieve cluster
		&RetrieveClusterFromFile("${outputDIR}1.Orthologs_Cluster.txt",\@clusters);
		chomp(@clusters);



##################
##
##  Distance based
##
##################

### caculate the distance between each two strains by cluster
		&ClusterProfil4Specie(\%tmpHash);   # caculate Clusters profile for each specie
		&DistanceMatrix($spnum,\%tmpHash);   # caculate distance matrix acoording to Clusters profile
	###output distance
		open(DIST,">${outputDIR}4.Species_Distance_Clusters_Based.txt");
		###header
		printf DIST "%5d", $spnum;
		print DIST "\n";
		foreach $i (0..($spnum-1)) 
			{
				$key="sp".$i."sp";
				printf DIST "%-10s",$key;
				foreach $j (0..($spnum-1)) 
					{
						printf DIST "  %8f",$tmpHash{$i."-".$j};
					}
				print DIST "\n";
			}
		close(DIST);

		%tmpHash=(); 

		### based on pan genome (distance)
		print "\nDraw pangenome based phylogenetic tree ...\n\n";

		&PanBasedTree("${outputDIR}4.Species_Distance_Clusters_Based.txt","${outputDIR}4.PanBased");

##################
##
##  SNP based
##
##################
		%tmpHash=();
		if (!(-e "${outputDIR}3.CDS.variation.for.evolution.txt")) 
			{
				print "Variation in core orthologs cluster is not found from ${outputDIR}3.CDS.variation.for.evolution.txt.\n";
				print "Maybe you have skipped CDS variation analysis.\n";
			}else
				{
					&ReadSequenceInToHash("${outputDIR}3.CDS.variation.for.evolution.txt",\%tmpHash);
					open(R,">mlst.aln");
					for ($i=0;$i<@species;$i++) 
						{
							print R ">sp${i}sp\n".$tmpHash{$species[$i]}."\n";
						}
					close(R);
					&fasta2phylip("mlst.aln","mlst.phylip");
					system("rm -rf mlst.aln");

					print "\nDraw SNP based phylogenetic tree ...\n\n";
					&SNPBasedTree("mlst.phylip","${outputDIR}4.SNPBased");
					system("rm -rf mlst.phylip")
				}

########
# replace speices name
########

			opendir(DIR,"${outputDIR}");
			@tmp=readdir(DIR);
			closedir(DIR);
			@tmp=grep(/^4/,@tmp);
			foreach $tmp (@tmp) 
				{
					&ReplaceName(\@species,"${outputDIR}$tmp");
				}
	}else
		{
			print "Evolution analysis is skipped.\n";
		}


	############################################
	#	section 4) Function analysis
	############################################

if ($run_function) 
	{
#### Retrieve cluster
		&RetrieveClusterFromFile("${outputDIR}1.Orthologs_Cluster.txt",\@clusters);
		chomp(@clusters);

#### prepare annotation file
		&PrepareTable(\@species,$inputDIR,".function"); ###prepare location file
		&ReadAnnotation(\@species,\%cog,\%description);


#### assign function 
		open(R,">${outputDIR}5.Orthologs_Cluster_Function.txt");
		print R "ClusterID\tConservation_Level\tCOG\tDescription\n";
		for ($i=0;$i<@clusters;$i++) 
			{
				@row=split(/\t/,$clusters[$i]);
				$ClusterID=$row[0];
				splice(@row,0,1);
				print R $ClusterID."\t".&CountSpeicesInCluster(join("\t",@row))."\t".&getCOG(\@row,\%cog)."\t".&getDescription(\@row,\%description)."\n";
			}
		close(R);

#### COG distribution

###Whole Clusters COG Distribution
&outputCOGStatistic("${outputDIR}5.Orthologs_Whole_Cluster_COG_Distribution.txt",&scanCOG("${outputDIR}5.Orthologs_Cluster_Function.txt",$spnum,1));

###Core Clusters COG Distribution
&outputCOGStatistic("${outputDIR}5.Orthologs_Core_Cluster_COG_Distribution.txt",&scanCOG("${outputDIR}5.Orthologs_Cluster_Function.txt",$spnum,$spnum));

###Dispensable Clusters COG Distribution
&outputCOGStatistic("${outputDIR}5.Orthologs_Dispensable_Cluster_COG_Distribution.txt",&scanCOG("${outputDIR}5.Orthologs_Cluster_Function.txt",($spnum-1),2));

###strains specifc Clusters COG Distribution
&outputCOGStatistic("${outputDIR}5.Orthologs_specifc_Cluster_COG_Distribution.txt",&scanCOG("${outputDIR}5.Orthologs_Cluster_Function.txt",1,1));

system("rm -rf *.function");

	}else
		{
			print "Function analysis is skipped.\n";
		}

sub outputCOGStatistic()
	{
		(my $file,my $subcogcount)=@_;
		my @cogcat=("J A K L B","D Y V T M N Z W U O","C G E F H I P Q","R S -");
		my @cogdesc=("INFORMATION STORAGE AND PROCESSING","CELLULAR PROCESSES AND SIGNALING","METABOLISM","POORLY CHARACTERIZED");
		my @subcogcat=qw(J A K L B D Y V T M N Z W U O C G E F H I P Q R S -);
		my @subcogdesc=("[J] Translation, ribosomal structure and biogenesis","[A] RNA processing and modification","[K] Transcription","[L] Replication, recombination and repair","[B] Chromatin structure and dynamics","[D] Cell cycle control, cell division, chromosome partitioning","[Y] Nuclear structure","[V] Defense mechanisms","[T] Signal transduction mechanisms","[M] Cell wall/membrane/envelope biogenesis","[N] Cell motility","[Z] Cytoskeleton","[W] Extracellular structures","[U] Intracellular trafficking, secretion, and vesicular transport","[O] Posttranslational modification, protein turnover, chaperones","[C] Energy production and conversion","[G] Carbohydrate transport and metabolism","[E] Amino acid transport and metabolism","[F] Nucleotide transport and metabolism","[H] Coenzyme transport and metabolism","[I] Lipid transport and metabolism","[P] Inorganic ion transport and metabolism","[Q] Secondary metabolites biosynthesis, transport and catabolism","[R] General function prediction only","[S] Function unknown","[-] Unclassified");
		my %subcogdesc;
		my $key;
		my @cog;
		my $i;
		my $cognum;

		for ($i=0;$i<@subcogcat;$i++) 
			{
				$subcogdesc{$subcogcat[$i]}=$subcogdesc[$i];
			}

		open(R,">$file");
		for ($i=0;$i<@cogcat;$i++) 
			{
				$cognum=0;
				foreach $key (split(" ",$cogcat[$i])) 
					{
						$cognum=$cognum+$$subcogcount{$key};
					}
				print R $cogdesc[$i]." ( ".$cognum." )\n";
				foreach $key (split(" ",$cogcat[$i])) 
					{
						printf R "%-6d   %s\n",$$subcogcount{$key},$subcogdesc{$key};
					}
				print R "\n";
				
			}
		close(R);

	}

sub scanCOG()
	{
		(my $file,my $max_orth,my $min_orth)=@_;
		my @row;
		my @subcogcat=qw(J A K L B D Y V T M N Z W U O C G E F H I P Q R S -);
		my %subcogcount;
		my $cog;
		my $key;

		foreach $key (@subcogcat) 
			{
				$subcogcount{$key}=0;
			}

		@subcogcat=qw(J A K L B D Y V T M N Z W U O C G E F H I P Q R S);

		open(F,"$file");
		$_=<F>;
		while (<F>) 
			{
				@row=split(/\t/,$_);
				if ($row[1]>=$min_orth and $row[1]<=$max_orth) 
					{
						if ($row[2] eq "-") 
							{
								$subcogcount{"-"}++;
							}else
								{
									$_=uc($row[2]);
									s/COG//gi;
									$cog=$_;
									foreach $key (@subcogcat) 
										{
											if ($cog=~/$key/) 
												{
													$subcogcount{$key}++;
												}
										}
								}
					}
			}
		close(F);

		return \%subcogcount;
	}




sub getCOG()
	{
		(my $data,my $coghash)=@_;
		my $cog="";
		my @cog;
		my $key;
		my %hash;
		my @gene=split(/\t|\,/,join("\t",@$data));
		@gene=grep(/^S/,@gene);

		foreach $key (@gene) 
			{
				if (($$coghash{$key} ne "-") and ($$coghash{$key} ne "")) 
					{
						$cog=$cog.",".$$coghash{$key};
					}
			}
		@cog=split(/,/,$cog);
		foreach $cog (@cog) 
			{
				if ($cog ne "") 
					{
						$hash{$cog}=1;
					}
			}

		$cog=join(",",(keys %hash));
		if ($cog eq "") 
			{
				$cog="-";
			}
		return $cog;
	}

sub getDescription()
	{
		(my $data,my $deschash)=@_;
		my $desc="";
		my $key;
		my @gene=split(/\t|\,/,join("\t",@$data));
		@gene=grep(/^S/,@gene);

		foreach $key (@gene) 
			{
				if ( ($$deschash{$key} ne "") and ($$deschash{$key} ne "-") and ($$deschash{$key}!~/hypothetical/)) 
					{
						$desc=$$deschash{$key};
					}
			}

		if ($desc eq "") 
			{
				$desc="hypothetical protein";
			}

		return $desc;
	}



sub ReadAnnotation()
	{
		(my $species,my $cog,my $description)=@_;
		my $i;
		my @row;

		for ($i=0;$i<@$species;$i++) 
			{
				open(F,"$$species[$i].function");
				while (<F>) 
					{
						chomp($_);
						@row=split(/\t/,$_);
						if (scalar(@row)>=2) 
							{
								$$cog{$row[0]}=$row[1];
							}else
								{
									$$cog{$row[0]}="-";
								}

						if (scalar(@row)>=3) 
							{
								$$description{$row[0]}=$row[2];
							}else
								{
									$$description{$row[0]}="hypothetical protein";
								}
					}
				close(F);
			}
	}


sub SNPBasedTree()
	{
		(my $infile,my $outfileprefix)=@_;
		my $tmpin=$infile;
		my $tmpout;


#### boootstrap
print "\n#### seqboot ...\n\n";
		open(R,">seqboot.cmd");
		#print R "$tmpin\n";
		print R "R\n";
		print R "$bootstrap\n";
		print R "Y\n";
		print R "1\n";
		close(R);

		system("cp $tmpin infile");
		system("$seqboot < seqboot.cmd");
		system("mv outfile 100dnaseq");
		system("rm -rf infile");
		system("rm seqboot.cmd");                # 100dnasseq

#### dnaml
print "\n#### dnaml ...\n\n";
		open(R,">dnaml.cmd");
		#print R "100dnaseq\n";
		print R "T\n";
		print R "25\n";
	if ($bootstrap>1) 
		{
		print R "M\n";
		print R "D\n";
		#print R "100\n";
		print R "$bootstrap\n";
		print R "1\n"; # Random number seed (must be odd)?
		print R "5\n"; # Number of times to jumble?
		}
		print R "Y\n";
		close(R);

		system("cp 100dnaseq infile");
		system("$dnaml < dnaml.cmd");
		system("rm -rf outfile");
		system("rm -rf infile");
		system("mv outtree 100dnaseqtree");     # 100dnaseq, 100dnaseqtree

#### consense
print "\n#### dnaml consense ...\n\n";

		open(R,">consense.cmd");
		#print R "100dnaseqtree\n";
		print R "Y\n";
		close(R);

		system("cp 100dnaseqtree intree");
		system("$consense < consense.cmd");
		system("mv outfile ${outfileprefix}.ML.outfile");
		system("mv outtree ${outfileprefix}.ML.tree");
		system("rm -rf infile");
		system("rm -rf 100dnaseqtree");                        # 100dnaseq

#### dnadist
print "\n#### dnadist ...\n\n";
		open(R,">dnadist.cmd");
		#print R "100dnaseq\n";
		print R "T\n";
		print R "25\n";
	if ($bootstrap>1) 
		{
		print R "M\n";
		print R "D\n";
		#print R "100\n";
		print R "$bootstrap\n";
		}
		print R "Y\n";
		close(R);

		system("cp 100dnaseq infile");
		system("$dnadist < dnadist.cmd");
		system("rm -rf 100dnaseq");
		system("rm -rf infile");
		system("mv outfile 100dnadist");                 # 100dnadist

#### Neighbor-joining tree
print "\n#### Neighbor-joining ...\n\n";
		open(R,">NJ.cmd");
	if ($bootstrap>1) 
		{
		#print R "100dnadist\n";
		print R "M\n";
		#print R "100\n";
		print R "$bootstrap\n";
		print R "1\n";
		}
		print R "Y\n";
		close(R);

		system("cp 100dnadist infile");
		system("$neighbor < NJ.cmd");
		system("mv outtree 100dnadistNJtree");
		system("rm outfile");
		system("rm -rf infile");
		system("rm -rf NJ.cmd");                      # 100dnadist,100dnadistNJtree

#### NJ-consense
print "\n#### NJ-consense ...\n\n";
		open(R,">NJ-consense.cmd");
		#print R "100dnadistNJtree\n";
		print R "Y\n";
		close(R);

		system("cp 100dnadistNJtree intree");
		system("$consense < NJ-consense.cmd");
		system("mv outfile ${outfileprefix}.Neighbor-joining.outfile");
		system("mv outtree ${outfileprefix}.Neighbor-joining.tree");
		system("rm -rf NJ-consense.cmd");
		system("rm -rf intree");
		system("rm -rf 100dnadistNJtree");


#### UPGMA tree
print "\n#### UPGMA ...\n\n";
		open(R,">UPGMA.cmd");
		#print R "100dnadist\n";
		print R "N\n";
	if ($bootstrap>1) 
		{
		print R "M\n";
		#print R "100\n";
		print R "$bootstrap\n";
		print R "1\n";
		}
		print R "Y\n";
		close(R);

		system("cp 100dnadist infile");
		system("$neighbor < UPGMA.cmd");
		system("mv outtree 100dnadistUPGMAtree");
		system("rm -rf outfile");
		system("rm -rf infile");
		system("rm -rf UPGMA.cmd");

#### UPGMA-consense
print "\n#### UPGMA-consense ...\n\n";
		open(R,">UPGMA-consense.cmd");
		#print R "100dnadistUPGMAtree\n";
		print R "Y\n";
		close(R);

		system("cp 100dnadistUPGMAtree intree");
		system("$consense < UPGMA-consense.cmd");
		system("mv outfile ${outfileprefix}.UPGMA.outfile");
		system("mv outtree ${outfileprefix}.UPGMA.tree");
		system("rm -rf UPGMA-consense.cmd");
		system("rm -rf 100dnadistUPGMAtree");
		system("rm -rf intree");

###CLEAN TMP FILE

		system("rm -rf *.cmd");
		system("rm -rf 100dnadist");
	}

sub PanBasedTree()
	{
		(my $infile,my $outfileprefix)=@_;
		my $tmpin;
		my $tmpout;

		$tmpin=$infile;
#### Neighbor-joining tree

		open(R,">NJ.cmd");
		#print R "$tmpin\n";
		print R "Y\n";
		close(R);

		system("cp $tmpin infile");
		system("$neighbor < NJ.cmd");
		system("mv outfile ${outfileprefix}.Neighbor-joining.outfile");
		system("mv outtree ${outfileprefix}.Neighbor-joining.tree");
		system("rm -rf  NJ.cmd");
		system("rm -rf infile");

#### UPGMA tree

		open(R,">UPGMA.cmd");
		#print R "$tmpin\n";
		print R "N\n";
		print R "Y\n";
		close(R);

		system("cp $tmpin infile");
		system("$neighbor < UPGMA.cmd");
		system("mv outfile ${outfileprefix}.UPGMA.outfile");
		system("mv outtree ${outfileprefix}.UPGMA.tree");
		system("rm -rf UPGMA.cmd");
		system("rm -rf infile");


###CLEAN TMP FILE
		system("rm -rf *.cmd");
	}

sub DistanceMatrix()
	{
		(my $spnum,my $hash)=@_;
		my $i;
		my $j;
		my $k;
		my $dist;
		my $ref;
		my $query;
		foreach $i (0..($spnum-1)) 
			{
				foreach $j ($i..($spnum-1)) 
					{
						$ref=$$hash{$i};
						$query=$$hash{$j};
						$dist=0;
						for ($k=0;$k<length($ref);$k++) 
							{
								if (substr($ref,$k,1) ne substr($query,$k,1)) 
									{
										$dist++;
									}
							}
						$$hash{$i."-".$j}=$dist;
						$$hash{$j."-".$i}=$dist;
					}
			}
	}


sub ClusterProfil4Specie
	{
		(my $hash)=@_;
		my @row;
		my $i;

		foreach  (0..($spnum-1))   #initialization Hash
			{
				$$hash{$_}="";
			}

		foreach  (@clusters) 
			{
				@row=split(/\t/,$_);
				splice(@row,0,1);
				if (&CountSpeicesInCluster(join("\t",@row))>1) 
					{
					for ($i=0;$i<@row;$i++) 
						{
							if ($row[$i] eq "-") 
								{
									$$hash{$i}=$$hash{$i}."0";
								}else
									{
										$$hash{$i}=$$hash{$i}."1";
									}
						}
					}
			}
	}

# &ExtractSNP4tree(\%tmpHash,\%nt4tree);

sub ExtractSNP4tree()
	{
		(my $hash,my $nt4treeRef)=@_;
		my $key;
		my @row;
		my $i;
		my $len;
		my @tribases;
		foreach $key (keys %$hash) 
			{
				$$hash{substr($key,0,index($key,"G"))}=$$hash{$key};
				delete($$hash{$key});
			}

		for ($i=0;$i<$spnum;$i++)
		{
			$nt4tree{"S".$i}=$nt4tree{"S".$i}.$$hash{"S".$i};
		}
	}


=pod
sub ExtractSNP4tree()
	{
		(my $hash,my $nt4treeRef)=@_;
		my $key;
		my @row;
		my $i;
		my $len;
		my @tribases;
		foreach $key (keys %$hash) 
			{
				$$hash{substr($key,0,index($key,"G"))}=$$hash{$key};
				delete($$hash{$key});
			}
		@_=(keys %$hash);
		$len=length($_[0]);
		for ($j=0;3*$j<$len;$j++) 
			{
##### scanning each codon
				for ($i=0;$i<$spnum;$i++) 
					{
						$tribases[$i]=substr($$hash{"S".$i},3*$j,3);
					}
##### checking each codon
				if (&IsTheSame(@tribases) ==0) 
					{
					for ($i=0;$i<@tribases;$i++) 
						{
							$nt4tree{"S".$i}=$nt4tree{"S".$i}.$tribases[$i];
						}
					}
			}
	}
=cut


sub pal2nal()
	{
		(my $pal,my $nuc, my $nal)=@_;
		my %aaAln=();
		my %ffn=();
		my %ntAln=();
		my %nt;
		my $dna;
		my $nt;
		my $key;
		my $flag=1;
		my $i=0;
		my $j;

### read protein aligment result
		&ReadAlignmentToHash("$pal",\%aaAln);
### read nt sequences
		&ReadSequenceInToHash("$nuc",\%ffn);
		foreach $key (keys %ffn) 
			{
				$dna=$ffn{$key};
				#if (int(length($nt{$key})/3)*3 ne length($nt{$key})) 
				if (int(length($dna)/3)*3 ne length($dna)) 
					{
						$flag=0;
						print "The length of nucleotide sequence is not 3 integer times.\n";
						last;
					}else
						{
							for ($i=0;$i<(length($dna)/3);$i++) 
								{
									$nt{$key."|".$i}=substr($dna,$i*3,3);
								}
						}
			}

		if ($flag==0) 
			{
				return 0;
			}else
				{
					foreach $key (keys %aaAln)  ### replace aa with corresponding nt
						{
							$nt="";
							$i=0;
							for ($j=0;$j<length($aaAln{$key});$j++) 
								{
									if (substr($aaAln{$key},$j,1) eq "-") 
										{
											$nt=$nt."---";
										}else
											{
												$nt=$nt.$nt{$key."|".$i};
												$i++;
											}
								}
							$ntAln{$key}=$nt;
						}

					### output 
						open(R,">$nal");
						foreach  (keys %ntAln) 
							{
								print R ">$_\n".$ntAln{$_}."\n";
							}
						close(R);

					return 1;
				}
	}




sub DetectSNP()
	{
		my %faa;
		my %ffn;
		my @row;
		my $count_gene;
		my $count_sp;
		my @genelist;
		my $i;
		my $j;
		my $pepalnlen;
		my @cdsvar=qw();
		my $cdi=0;
		my @tribases;
		my @bases;
		my @aa;


### fetch gene list
		open(F,"$ClusterID.pep");
		@genelist=<F>;
		close(F);
		@genelist=grep(/^>/,@genelist);
		chomp(@genelist);
		$_=join("\t",@genelist);
		s/>//g;
		@genelist=split(/\t/,$_);

### count gene number and species number
		@row=split(/\t/,$clusters[$ClusterID-1]);
		splice(@row,0,1);
		$count_sp=&CountSpeicesInCluster(join("\t",@row));
		$count_gene=&CountGeneInCluster(join("\t",@row));

### read alignment sequences
		&ReadAlignmentToHash("$ClusterID.pal",\%faa);
		&ReadAlignmentToHash("$ClusterID.nal",\%ffn);

@_=(keys %faa);
$pepalnlen=length($faa{$_[0]});
### scan SNP
		for ($i=1;$i<=$pepalnlen;$i++) 
			{
				@tmp=qw();
				@tribases=qw();
				for ($j=0;$j<@genelist;$j++) ### fetch triplet codon
					{
						$tribases[$j]=substr($ffn{$genelist[$j]},3*($i-1),3);
					}
				if (&IsTheSame(@tribases) ==0) ### if triplet codon is not consistent
					{
						@aa=qw();
						for ($j=0;$j<@genelist;$j++) 
							{
								$aa[$j]=substr($faa{$genelist[$j]},($i-1),1);
							}
						if (&IsTheSame(@aa) ==0) ### aa is not consistent
							{
								if (join("",@aa) =~/-/) 
									{
										$cdsvar[$cdi++]=$ClusterID."\t".$count_sp."\t".$count_gene."\t".$i."\t".&CharType(\@aa)."\t-\t-\tInDel\n";
									}else
										{
									#### base 1
											for ($j=0;$j<@genelist;$j++) 
												{
													$bases[$j]=substr($ffn{$genelist[$j]},3*($i-1),1);
												}
											if (&IsTheSame(@bases) ==0) 
												{
													$cdsvar[$cdi++]=$ClusterID."\t".$count_sp."\t".$count_gene."\t".($i+0.1)."\t".&CharType(\@aa)."\t".&CharType(\@bases)."\t".join("",@bases)."\tNonsynonymous mutation\n";
												}
									#### base 2
											for ($j=0;$j<@genelist;$j++) 
												{
													$bases[$j]=substr($ffn{$genelist[$j]},3*($i-1)+1,1);
												}
											if (&IsTheSame(@bases) ==0) 
												{
													$cdsvar[$cdi++]=$ClusterID."\t".$count_sp."\t".$count_gene."\t".($i+0.2)."\t".&CharType(\@aa)."\t".&CharType(\@bases)."\t".join("",@bases)."\tNonsynonymous mutation\n";
												}
									#### base 3
											for ($j=0;$j<@genelist;$j++) 
												{
													$bases[$j]=substr($ffn{$genelist[$j]},3*($i-1)+2,1);
												}
											if (&IsTheSame(@bases) ==0) 
												{
													$cdsvar[$cdi++]=$ClusterID."\t".$count_sp."\t".$count_gene."\t".($i+0.3)."\t".&CharType(\@aa)."\t".&CharType(\@bases)."\t".join("",@bases)."\tNonsynonymous mutation\n";
												}
										}
							}else
								{
									#### base 1
											for ($j=0;$j<@genelist;$j++) 
												{
													$bases[$j]=substr($ffn{$genelist[$j]},3*($i-1),1);
												}
											if (&IsTheSame(@bases) ==0) 
												{
													$cdsvar[$cdi++]=$ClusterID."\t".$count_sp."\t".$count_gene."\t".($i+0.1)."\t".&CharType(\@aa)."\t".&CharType(\@bases)."\t".join("",@bases)."\tSynonymous mutation\n";
												}
									#### base 2
											for ($j=0;$j<@genelist;$j++) 
												{
													$bases[$j]=substr($ffn{$genelist[$j]},3*($i-1)+1,1);
												}
											if (&IsTheSame(@bases) ==0) 
												{
													$cdsvar[$cdi++]=$ClusterID."\t".$count_sp."\t".$count_gene."\t".($i+0.2)."\t".&CharType(\@aa)."\t".&CharType(\@bases)."\t".join("",@bases)."\tSynonymous mutation\n";
												}
									#### base 3
											for ($j=0;$j<@genelist;$j++) 
												{
													$bases[$j]=substr($ffn{$genelist[$j]},3*($i-1)+2,1);
												}
											if (&IsTheSame(@bases) ==0) 
												{
													$cdsvar[$cdi++]=$ClusterID."\t".$count_sp."\t".$count_gene."\t".($i+0.3)."\t".&CharType(\@aa)."\t".&CharType(\@bases)."\t".join("",@bases)."\tSynonymous mutation\n";
												}
								}
					}
			}
		return @cdsvar;
	}


sub VarAnalysis()
	{
		(my $data)=@_;
		my @data=@$data;
		my $indel=0;
		my $syn=0;
		my $nonsyn=0;
		my @tmp;
		$indel=scalar(grep(/InDel$/,@data));
		$nonsyn=scalar(grep(/Nonsynonymous mutation$/,@data));;
		$syn=scalar(grep(/Synonymous mutation$/,@data));
		return "$indel\t$nonsyn\t$syn";
	}



sub CharType()
	{
		(my $str)=@_;
		my %hash;
		my @data=@$str;
		foreach  (@data) 
			{
				$hash{$_}=1;
			}
		return join(",",(keys %hash));
	}

sub IsTheSame()
	{
		(my @data)=@_;
		my %hash;
		foreach  (@data) 
			{
				$hash{$_}=1;
			}
		if (scalar(keys %hash) ==1) 
			{
				return 1;
			}else
				{
					return 0;
				}
	}



sub FormatClusterOutPut()
	{
		(my $speices,my $file,my $cluster)=@_;
		my @row;
		my $gid=1;
		my $key;
		my %hash;
		my $gene;
		my @tmp;
		my $i;
		my $j;
		open(R,">$file");
		print R "ClutserID\t".join("\t",@$speices)."\n";
		foreach $key (@$cluster) 
			{
				@row=split(/\t/,$key);
				for ($i=0;$i<@row;$i++) 
					{
						if ($row[$i] ne "-") 
							{
								@tmp=split(/,/,$row[$i]);
								for ($j=0;$j<@tmp;$j++) 
									{
										$_=$tmp[$j];
										s/^S[0-9]+G//;
										$tmp[$j]=$_;
									}
								$row[$i]=join(",",@tmp);
							}
					}
				print R $gid."\t".join("\t",@row)."\n";
				$gid++;
			}
		close(R);
	}
sub RetrieveClusterFromFile()
	{
		(my $file,my $clusters)=@_;
		my @content;
		my @row;
		my $spid;
		my $line=0;
		my $i=0;
		my $j;
		my @tmp;
		open(F,$file) or die "Could open $file\n";
		@content=<F>;
		close(F);
		splice(@content,0,1);
		chomp(@content);
		foreach  (@content) 
			{
				@row=split(/\t/,$_);
				$$clusters[$line]=$row[0];
				splice(@row,0,1);
				for ($i=0;$i<@row;$i++) 
					{
						if ($row[$i] ne "-") 
							{
								@tmp=split(/,/,$row[$i]);
								for ($j=0;$j<@tmp;$j++) 
									{
										$tmp[$j]="S${i}G".$tmp[$j];
									}
								$row[$i]=join(",",@tmp);
							}
					}
				$$clusters[$line]=$$clusters[$line]."\t".join("\t",@row)."\n";
				$line++;
			}
	}




sub GeneDistribution()
	{
		(my $clusters,my $hash)=@_;
		my @row;
		my $spid;
		my $orth;
		my $key;
		foreach (@$clusters) 
			{
				@row=split(/\t/,$_);
				splice(@row,0,1);
				$orth=&CountSpeicesInCluster(join("\t",@row));
				@row=split(/\t|\,/,join("\t",@row));
				foreach $key (@row) 
					{
						if ($key ne "-") 
							{
							$spid=substr($key,1,(index($key,'G')-1)); ###extract strains id
							if (exists($$hash{$orth."|".$spid})) 
								{
									$$hash{$orth."|".$spid}++;
								}else
									{
										$$hash{$orth."|".$spid}=1;
									}
							}
					}
			}
	}

sub CountSpeicesInCluster()
	{
		(my $str)=@_;
		chomp($str);
		my @list=split(/\t/,$str);
		my $key;
		my $count=0;

		foreach $key (@list) 
			{
				if ($key ne "-") 
					{
						$count++;
					}
			}
		return $count;
	}

sub CountGeneInCluster()
	{
		(my $str)=@_;
		chomp();
		my @list=split(/\t|\,/,$str);
		my $key;
		my $count=0;
		foreach $key (@list) 
			{
				if ($key ne "-") 
					{
						$count++;
					}
			}
		return $count;
	}




sub GF()
	{
		&PrepareFasta(\@species,$inputDIR,".pep"); ###prepare pep file
		system("cat ".join(".pep ",@species).".pep > All.pep");
		system("grep '>' All.pep > genelist");
		system("$formatdb -p T -i All.pep");
		system("$blastall -p blastp -i All.pep -d All.pep -M BLOSUM45 -m9 -e $evalue -o All.blastp -a $thread");
		system("perl ./Blast_Filter.pl All.blastp All.pep $coverage $identity $score | $mcl - --abc -I 2.0 -o All.cluster");
		&FormatCluster("All.cluster","genelist",$spnum,\@clusters);
		#system("rm -rf *.pep* All.blastp All.cluster genelist");
	}

sub MP()
	{
#		(my $species,my $inputDIR,my $thread,my $evalue,my $score,my $coverage,my $identity)=@_;
		my $i;
		my $j;
		&PrepareFasta(\@species,$inputDIR,".pep"); ###prepare pep file
		system("cat ".join(".pep ",@species).".pep > All.pep");
		system("grep '>' All.pep > genelist");
		system("rm -rf All.pep");
		for ($i=0;$i<$spnum;$i++) 
			{
				for ($j=$i+1;$j<$spnum;$j++) 
					{
						system("perl ./inparanoid.pl $blastall $thread $formatdb $score $global $local $species[$i].pep $species[$j].pep");
					}
			}
		system("perl ./multiparanoid.pl -species ".join(".pep+",@species).".pep -unique 1");
###convert the MP result to table list based on gene
		&MP_Result_to_Table("MP.Cluster","All.cluster");
		&FormatCluster("All.cluster","genelist",$spnum,\@clusters);
		system("rm -rf sqltable.* *.pep* MP.Cluster genelist");
	}

sub fasta2phylip()
	{
		(my $input,my $output)=@_;
		use Bio::AlignIO;
		my $inputfilename = "10.aln";
		my $in= Bio::AlignIO->new(-file   => $input ,
		                          -format => 'fasta');
		my $out = Bio::AlignIO->new(-file   => ">$output" ,
		                          -format => 'phylip');
		while ( my $aln = $in->next_aln() ) 
			{
			$out->write_aln($aln);
			}
	}


sub RemoveHeadGap()
	{
		(my $nal,my $hash)=@_;
		my %aln;
		my $key;
		my $gaplength=0;
		my $len1;
		my $len2;
		&ReadSequenceInToHash("$nal",\%aln);
		foreach $key (keys %aln) 
			{
				$len1=length($aln{$key});
				$_=$aln{$key};
				s/^-+//;
				$len2=length($_);
				if (($len1-$len2)>$gaplength) 
					{
						$gaplength=$len1-$len2;
					}
			}
		foreach $key (keys %aln) 
			{
				$$hash{$key}=$$hash{$key}.substr($aln{$key},$gaplength,(length($aln{$key})-$gaplength));
			}
	}

sub PrepareFasta()
	{
		(my $species,my $inputDIR,my $extention)=@_;
		my $sp;
		my $file;
		my $i;
		my %hash;
		my $key;
		for ($i=0;$i<@$species;$i++) 
			{
				$file=$inputDIR.$$species[$i].$extention;
				%hash=();
				&ReadSequenceInToHash($file,\%hash);
				open(R,">$$species[$i]${extention}") or die "Could write into $file\n";
				foreach $key (keys %hash) 
					{
						print R ">S${i}G$key\n";
						print R $hash{$key}."\n";
					}
				close(R);
			}
	}

sub PrepareTable()
	{
		(my $species,my $inputDIR,my $extention)=@_;
		my @content;
		my $i;
		my @row;
		my $file;
		for ($i=0;$i<@$species;$i++) 
			{
				$file=$inputDIR.$$species[$i].$extention;
				open(F,$file) or die "Could open $file\n";
				@content=<F>;
				close(F);
				chomp(@content);
				open(R,">$$species[$i]${extention}") or die "Could write into $file\n";
				foreach  (@content) 
					{
						@row=split(/\t/,$_);
						$row[0]="S${i}G$row[0]";
						if ($extention eq ".location") 
							{
								$row[0]=$row[0]."\t".$row[0];
							}
						print R join("\t",@row)."\n";
					}
				close(R);
			}
	}

sub CheckExtraProgram
	{
		#(my $section, my $method, my $tmparray)=@_;
		my @error;
		my $ei=0;

#####cluster gene 
		if (substr($section,0,1) eq "1") 
			{
###MP: blastall formatdb 
###GF: blastall formatdb mcl

			if (!(-e $formatdb)) 
				{
					$error[$ei++]="formatdb is not found at $formatdb\n";
				}

			if (!(-X $formatdb)) 
				{
					$error[$ei++]="there is not premission to execute $formatdb\n";
				}

			if (!(-e $blastall)) 
				{
					$error[$ei++]="blastall is not found at $blastall\n";
				}

			if (!(-X $blastall)) 
				{
					$error[$ei++]="there is not premission to execute $blastall\n";
				}

				if ($method eq "GF") 
					{
						if (!(-e $mcl)) 
							{
								$error[$ei++]="mcl is not found at $mcl\n";
							}
						if (!(-X $mcl)) 
							{
								$error[$ei++]="there is not premission to execute $mcl\n";
							}
					}
			}

#####CDS variation
		if (substr($section,2,1) eq "1") 
			{
				if (!(-e $mafft)) 
					{
						$error[$ei++]="mafft is not found at $mafft\n";
					}
				if (!(-X $mafft)) 
					{
						$error[$ei++]="there is not premission to execute $mafft\n";
					}
			}

#####CDS variation
		if (substr($section,3,1) eq "1") 
			{
				if (!(-e $mafft)) 
					{
						$error[$ei++]="mafft is not found at $mafft\n";
					}
				if (!(-X $mafft)) 
					{
						$error[$ei++]="there is not premission to execute $mafft\n";
					}
			}
#####Evolution analysis
		if (substr($section,3,1) eq "1") 
			{
				if (-e $seqboot) 
					{
						$error[$ei++]="there is not premission to execute $seqboot\n" if(!(-X $seqboot));
					}else
						{
							$error[$ei++]="seqboot is not found at $seqboot\n";
						}
				if (-e $dnaml) 
					{
						$error[$ei++]="there is not premission to execute $dnaml\n" if(!(-X $dnaml));
					}else
						{
							$error[$ei++]="dnaml is not found at $dnaml\n";
						}
				if (-e $dnadist) 
					{
						$error[$ei++]="there is not premission to execute $dnadist\n" if(!(-X $dnadist));
					}else
						{
							$error[$ei++]="dnadist is not found at $dnadist\n";
						}
				if (-e $neighbor) 
					{
						$error[$ei++]="there is not premission to execute $neighbor\n" if(!(-X $neighbor));
					}else
						{
							$error[$ei++]="neighbor is not found at $neighbor\n";
						}
				if (-e $consense) 
					{
						$error[$ei++]="there is not premission to execute $consense\n" if(!(-X $consense));
					}else
						{
							$error[$ei++]="consense is not found at $consense\n";
						}
				if (-e $dnapars) 
					{
						$error[$ei++]="there is not premission to execute $dnapars\n" if(!(-X $dnapars));
					}else
						{
							$error[$ei++]="dnapars is not found at $dnapars\n";
						}
			}
		#@$tmparray=(@$tmparray,@error);
		@tmp=(@tmp,@error);
	}


sub CheckInputFile()
	{
		(my $species,my $inputDIR,my $section,my $method,my $tmparray)=@_;
####cluster
		if (substr($section,0,1) eq "1") 
			{
				if ($method eq "MM") 
					{
#						@$tmparray=(@$tmparray,&chk2SEQ($species,$inputDIR)); ### check pep and nuc
						@$tmparray=(@$tmparray,&chktab($species,$inputDIR,".location"));### chk  pep nuc location
					}else
						{
						@$tmparray=(@$tmparray,&chk1SEQ($species,$inputDIR));
						}
			}
###CDS variation
		if (substr($section,2,1) eq "1")
			{
#				@$tmparray=(@$tmparray,&chk2SEQ($species,$inputDIR));
			}
###function analysis
		if (substr($section,4,1) eq "1")
			{
				@$tmparray=(@$tmparray,&chktab($species,$inputDIR,".function"));
			}
	}


sub chk1SEQ()
	{
		(my $species,my $inputDIR)=@_;
		my @error;
		my $ei=0;
		my $sp;
		my $pepfile;
		my %pep;
		foreach $sp (@$species) 
			{
				%pep=();
				$pepfile=$inputDIR.$sp.".pep";
				&ReadSequenceInToHash($pepfile,\%pep);
				if (scalar(keys %pep)<2) 
					{
						$error[$ei++]="format error in $pepfile\n";
					}
			}
		return @error;
	}

sub chk2SEQ()
	{
		(my $species,my $inputDIR)=@_;
		my $sp;
		my %pep;
		my %nuc;
		my $pepfile;
		my $nucfile;
		my $key;
		my @error;
		my $ei=0;
		foreach $sp (@$species) 
			{
				$pepfile=$inputDIR.$sp.".pep";
				$nucfile=$inputDIR.$sp.".nuc";
				%pep=();
				%nuc=();
				&ReadSequenceInToHash("$pepfile",\%pep);
				&ReadSequenceInToHash("$nucfile",\%nuc);
				if (scalar(keys %pep) ne scalar(keys %nuc)) 
					{
						$error[$ei++]="Sequences number is not consistent in the following two file:\n\t$pepfile\n\t$nucfile\n";
					}else
						{
							foreach $key (keys %pep) 
								{
									if (exists($nuc{$key})) 
										{
											if (length($nuc{$key}) ne ((length($pep{$key})+1)*3)) 
												{
													$error[$ei++]="the length of $key in $nucfile is not consistent with its corresponding protein length\n";
												}
										}else
											{
												$error[$ei++]="$key lost in $nucfile\n";
											}
								}

							foreach $key (keys %nuc) 
								{
									if (!exists($pep{$key})) 
										{
											$error[$ei++]="1048 $key lost in $pepfile\n";
										}
								}
						}
			}
		return @error;
	}


sub chktab()
	{
		(my $species,my $inputDIR,my $extention)=@_;
		my %pep;
		my @row;
		my $key;
		my %tab;
		my @error;
		my $ei=0;
		my $sp;
		my $tabfile;
		my $pepfile;
		foreach $sp (@$species) 
			{
				%tab=();
				%pep=();
				$tabfile=$inputDIR.$sp.$extention;
				open(F,"$tabfile");
				while (<F>) 
					{
						chomp();
						@row=split(/\t/,$_);
						if (scalar(@row)<3) 
							{
								$error[$ei++]="format error in $tabfile\n";
							}else
								{
									$tab{$row[0]}=$row[1];
								}
					}
				close(F);
				$pepfile=$inputDIR.$sp.".pep";
				&ReadSequenceInToHash($pepfile,\%pep);
				foreach $key (keys %pep) 
					{
						if (!exists($tab{$key})) 
							{
								$error[$ei++]="sequence $key lost infomation in $tabfile\n";
							}
					}
			}
		return @error;
	}



sub ReadSequenceInToHash()
	{
		use Bio::SeqIO;
		(my $file,my $hash)=@_;
		my $seq;
		my $in=Bio::SeqIO->new(-file=>"$file",-format=>"fasta");
		while ($seq=$in->next_seq()) 
			{
				#$$hash{$id."|".$seq->id}=$seq->seq();
				$$hash{$seq->id}=$seq->seq();
			}
	}

sub ReadAlignmentToHash()
	{
		(my $file,my $hash)=@_;
		my $name="";
		my $seq="";
		my @content;
		my $line;
		open(F,"$file");
		@content=<F>;
		close(F);
		chomp(@content);
		for ($line=0;$line<@content;$line++) 
			{
				if ($content[$line]=~/^>/) 
					{
						if ($line>0) 
							{
								$$hash{$name}=$seq;
								$name="";
							}

						$_=$content[$line];
						s/^>//;
						$name=$_;
						$seq="";
					}else
						{
							if ($name ne "") 
								{
								$seq=$seq.$content[$line];
								}
						}
			}
		$$hash{$name}=$seq;
	}



sub Combination()
	{
		(my $m,my $n,my $comRef)=@_;
		my $str="";
		my %hash;
		my $fpos;
		my $num0;
		my $rest;
		my $tmp;
		my $i;
		my $j;
		my $key;
		#my $m=scalar(@$array);
		my @combination;

		for ($i=1;$i<=$n;$i++) 
			{
				$str="1".$str;
			}

		for ($i=1;$i<=($m-$n);$i++) 
			{
				$str=$str."0";
			}

		$hash{$str}=1;
		while ($str=~/10/) 
			{
				$fpos=index($str,"10");
				$_=$str;
				s/10/01/;
				$str=$_;
				$tmp=substr($str,0,$fpos);
				$_=$tmp;
				s/0//g;
				$rest=$_;
				$num0=$fpos-length($_);
				for ($i=1;$i<=$num0;$i++) 
					{
						$rest="$rest"."0";
					}
				$str="$rest".substr($str,$fpos,$m-$fpos);
				$hash{$str}=1;
			}
		$j=0;
		foreach $key (keys %hash) 
			{
				$combination[$j]="";
				for ($i=0;$i<$m;$i++) 
					{
						if (substr($key,$i,1) eq "1") 
							{
							if ($combination[$j] ne "") 
								{
									#$combination[$j]=$combination[$j]."\t".$$array[$i];
									$combination[$j]=$combination[$j]."\t".$i;
								}else
									{
										#$combination[$j]=$$array[$i];  ### For return species ID
										$combination[$j]=$i;
									}
							}
					}
				$j++;
			}
		@$comRef=@combination; ### update the data through the physic address
	}

sub ChkCombinationValue()
{
	(my $m,my $n)=@_;
	my %hash;
	my %vhash;
	my $value=0;
	my $key;
	my @row;
	my @sdA;
	my @sdB;

	### initialization
	$hash{$m."-".$n}=1;

	### split combination
	while (scalar(keys %hash)>0 and $value<=$sampleSize)
	{
		foreach $key (keys %hash)
		{
			if ($value > $sampleSize) ### threshold
			{
				last;
			}
			if (!exists($hash{$key}))
			{
				next;
			}
			@row=split(/-/,$key);
			#print $row[0]."|".$row[1]."\n";
			if ($row[0] eq $row[1])
			{
				$value=$value+$hash{$key};
			}else
			{
				##split
				$sdA[0]=$row[0]-1;
				$sdA[1]=$row[1];
				$sdB[0]=$row[0]-1;
				$sdB[1]=$row[1]-1;
				##storing A
				if (($sdA[0] eq $sdA[1]) or $sdA[1] ==0)
				{
					$value=$value+$hash{$key};
				}else
				{
					if (exists($hash{$sdA[0]."-".$sdA[1]}))
					{
						$hash{$sdA[0]."-".$sdA[1]}=$hash{$sdA[0]."-".$sdA[1]}+$hash{$key};
					}else
					{
						$hash{$sdA[0]."-".$sdA[1]}=$hash{$key};
					}
				}

				##storing B
				if (($sdB[0] eq $sdB[1]) or $sdB[1]==0)
				{
					$value=$value+$hash{$key};
				}else
				{
					if (exists($hash{$sdB[0]."-".$sdB[1]}))
					{
						$hash{$sdB[0]."-".$sdB[1]}=$hash{$sdB[0]."-".$sdB[1]}+$hash{$key};
					}else
					{
						$hash{$sdB[0]."-".$sdB[1]}=$hash{$key};
					}
				}
			}
			#delete original combination
			delete($hash{$key});
		}
	}

	if ($value>$sampleSize)
	{
		return 0;
	}else
	{
		return $value;
	}
}


sub SampleCombination()
{
	(my $m,my $n,my $comRef)=@_;
	my %hash;
	my $sampleTimes=0;
	my @randNum;
	my @sortID;
	my $i;
	my $j;
	my $tmp;
	while ( scalar(keys %hash)<$sampleSize and $sampleTimes<($sampleSize*2))
	{
		for ($i=0;$i<$m;$i++) # generate random data
		{
			$randNum[$i]=int(100000 * rand(100));
			$sortID[$i]=$i;
		}

		for ($i=0;$i<$m;$i++) # sorting random data
		{
			for ($j=0;$j<$m;$j++)
			{
				if ($randNum[$sortID[$i]]<$randNum[$sortID[$j]])
				{
					$tmp=$sortID[$i];
					$sortID[$i]=$sortID[$j];
					$sortID[$j]=$tmp;
				}
			}
		}

		#storing data
		$tmp=join("\t",sort {$a<=>$b} (splice(@sortID,0,$n)));
		$hash{$tmp}=1;
		$sampleTimes++;
	}
	@$comRef=keys %hash;
}


sub PanGenomeNumber()
	{
		(my $spID)=@_;
		my $pan=0;    
		my $core=0;
		my $count; #### counter;
		my @row;

		foreach  (@clusters) 
			{
				$count=0;
				@row=split(/\t/,$_);

				foreach  (@$spID) 
					{
						$count=$count+$row[$_];
					}
				
				if ($count>0) 
					{
						$pan++;
						if ($count == scalar(@$spID)) 
							{
								$core++;
							}
					}
			}
		return $pan."\t".$core;
	}

sub fit_model_A()
	{
### model y = A * x**B + C
		(my $xdata,my $ydata)=@_;
		my $i;
		my $b;
		my $max_B=0;
		my $max_R=0;
		my $max_A=0;
		my $max_A_interval;
		my $max_C=0;
		my $max_C_interval;
		my $R=1e-100;
		my $start;
		my $end;
		my $step;
		my @xValues;
		my @yValues;

		$start=1;
		$step=0.001;
		$b=$start;
		$max_R=0;
		$R=1e-100;

		use Statistics::LineFit;
		use Statistics::Distributions;

		while ($max_R<=$R) 
			{
				if (($b < 0.02) and ($b >-0.02)) 
					{
						$b=-0.02;
					}

				for ($i=0;$i<@$xdata;$i++) 
					{
						$xValues[$i]=$$xdata[$i]**$b;
					}
				@yValues=@$ydata;
				my $lineFit = Statistics::LineFit->new();
				$lineFit->setData (\@xValues, \@yValues) or die "Invalid data";
				(my $intercept, my $slope) = $lineFit->coefficients();
				my $rSquared = $lineFit->rSquared();
				my $meanSquaredError = $lineFit->meanSqError();
				my $durbinWatson = $lineFit->durbinWatson();
				my $sigma = $lineFit->sigma();
				(my $tStatIntercept, my $tStatSlope) = $lineFit->tStatistics();
				(my $varianceIntercept,my $varianceSlope) = $lineFit->varianceOfEstimates();

				$max_R=$R;
				$R=$rSquared;
				if ($max_R<=$R) 
					{
						$max_R=$R;
						($max_C,$max_A)=$lineFit->coefficients();
						$max_A_interval=Statistics::Distributions::tdistr (($spnum-2),.025)*sqrt($varianceSlope);
						$max_C_interval=Statistics::Distributions::tdistr (($spnum-2),.025)*sqrt($varianceIntercept);
					}
				$b=$b-$step;
			}
		$max_B=$b;
		return ($max_R,$max_A,$max_A_interval,$max_B,$max_C,$max_C_interval);

	}

sub fit_model_B()
	{
### model y = A * exp(x*B) + C
		(my $xdata,my $ydata)=@_;
		my $i;
		my $b;
		my $max_B=0;
		my $max_R=0;
		my $max_A=0;
		my $max_A_interval;
		my $max_C=0;
		my $max_C_interval;
		my $R=1e-100;
		my $start;
		my $end;
		my $step;
		my @xValues;
		my @yValues;

		$start=0;
		$step=0.001;
		$b=$start;
		$max_R=0;
		$R=1e-100;

		use Statistics::LineFit;
		use Statistics::Distributions;

		while ($max_R<=$R) 
			{
				if (($b < 0.02) and ($b >-0.02)) 
					{
						$b=-0.02;
					}

				for ($i=0;$i<@$xdata;$i++) 
					{
						$xValues[$i]=exp($$xdata[$i]*$b);
					}
				@yValues=@$ydata;
				my $lineFit = Statistics::LineFit->new();
				$lineFit->setData (\@xValues, \@yValues) or die "Invalid data";
				(my $intercept, my $slope) = $lineFit->coefficients();
				my $rSquared = $lineFit->rSquared();
				my $meanSquaredError = $lineFit->meanSqError();
				my $durbinWatson = $lineFit->durbinWatson();
				my $sigma = $lineFit->sigma();
				(my $tStatIntercept, my $tStatSlope) = $lineFit->tStatistics();
				(my $varianceIntercept,my $varianceSlope) = $lineFit->varianceOfEstimates();

				$max_R=$R;
				$R=$rSquared;
				if ($max_R<=$R) 
					{
						$max_R=$R;
						($max_C,$max_A)=$lineFit->coefficients();
						$max_A_interval=Statistics::Distributions::tdistr (($spnum-2),.025)*sqrt($varianceSlope);
						$max_C_interval=Statistics::Distributions::tdistr (($spnum-2),.025)*sqrt($varianceIntercept);
					}
				$b=$b-$step;
			}
		$max_B=$b;
		return ($max_R,$max_A,$max_A_interval,$max_B,$max_C,$max_C_interval);
	}


sub ReadData2Array()
	{
		(my $file, my $array1,my $col1,my $array2,my $col2)=@_;
		my $i=0;
		open(F,$file);
		$_=<F>;
		while (<F>) 
			{
				chomp();
				@_=split(/\t/,$_);
				$$array1[$i]=$_[$col1];
				$$array2[$i]=$_[$col2];
				$i++;
			}
		close(F);
	}

sub SumData()
	{
		(my $xdata,my $ydata,my $SumMethod)=@_;
		my %hash;
		my $i;
		my $key;
		my $max=0;
		for ($i=0;$i<@$xdata;$i++) 
			{
				if (exists($hash{$$xdata[$i]})) 
					{
						$hash{$$xdata[$i]}=$hash{$$xdata[$i]}." ".$$ydata[$i];
					}else
						{
							$hash{$$xdata[$i]}=$$ydata[$i];
							if ($$xdata[$i]>$max) 
								{
									$max=$$xdata[$i];
								}
						}
			}
		@$xdata=qw();
		@$ydata=qw();
		$i=0;
		foreach $i (1..$max) 
			{
				$$xdata[$i-1]=$i;
				if ($SumMethod eq "median") 
					{
					$$ydata[$i-1]=&median($hash{$i});
					}else
						{
							$$ydata[$i-1]=&mean($hash{$i});
						}
			}
		#print join(",",@$xdata)."\n";
		#print join(",",@$ydata)."\n";
	}

sub median()
	{
		(my $data)=@_;
		my @data=split(/ /,$data);
		my $arraylen=scalar(@data);
		@data=sort{$a<=>$b} @data;
		if (int($arraylen/2)*2 == $arraylen) 
			{
				return ($data[$arraylen/2]+$data[$arraylen/2-1])/2;
			}else
				{
					return $data[int($arraylen/2)];
				}
	}

sub mean()
	{
		(my $data)=@_;
		my @data=split(/ /,$data);
		my $sum=0;
		foreach  (@data) 
			{
				$sum=$sum+$_;
			}
		return int(($sum/scalar(@data))*1000)/1000;
	}

sub ReplaceName()
	{
		(my $sp,my $file)=@_;
		my @content;
		my $line;
		my $i;
		my $target;
		open(F,$file);
		@content=<F>;
		close(F);
		for ($line=0;$line<@content;$line++) 
			{
				for ($i=0;$i<@$sp;$i++) 
					{
						$_=$content[$line];
						$target="sp".$i."sp";
						s/$target/$$sp[$i]/;
						$content[$line]=$_;
					}
			}
		open(R,">$file");
		print R @content;
		close(R);
	}

sub MP_Result_to_Table()
	{
		(my $MPresult, my $outputfile)=@_;
		my %hash;
		my $maxid=0;
		my $i;
		my @row;

		open(F,"$MPresult");
		$_=<F>;
		while (<F>) 
			{
				@row=split(/\t/,$_);
				if (exists($hash{$row[0]})) 
					{
						$hash{$row[0]}=$hash{$row[0]}."\t".$row[2];
					}else
						{
							$hash{$row[0]}=$row[2];
							if ($row[0]>$maxid) 
								{
									$maxid=$row[0];
								}
						}
			}
		close(F);

		open(R,">$outputfile");
		foreach $i (1..$maxid) 
			{
				print R $hash{$i}."\n";
			}
		close(R);
	}



sub FormatCluster()
	{
		(my $infile,my $genelist,my $spnum,my $cluster)=@_;
		my %hash;
		my %gene;
		my $key;
		my @row;
		my $sp;
		my $line;
		my $i=0;
		my $j=0;
		my @content;

### record gene in clusters
		open(F,"$infile");
		@content=<F>;
		close(F);
		chomp(@content);
		for ($line=0;$line<@content;$line++) 
			{
				@row=split(/\t/,$content[$line]);
				foreach $key (@row) 
					{
						$gene{$key}=1;
					}
			}
###retrieves gene which is not in clutsers

		open(F,"$genelist");
		while ($key=<F>) 
			{
				if ($key=~/^>/) 
					{
						chomp($key);
						$_=$key;
						s/^>//;
						$key=$_;
						if (!exists($gene{$key})) 
							{
								$content[$line]=$key;
								$line++;
							}
					}
			}
		close(F);

#### initialization @cluster
		@$cluster=qw();
		$j=0;

		foreach $line (@content)
			{
			if ($line ne "") 
				{
					%hash=();
					@row=split(/\t/,$line);
					foreach $key (@row) 
						{
							$sp=substr($key,0,index($key,"G"));
							$gene{$key}=1;
							if (exists($hash{$sp})) 
								{
									$hash{$sp}=$hash{$sp}.",".$key;
								}else
									{
										$hash{$sp}=$key;
									}
						}

					$i=0;
					@row=qw();
					
					foreach $i (0..($spnum-1)) 
						{
							if (exists($hash{"S$i"})) 
								{
									$row[$i]=$hash{"S$i"};
								}else
									{
										$row[$i]="-";
									}
						}
					$$cluster[$j++]=join("\t",@row);
				}
			}
	}


