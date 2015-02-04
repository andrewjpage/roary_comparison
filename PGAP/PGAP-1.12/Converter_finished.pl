#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

my %opt;
getopts('S:I:O:',\%opt);

my @usage=qq(
Usage:   perl Converter_finished.pl [options]

Options: 

  -S String    Input the strains nickname, 
               If 2 or more, join them with '+', 
               For example: CT18+NC_011834+SPA
  -I String    Input file directory
  -O String    Output file directory
);

if (!scalar(keys %opt)) 
	{
		print join("\n",@usage)."\n";
		exit;
	}

my @sp;
if (exists($opt{"S"})) 
	{
		@sp=split(/\+/,$opt{"S"});
	}else
		{
			print "-S could not be empty!";
			print join("\n",@usage)."\n";
			exit;
		}

my $output;
if (exists($opt{"O"})) 
	{
		$output=$opt{"O"};
	}else
		{
			print "-O could not be empty!";
			print join("\n",@usage)."\n";
			exit;
		}

my $input;
if (exists($opt{"I"})) 
	{
		$input=$opt{"I"};
	}else
		{
			print "-I could not be empty!";
			print join("\n",@usage)."\n";
			exit;
		}

my $sp;
my $line;
my @row;
my @tmp;
my %hash;
my $flag;
my $file;


if ((-e $output) and ((-d $output))) 
	{
	}else
	{
		mkdir($output);
	}

if ($input!~/\/$/) 
	{
		$input=$input."/";
	}


if ($output!~/\/$/) 
	{
		$output=$output."/";
	}


foreach $sp (@sp) 
	{
		%hash=();
		$file=$input.$sp.".faa";
		open(F,$file) or die "could not open $file";
		open(R,">$output$sp.pep");
		while ($line=<F>) 
			{
				if ($line=~/^>/) 
					{
						@row=split(/\|/,$line);
						print R ">$row[1]\n";
					}else
						{
							print R $line;
						}
			}
		close(F);
		close(R);

		$file=$input.$sp.".ptt";
		open(F,"$file") or die "could not open $file";
		open(R,">$output$sp.function");
		$_=<F>;
		$_=<F>;
		$_=<F>;
		while ($line=<F>) 
			{
				chomp($line);
				@row=split(/\t/,$line);
				print R $row[3]."\t".$row[7]."\t".$row[8]."\n";
				@tmp=split(/\.\./,$row[0]);
				if ($row[1] eq "+") 
					{
						$hash{$tmp[0]."-".$tmp[@tmp-1]}=$row[3];
					}else
						{
							$hash{"c".$tmp[@tmp-1]."-".$tmp[0]}=$row[3];
						}
			}
		close(R);
		close(F);

		$file=$input.$sp.".ffn";
		open(F,"$file") or die "could not open $file";;
		open(R,">$output/$sp.nuc");
		while ($line=<F>) 
			{
				if ($line=~/^>/) 
					{
						my $key=&getKey($line);
						if (exists($hash{$key})) 
							{
								$flag=1;
								print R ">$hash{$key}\n";
							}else
								{
									$flag=0;
								}
					}else
						{
							if ($flag) 
								{
									print R $line;
								}
						}
			}
		close(R);
		close(F);
	}

sub getKey()
{
	(my $line)=@_;
	my @tmp;
	my $strand;
	chomp($line);
	@tmp=split(/ /,$line);
	@tmp=split(/\:/,$tmp[0]);

	if($tmp[@tmp-1]=~/c/)
	{
		$strand="-";
	}else
	{
		$strand="+";
	}
	$_=$tmp[@tmp-1];
	s/c//g;
	s/ //g;
	@tmp=split(/\,|-/,$_);
	@tmp=sort{$a<=>$b} @tmp;
	if($strand eq "-")
	{
		return "c".$tmp[@tmp-1]."-".$tmp[0];
	}else
	{
		return $tmp[0]."-".$tmp[@tmp-1];
	}
}
