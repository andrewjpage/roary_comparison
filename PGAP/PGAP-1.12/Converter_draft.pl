#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

my %opt;
getopts('N:I:O:',\%opt);

my @usage=qq(
Usage:   perl Converter_draft.pl [options]

Options: 

  -N String    Input the strain nickname
  -I String    Input file directory
  -O String    Output file directory
);

if (!scalar(keys %opt)) 
	{
		print join("\n",@usage)."\n";
		exit;
	}

my $prefix;
if (exists($opt{"N"})) 
	{
		$prefix=$opt{"N"}
	}else
		{
			print "-N could not be empty!";
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
my $list;
my @list;
my $pttlost;
my $gi;

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

opendir(DIR,"$input") || die "The input directory ( $input ) is not exists!\n";
@list=grep(/faa$/,readdir(DIR));
closedir(DIR);
$_=join("\t",@list);
s/.faa//g;
@list=split(/\t/,$_);

open(PEP,">$output$prefix.pep");
open(NUC,">$output$prefix.nuc");
open(FUN,">$output$prefix.function");

foreach $list (@list) 
	{
		%hash=();
		if (!(-e $input.$list.".faa")) 
			{
				print $input.$list.".faa is not exists!\n$list.faa, $list.ffn and $list.ptt are skipped!\n";
				next;
			}

		if (!(-e $input.$list.".ffn")) 
			{
				print $input.$list.".ffn is not exists!\n$list.faa, $list.ffn and $list.ptt are skipped!\n";
				next;
			}

		if (!(-e $input.$list.".ptt")) 
			{
				open(F,$input.$list.".faa");
				@_=<F>;
				close(F);
				@_=grep(/^>/,@_);
				if (scalar(@_)>1) 
					{
					print $input.$list.".ptt is not exists!\n";
					print "There are more than 1 sequence in $list.faa and $list.ffn, So $list.faa, $list.ffn and $list.ptt are skipped!\n";
					next;
					}
				$pttlost=1;
			}else
				{
					$pttlost=0;
				}

		$file=$input.$list.".faa";
		open(F,$file) or die "could not open $file";
		while ($line=<F>) 
			{
				if ($line=~/^>/) 
					{
						@row=split(/\|/,$line);
						print PEP ">$row[1]\n";
						if ($pttlost ==1) 
							{
								$gi=$row[1];
							}
					}else
						{
							print PEP $line;
						}
			}
		close(F);

		if ($pttlost ==1) 
			{
				print FUN "$gi\t-\thypothetical protein\n";
			}else
				{
				$file=$input.$list.".ptt";
				open(F,"$file") or die "could not open $file";
				$_=<F>;
				$_=<F>;
				$_=<F>;
				while ($line=<F>) 
					{
						chomp($line);
						@row=split(/\t/,$line);
						print FUN $row[3]."\t".$row[7]."\t".$row[8]."\n";
						@tmp=split(/\.\./,$row[0]);
						if ($row[1] eq "+") 
							{
								$hash{$tmp[0]."-".$tmp[@tmp-1]}=$row[3];
							}else
								{
									$hash{"c".$tmp[@tmp-1]."-".$tmp[0]}=$row[3];
								}
					}
				close(F);
				}



		$file=$input.$list.".ffn";
		open(F,"$file") or die "could not open $file";;
		while ($line=<F>) 
			{
				if ($line=~/^>/) 
					{
						if ($pttlost==1) 
							{
								print NUC ">$gi\n";
								$flag=1;
							}else
								{
									my $key=&getKey($line);
									if (exists($hash{$key})) 
										{
											$flag=1;
											print NUC ">$hash{$key}\n";
										}else
											{
												$flag=0;
											}
								}
					}else
						{
							if ($flag) 
								{
									print NUC $line;
								}
						}
			}
		close(F);
	}

close(PEP);
close(NUC);
close(FUN);

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