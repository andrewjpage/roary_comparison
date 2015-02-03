#!/usr/bin/env perl

use Moose;
use Getopt::Long;
use Parallel::ForkManager;
use Bio::Tools::GFF;
my ( $help,$parallel_processes );

GetOptions( 'h|help'  => \$help,'parallel_processes=i' => \$parallel_processes,);

if($help)
{
  print "Script to run panOCT from a GFF3 file\nrun_pan_oct.pl *.gff";
}
$parallel_processes ||=8;
my @gff_files = @ARGV;

my $pm = new Parallel::ForkManager( $parallel_processes );
for my $file (@gff_files)
{
  $pm->start and next;    # fork here
  # Extract protein sequences from the GFF files of nucleotides
  my $cmd = "extract_proteome_from_gff -o faa $file";
  system($cmd);
  
  # Extract annotation for each sequence
  open(my $gene_att_fh, '>', $file.'.gene_att');
  my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => 3);
  my $feature;
  while($feature = $gffio->next_feature()) {
      next unless ( $feature->has_tag('product') );
      next unless ( $feature->has_tag('ID') );
      my ($gene_id,$product, @junk)
      ($gene_id, @junk) = $feature->get_tag_values('ID');
      ($product, @junk) = $feature->get_tag_values('product');
      my $gene_start = $feature->start;
      my $gene_end   = $feature->end;
      my $chromosome = $feature->seq_id
      print {$gene_att_fh} join("\t",($chromosome, $gene_id, $gene_start,$gene_end, $product, $file))."\n";
  }
  $gffio->close();
  close($gene_att_fh);

  $pm->finish;
}
$pm->wait_all_children;

# Create tags file - just use the file names
open(my $tags_fh, '>', 'run_tags.txt');
for my $file (@gff_files)
{
  print {$tags_fh} $file."\n";
}
close($tags_fh);

# Merge  files together
my $faa_files = '';
my $att_files = '';
for my $file (@gff_files)
{
  $faa_files .= $file.".faa ";
  $att_files .= $file.".gene_att ";
}
system("cat $faa_files > run.pep");
system("rm $faa_files");
system("cat $att_files > run.gene_att");
system("rm $att_files");

# Run parallel blast
system("parallel_all_against_all_blastp -p $parallel_processes -j Parallel -o run_blast.txt run.pep");


system("perl PanOCT.pl -t run_blast.txt -f run_tags.txt -g run.gene_att -P run.pep ");
