#!/usr/bin/env perl

use Moose;
use Getopt::Long;
use Parallel::ForkManager;
use Bio::Tools::GFF;
use File::Slurp;
use Data::Dumper;
use File::Basename;
use Cwd;
use Bio::PanGenome::ExtractProteomeFromGFF;
my ( $help,$parallel_processes );

my $start_time = time();

GetOptions( 'h|help'  => \$help,'parallel_processes=i' => \$parallel_processes,);

if($help)
{
  print "Script to run PGAP from a GFF3 file\nrun_pgap.pl *.gff";
}
$parallel_processes ||=1;
my @gff_files = @ARGV;

print "Extracting proteomes, nuc and annotation\n";
my $pm = new Parallel::ForkManager( $parallel_processes );
for my $file (@gff_files)
{
  $pm->start and next;    # fork here
  # Extract protein sequences from the GFF files of nucleotides
  
  my($base_filename, $dirs, $suffix) = fileparse($file, qr/\.[^.]*/);
  my $original_base_filename = $base_filename;
  $base_filename =~ s![^a-zA-Z0-9]!!g;
  
  my $cmd = "extract_proteome_from_gff -o faa $file";
  system($cmd);
  system("mv ${original_base_filename}${suffix}.faa ${base_filename}.pep");
  
  # Extract annotation for each sequence
  open(my $gene_att_fh, '>', $base_filename.'.function');
  my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => 3);
  my $feature;
  while($feature = $gffio->next_feature()) {
      next unless ( $feature->has_tag('product') );
      next unless ( $feature->has_tag('ID') );
      my ($gene_id,$product, @junk);
      ($gene_id, @junk) = $feature->get_tag_values('ID');
      ($product, @junk) = $feature->get_tag_values('product');
      print {$gene_att_fh} join("\t",($gene_id, '-', $product))."\n";
  }
  $gffio->close();
  close($gene_att_fh);
  
   my $obj = Bio::PanGenome::ExtractProteomeFromGFF->new(
       gff_file        => $file,
     );
   $obj->_extract_nucleotide_regions();
   system("mv ".$obj->_extracted_nucleotide_fasta_file_from_bed_filename." ${base_filename}.nuc");
   

  $pm->finish;
}
$pm->wait_all_children;

my @strain_names ;
for my $file (@gff_files)
{
  my($base_filename, $dirs, $suffix) = fileparse($file, qr/\.[^.]*/);
  $base_filename =~ s![^a-zA-Z0-9]!!g;
  push(@strain_names, $base_filename);
}
my $strain_name_input_str = join('+',@strain_names);
my $cwdir = getcwd;
my $cmd = "perl PGAP.pl -strains $strain_name_input_str -input $cwdir -output output_directory --cluster --pangenome --variation --evolution -thread $parallel_processes --function -method MP";
# Also GP method
print $cmd."\n";
system($cmd);
