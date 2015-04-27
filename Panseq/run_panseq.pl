

# take in gff files
strip out fasta file
change the header of each sequence
put into query directory
create & update settings file


sed -n '/##FASTA/,//p' 10209_5#1.gff  | grep -v '##FASTA' > 10209_5#1.fa
