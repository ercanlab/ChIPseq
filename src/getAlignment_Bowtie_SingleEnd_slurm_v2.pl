# this script aligns single-end reads in fastq format using bowtie
# usage
# go to the directory where you have the required files
# perl getAlignment_Bowtie_SingleEnd.pl
# updated by Anna 11/30/11

use Cwd;
use strict;
use warnings;
my $dir = getcwd;

my $file = shift();

	my @temp = split(/\./, $file);
	my $id = $temp[0];
	$id =~ s/.fastq//;
	my $log2 = $dir."/Read_Alignment_module_log_".$id.".txt";
	open(LOG2,">$log2");
	print LOG2 "Input Directory: $dir\n";
	print LOG2 "Input Sequenced Reads .fastq file: $file\n\n\n";

	my $indexfile = "/scratch/cgsb/ercan/annot/forBowtie/c_elegans.WS220";
	my $samfile = $file;
	$samfile =~ s/.fastq/.sam/;
	my $bamfile = $samfile;
	$bamfile =~ s/.sam/.bam/;
	my $bresult = $id."_bresult.txt";
	print LOG2 "Index file is: $indexfile\n";
	print LOG2 "Output file will be: $bamfile\n";
	print LOG2 "Running BOWTIE issuing the command:";
	my $bowtie_cmd = "bowtie2 -p 8 -x $indexfile -U $file -S $samfile";
	print LOG2 "$bowtie_cmd\n";
	my $start = time();
	my $bowtie = `bowtie2 -p 8 -x $indexfile -U $file -S $samfile 2>> $bresult`;
	my $time = int((time() - $start)/60);
	print LOG2 "Finished running BOWTIE ($time min)\n\n\n";
	print LOG2 "Converting $samfile to $bamfile using the command:\n";
	my $sam_to_bam_cmd = "samtools view -bS $samfile > $bamfile";
	print LOG2 "$sam_to_bam_cmd\n";
	my $sam_to_bam = `samtools view -bS $samfile > $bamfile`;
	print LOG2 "Converted .sam file to .bam file\n\n\n";

	my $r = `rm $samfile`;
	print LOG2 "Removed $samfile\n";
	print LOG2 "END OF RUN\n\n\n";

	print LOG2 "BOWTIE RESULT\n";
	my $append_bresult = `cat $bresult >>$log2`;
	my $rm = `rm $bresult`;
