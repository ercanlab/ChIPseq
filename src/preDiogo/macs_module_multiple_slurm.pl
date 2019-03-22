#!/usr/bin/perl -w
# This script takes in treat.bam and its matching input (control.bam) files for multiple samples
# performs MACS generating variable step wiggle files, converts it to fixed step start to end wiggle files
# performs normalization using wignorm generating variable step normalized wiggle files, converts it to fixed step start to end wiggle files
# usage perl macs_module.pl file containing treat.bam control.bam identifier
# change by Anna 04/19/12

use strict;
use Cwd;

my $mdir = getcwd;
my $dir = getcwd;

	chdir($dir);
	my $treatbam = shift();
	my $controlbam = shift();
	my $identifier = shift();
	print "Running MACS for $treatbam...\n";
	my $log3 = $mdir."/MACS_module_log_$identifier.txt"; 
	open (LOG3,">$log3");
	print LOG3 "Input Directory: $mdir\n";
	print LOG3 "Input treat.bam file: $treatbam\n";
	print LOG3 "Input control.bam file: $controlbam\n";
	print LOG3 "MACS output file identifier: $identifier\n\n\n";

	my $output_dir = $treatbam;
	$output_dir =~ s/.bam//;
	mkdir($output_dir);
	chdir($output_dir);

	unless (($treatbam=~/.bam/) && ($controlbam=~/.bam/)) {
		print LOG3 "Error: Please enter Input files in .bam format\n";
		exit(0);
	}	

	if ($treatbam =~ /gz/) {
		my $cmd = "gzip -d $treatbam";
		chomp($cmd);
		system($cmd);
		$treatbam =~ s/\.gz//g;
		print LOG3 "Unzipped file: $treatbam\n";	
	}
	if ($controlbam =~ /gz/) {
		my $cmd = "gzip -d $controlbam";
		chomp($cmd);
		system($cmd);
		$controlbam =~ s/\.gz//g;
		print LOG3 "Unzipped file: $controlbam\n";	
	}

	mkdir("MACS_pval1e-10") || print $!;
	mkdir("MACS_pval1e-5") || print $!;
	my $macs5folder = $dir."/".$output_dir."/MACS_pval1e-5";
	my $macs10folder = $dir."/".$output_dir."/MACS_pval1e-10";

	chdir($macs5folder)|| print $!;
	print LOG3 "In folder $macs5folder to run MACS at p-value 1e-5 and generate treat and control variable step wiggle files \n";
	print LOG3 "Running MACS1.4 issuing the command:";
	my $macs5_cmd = "macs -t $treatbam -c $controlbam -n MACS141_e5_$identifier -g ce -f BAM -w -S --space=1";
	print LOG3 "$macs5_cmd\n";
	# the following line is the command line running macs e-5 is default, thus is not indicated.
	my $macs5 = `macs -t $dir/$treatbam -c $dir/$controlbam -n MACS141_e5_$identifier -g ce -f BAM -w -S --space=1`;
	print LOG3 "Finished running MACS1.4 at p-value 1e-5\n\n";

	chdir($macs10folder)|| print $!;
	print LOG3 "In folder $macs10folder to run MACS at p-value 1e-10\n";
	print LOG3 "Running MACS1.4 issuing the command:";
	my $macs10_cmd = "macs -t $treatbam -c $controlbam -n MACS141_e10_$identifier -g ce -f BAM -p 1e-10";
	print LOG3 "$macs10_cmd\n";
	my $macs10 = `macs -t $dir/$treatbam -c $dir/$controlbam -n MACS141_e10_$identifier -g ce -f BAM -p 1e-10`;
	print LOG3 "Finished running MACS1.4 at p-value 1e-10\n\n\n";

	### wignorm
	chdir($mdir."/".$output_dir)|| print $!;
	print LOG3 "In folder $mdir/$output_dir to perform normalization using WIGNORM\n";
	print LOG3 "Unzipping the wiggle files from MACS1.4\n";
	my $treatwig = $macs5folder."/MACS141_e5_".$identifier."_MACS_wiggle/treat"."/MACS141_e5_".$identifier."_treat_afterfiting_all.wig.gz";
	my $cmd1 = "gzip -d $treatwig";
	chomp($cmd1);
	system($cmd1);
	$treatwig =~ s/\.gz//g;
	print LOG3 "\nUnzipped treat wiggle file: $treatwig";	
	
	#$controlwig=$mdir."/MACS_pval1e-5/MACS5_MACS_wiggle/control/MACS5_control_afterfiting_all.wig.gz";
	my $controlwig = $macs5folder."/MACS141_e5_".$identifier."_MACS_wiggle/control"."/MACS141_e5_".$identifier."_control_afterfiting_all.wig.gz";
	my $cmd2 = "gzip -d $controlwig";
	chomp($cmd2);
	system($cmd2);
	$controlwig =~ s/\.gz//g;
	print LOG3 "\nUnzipped control wiggle file: $controlwig\n";
	
	my $wignorm_cmd = "wignorm -t $treatwig -c $controlwig -n treat_control";
	print LOG3 "Running wignorm issuing the command: $wignorm_cmd\n";
	my $wignorm = `wignorm -t $treatwig -c $controlwig -n $identifier"_wignorm"`;
	print LOG3 "Finished running wignorm\n\n\n";
	# my $rm_scores_bed = `rm $mdir/$output_dir/treat_control_wignorm_peaks.bed`;
	
	### convert the variable step wiggle files to fixed step wiggle files
	my $wigfile_t = $macs5folder."/MACS141_e5_".$identifier."_MACS_wiggle/treat"."/MACS141_e5_".$identifier."_treat_afterfiting_all.wig";
	my $wigfile_c = $macs5folder."/MACS141_e5_".$identifier."_MACS_wiggle/control"."/MACS141_e5_".$identifier."_control_afterfiting_all.wig";
	my $wigfile = $mdir."/".$output_dir."/".$identifier."_wignorm_scores.wig";
	my @wfiles = ($wigfile_t,$wigfile_c,$wigfile);

	my %chr_len = ("CHROMOSOME_I" => 15072423, "CHROMOSOME_II" => 15279345, "CHROMOSOME_III" => 13783700, "CHROMOSOME_IV" => 17493793, "CHROMOSOME_MtDNA" => 13794, "CHROMOSOME_V" => 20924149, "CHROMOSOME_X" => 17718866, "plasmid+rexwithSNPs" => 7763);
	
	foreach my $f (@wfiles) {
		print LOG3 "Converting variable step wiggle file $f to fixed step wiggle file\n";

		my $end = 0;	
		my $cnt = 1;
		my $flag = 0;
		my $out = $f;
		$out =~ s/\.wig/\_FIXED\.wig/g;
		open (OUT,">$out");
		my $currchrom = "none";
		# print STDERR "\nProcessing wig file $f.\nPrinting output to $out.\n";
		open (WIG,$f);
		while (<WIG>) {
			chomp();
		
			if ($flag == 0) {
				if ($_ =~ 'track') {
					print OUT "$_ windowingFunction=mean smoothingWindow=8 viewLimits=-2:5\n";
				}
				if ($_ =~ 'variableStep') {
					my ($variable,$chrom,$span) = split(" ",$_);
					$chrom =~ s/chrom\=//g;
					$span =~ s/span\=//;
					if($currchrom ne $chrom) {
						# print STDERR "\tProcessing chrom $chrom\n";
						$currchrom = $chrom;
						$end = $chr_len{$currchrom};
						print OUT "fixedStep chrom=$chrom start=1 step=1 span=$span\n";
						$flag=1;
					}	
				}
			} else {
				my ($pos,$val) = split(" ",$_);
				if($pos ne 'variableStep') {
					if ($cnt <= $pos && $pos <= $end) {
						for (my $i = $cnt; $i < $pos; $i++) {
							print OUT " \n";
							$cnt++;  
						}
						if ($pos <= $end) {
							print OUT $val;
							print OUT "\n";
							$cnt++;
						}
					}
				} else {
					for (my $i = $cnt; $i <= $end; $i++) {
						print OUT " \n";  
					}   
					$cnt= 1;
					my ($variable,$chrom,$span) = split(" ", $_);
					$chrom =~ s/chrom\=//g;
					$span =~ s/span\=//;
					if ($currchrom ne $chrom) {
						# print STDERR "\tProcessing chrom $chrom\n";
						$currchrom = $chrom;
						$end = $chr_len{$currchrom};
						print OUT "fixedStep chrom=$chrom start=1 step=1 span=$span\n";
						$flag = 1;		
					}
				}
			}
		} # end of while
	  
		for (my $i = $cnt; $i <= $end; $i++) {
			print OUT " \n";
		}     

		print LOG3 "Check fixed step normalized wiggle file: $out\n";
	
		my $cmd = `gzip $f`;
		print LOG3 "Gzipped file: $f\n";
	
		$cmd = `gzip $out`;
		print LOG3 "Gzipped file: $out\n\n\n";
	}

	### Check peak bed file for peak regions smaller than chromosome start or larger than chromosome length and move these to 0/chromsome length
	my @bedfiles = ($macs5folder."/MACS141_e5_".$identifier."_peaks.bed", $macs10folder."/MACS141_e10_".$identifier."_peaks.bed", $mdir."/".$output_dir."/".$identifier."_wignorm_peaks.bed");

	for my $bedfile_cur (@bedfiles) {
		my $bedfile_new = $bedfile_cur;
		$bedfile_new =~ s/.bed/_corrected.bed/;
		
		print "Total number of peaks in $bedfile_cur:\n";
                system("wc -l $bedfile_cur");
		print LOG3 "Checking peak bed file $bedfile_cur for annotations outside the chromosome regions.\n";
	
		open (INPUT, "<", $bedfile_cur) || print "$!\n";
		open (OUTPUT, ">", $bedfile_new) || print "$!\n";
		while (<INPUT>) {
			chomp();
			my ($chr, $start, $end, $peak, $value) = split("\t", $_);
			$start = 0 if ($start < 0);
			$end = $chr_len{$chr} if ($end > $chr_len{$chr});
			if (defined($peak) && defined($value)) {
            			print OUTPUT "$chr\t$start\t$end\t$peak\t$value\n";
        		} else {
          			print OUTPUT "$chr\t$start\t$end\n";
        		}
		}	
	
	close (INPUT);
		close (OUTPUT);
	
		system("rm $bedfile_cur");
		rename($bedfile_new, $bedfile_cur);
	}


