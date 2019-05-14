#!/usr/bin/perl -w

################################################
#################################################
## Script to generate bigwig and bigbed files
## Run script on mercer. module load kent/328
## Date: 01/09/2017
## Author: Sarah Albritton
## Email: sea283@nyu.edu or sarahea283@gmail.com
#################################################
#################################################



# This script takes in a wiggle file and converts it to bigwig
# Usage:
# perl getBW.pl files.txt


my %chr_len = ("chrI" => 15072423, "chrII" => 15279345, "chrIII" => 13783700, "chrIV" => 17493793, "chrM" => 13795, "chrV" => 20924149, "chrX" => 17718866); # ce10

use strict;
use Cwd;

#my $mdir = getcwd;
my $dir = getcwd;

#chdir($dir);

my $file = shift();
my $priority = 1;

    if ($file =~ /.*gz/) {
        system("gunzip $file");
        $file =~ s/\.gz//g;
    }

    if ($file =~ /wig$/) {
        loadWig("$file");
    } elsif ($file =~ /bed$/) {
        loadBed("$file");
    }


#### Function declarations ####

sub loadWig {

    my $temp = $file;
    $temp =~ s/.wig/_temp.wig/;
    my $desc = $file;
    $desc =~ s/.wig//;
    my $chrom;
    my $count = 1;

# 1. check for empty lines -> change to zero for uploading
    open (INPUT, "<", $file) || die "Couldn't open file $file: $!.\n";
    open (OUTPUT, ">", $temp) || die "Couldn't open file $temp: $!.\n";
    while (<INPUT>) {
        chomp();
        if (/liftover/) {
            next;
        } elsif (/^track/) {
                print OUTPUT "track type=wiggle_0 name=\"$desc\" description=\"$desc\"\n";
		} elsif (/^fixedStep/) {
                if (/fixedStep chrom=([a-zA-Z0-9_]*?) start=1 step=1 span=1/) {
                        $chrom = $1;
                        #print "Processing $chrom...\n";
                        if ($chrom =~ /CHROMOSOME/) {
                            $chrom =~ s/CHROMOSOME_/chr/;
                        }
                        if ($chrom =~ /MtDNA/) {
                            $chrom =~ s/MtDNA/M/;
                        }
                        if ($chrom =~ /^[IVX]+/) {
                            $chrom = "chr".$chrom;
                        }
                        print OUTPUT "variableStep chrom=$chrom\n";
                        $count = 1;
                } else {
                        print "Problem: $_\n";
                        exit(1);
                }
		} elsif (/^ $/) {
                $count++;
        } else {
                if ($count > $chr_len{$chrom}) {
                    print "Problem for count: $count in line $. for $chrom.\n";
                }
                print OUTPUT "$count $_\n";
                $count++;
        }
    } # end of while

    close INPUT;
    close OUTPUT;

# 2. convert wig to bigwig
    my $bwfile = $file;
    $bwfile =~ s/.wig/.bw/;
    system("wigToBigWig  $temp /scratch/cgsb/ercan/annot/ce10_chromInfo.txt $bwfile");

# 3. delete wig file
    system("rm $temp");
    #system("rm $wig");
}




sub loadBed {
    my $file = shift();
    # 1. convert to USCS chromosome format and sort bed file
    my $temp = $file;
    $temp =~ s/.bed/_temp.bed/;
    open (INPUT, "<", $file) || die "Couldn't open file $file: $!.\n";
    open (OUTPUT, ">", $temp) || die "Couldn't open file $temp: $!.\n";
    while (<INPUT>) {
        chomp();
	if ($_ =~ /liftover/ || $_ =~ /#/) {
		next;
	}
        if ($_ =~ /CHROMOSOME/) {
            s/CHROMOSOME_/chr/;
            s/MtDNA/M/ if (/MtDNA/);
        }
        print OUTPUT $_."\n";
    }
    close INPUT;
    close OUTPUT;
	system("sort -k1,1 -k2,2n $temp > $file");
    # 2. convert bed to bigbed
    my $bbfile = $file;
    $bbfile =~ s/.bed/.bb/;
        system("bedToBigBed -as=/scratch/cgsb/ercan/scripts/trackhubs/narrowPeak.as -type=bed6+4 $file /scratch/cgsb/ercan/annot/ce10_chromInfo.txt $bbfile");
    system("rm $temp");
}
