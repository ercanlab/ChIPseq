#! /usr/bin/perl -w
use strict;

#################################################################
#################################################################
# Script to load tracks into UCSC track hub (reading files from file)
# Run script within directory with your files to upload
# Date: 01/09/17
# Edited by: Sarah Albritton       (original author: Anna-Lena Kranz )
# Email : sea283@nyu.edu or sarahea283@gmail.com
#################################################################
#################################################################

my $files = shift();
my $priority = 1;
my $trackDb = "trackDb.txt";


my %chr_len = ("chrI" => 15072423, "chrII" => 15279345, "chrIII" => 13783700, "chrIV" => 17493793, "chrM" => 13794, "chrV" => 20924149, "chrX" => 17718866); # ce10
#my %chr_len = ("chrI" => 15455979, "chrI_random" => 43675, "chrII" => 16627154, , "chrIII" => 14578851, "chrIII_random" => 105344, "chrIV" => 17485439, "chrIV_random" => 25117, "chrV" => 19495157, "chrV_random" => 85292, "chrX" => 21540570, "chrX_random" => 28673, "chrun" => 2948414); # briggsae WS228

# check if trackDb.txt already exists. If so, delete
#if (-e "trackDb.txt") {
#    system("rm trackDb.txt");
#}

open (INFILE, "<", $files) || die "Couldn't open file $files: $!.\n";
while (<INFILE>) {
    chomp();
    my $file = $_;
    print "Loading $file...\n";
    # 2. check file type and write to track.txt
    if ($file =~ /bw$/) {
        loadBigWig("$file");
    } elsif ($file =~ /bb$/) {
        loadBigBed("$file");
    } elsif ($file =~ /bedgraph$/) {
        loadBedGraph("$file");
    }
}
close INFILE;

#### Function declarations ###


sub loadBedGraph {
    my $file = shift();
    # 1. check for empty lines -> change to zero for uploading
    my $temp = $file;
    $temp =~ s/.bedgraph/_temp.bg/;
    open (INPUT, "<", $file) || die "Couldn't open file $file: $!.\n";
    open (OUTPUT, ">", $temp) || die "Couldn't open file $temp: $!.\n";
    while (<INPUT>) {
        chomp();
        if ($_ eq " ") {
            print OUTPUT "0\n";
        } elsif ($_ =~ /CHROMOSOME/) {
            s/CHROMOSOME_/chr/;
            s/MtDNA/M/ if (/MtDNA/);
            print OUTPUT $_."\n";
        } else {
            print OUTPUT $_."\n";
        }
    }
    close INPUT;
    close OUTPUT;
    # 2. sort bedgraph according to bedGraphToBigWig instructions
    my $sorted = $temp;
    $sorted =~ s/_temp.bg/_sorted.bg/;
    system("sort -k1,1 -k2,2n $temp > $sorted");
    # 3. convert bedgraph to bigwig
    my $bwfile = $file;
    $bwfile =~ s/.bedgraph/.bw/;
    system("/Users/Sarah/Desktop/Computing_Resources/CommandLineTool/bedGraphToBigWig $sorted /Users/Sarah/Desktop/Data_Analysis/annotation/ce10_chromInfo.txt  $bwfile");
    # 4. delete wig file
    system("rm $file");
    system("rm $temp");
    system("rm $sorted");
    # 5. get name and labels
    my ($name, $shortLabel, $longLabel);
    if ($bwfile =~ /MACS141_e5_(.*?)_average.bw/) {
        $name = $1;
        $longLabel = $file;
        $shortLabel = $1;
    }  elsif ($bwfile =~ /(.*?).bw/) {
        $name = $1;
        $longLabel = $file;
        $shortLabel = $1;
    } else {
        $name = "SEA54_TY1916";
        $longLabel = $file;restart
        $shortLabel = "SEA54_TY1916_normed_subt";
        #print "Problem identifying name, shortLabel and longLabel!\n";
        #exit(1);
    }
    # 6. add to local trackDb.txt
    open (OUTPUT, ">>", $trackDb) || die "Couldn't open file $trackDb: $!.\n";
    print OUTPUT "track $name\n";
    print OUTPUT "type bigWig -2.0 5.0\n";
    print OUTPUT "visibility full\n";
    print OUTPUT "group user\n";
    print OUTPUT "priority $priority\n";
    print OUTPUT "shortLabel $shortLabel\n";
    print OUTPUT "longLabel $longLabel\n";
    print OUTPUT "bigDataUrl $bwfile\n";
    print OUTPUT "windowingFunction mean\n";
    print OUTPUT "smoothingWindow 8\n\n";
    close OUTPUT;
    $priority++;
}

sub loadBigWig {
    my $file = shift();
    # 1. get name and labels
    my ($name, $shortLabel, $longLabel);
    if ($file =~ /MACS141_e5_(.*?)_treat_afterfiting.*.bw/) {
        $name = $1."_wiggle";
        $longLabel = $file;
        $shortLabel = $1."_wiggle";
    } elsif ($file =~ /MACS141_e5_(.*?)_average.bw/) {
        $name = $1."_wiggle";
        $longLabel = $file;
        $shortLabel = $1."_wiggle";
    } elsif ($file =~ /..\/..\/data\/(.*?).bw/) {
        $name = $1."_wiggle";
        $longLabel = $1."_wiggle";
        $shortLabel = $1."_wiggle";
    } elsif ($file =~ /(.*?).bw/) {
        $name = $1."_wiggle";
        $longLabel = $file;
        $shortLabel = $1."_wiggle";
    } else {
        print "Problem identifying name, shortLabel and longLabel!\n";
        exit(1);
    }
    # 2. add to local trackDb.txt
    open (OUTPUT, ">>", $trackDb) || die "Couldn't open file $trackDb: $!.\n";
    print OUTPUT "track $name\n";
    print OUTPUT "type bigWig\n";
    print OUTPUT "visibility full\n";
    print OUTPUT "group user\n";
    print OUTPUT "priority $priority\n";
    print OUTPUT "shortLabel $shortLabel\n";
    print OUTPUT "longLabel $longLabel\n";
    print OUTPUT "bigDataUrl $file\n";
    print OUTPUT "windowingFunction mean\n";
    print OUTPUT "smoothingWindow 8\n";
    print OUTPUT "autoScale on\n\n";
    close OUTPUT;
    $priority++;
}



sub loadBigBed {
    my $file = shift();
    # 1. get labels and name
    my ($name, $shortLabel, $longLabel);
    if ($file =~ /..\/..\/data\/(.*?).bb/) {
        $name = $1."_peaks";
        $longLabel = $file;
        $shortLabel = $1."_peaks";
    } elsif ($file =~ /MACS141_e10_(.*?)_peaks/) {
        $name = $1."_peaks";
        $longLabel = $file;
        $shortLabel = $1."_peaks";
    } elsif ($file =~ /(.*?).bb/) {
        $name = $1."_peaks";
        $longLabel = $file;
        $shortLabel = $1."_peaks";
    } else {
        print "Problem identifying name, shortLabel and longLabel!\n";
        exit(1);
    }
    # 2. add to local trackDb.txt
    open (OUTPUT, ">>", "$trackDb") || die "Couldn't open file $trackDb: $!.\n";
    print OUTPUT "track $name\n";
    print OUTPUT "type bigBed 4 .\n";
    print OUTPUT "visibility dense\n";
    print OUTPUT "group user\n";
    print OUTPUT "priority $priority\n";
    print OUTPUT "shortLabel $shortLabel\n";
    print OUTPUT "longLabel $longLabel\n";
    print OUTPUT "bigDataUrl $file\n\n";
    close OUTPUT;
    $priority++;
}

sub NummernSort {
    if ($a < $b) {
        return -1;
    } elsif($a == $b) {
        return 0;
    } else {
        return 1;
    }
}
