#!/usr/bin/perl -w
use strict;

my $input = $ARGV[0];
my $files = $ARGV[1];
my $outfile = $input;
$outfile =~ s/.bed/_final.bed/;

my @files = split(",", $files);

my %cutoff = (1 => 1, 2 => 2, 3 => 2, 4 => 3, 5 => 3, 6 => 4, 7 => 4);
my %peaks;
my %results;

open (INPUT, "<", $input) || die "Couldn't open file $input: $!.\n";
while (<INPUT>) {
    chomp();
    my @line = split("\t", $_);
    unless ($results{$line[3]}) {
        $results{$line[3]} = $_;
    } else {
        print "Problem with peak $line[3]!\n";
    }
}
close INPUT;

my $output = sprintf("%s_temp.txt", $outfile);
foreach my $file (@files) {
    print "Comparing $input to $file...\n";
    system("intersectBed -wa -u -a $input -b $file > $output");
    sleep(10);
    open(INPUT, "<", $output) || die "Couldn't open intersect file: $!.\n";
    while (<INPUT>) {
        chomp();
        my @line = split("\t", $_);
        if ($peaks{$line[3]}) {
            $peaks{$line[3]}++;
        } else {
		$peaks{$line[3]} = 1;
        }
    }
    close INPUT;
}

system("rm $output");

my $n_files = scalar(@files);
my $cutoff = $cutoff{$n_files};
print "Number of files: $n_files\n";
print "Cutoff: $cutoff\n";
my $count = 0;
open(OUTPUT, ">", $outfile) || die "Couldn't open file $outfile: $!.\n";
foreach my $peak (keys %peaks) {
    if ($peaks{$peak} >= $cutoff) {
        print OUTPUT $results{$peak}."\n";
        $count++;
    }
}
close OUTPUT;

print $count." overlapping peaks found in replicates.\n";
