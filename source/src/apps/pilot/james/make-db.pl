#!/usr/bin/perl -w

use strict;

my $format_string = "%3d";

foreach my $file (@ARGV) {
	open FILE, "<$file" or die $!;
	my @lines = <FILE>;
	close FILE or die $!;
	
	my $id = $file;
	my $sequence = '';
	my @profile_lines;
	foreach my $line (@lines) {
		chomp $line;
		if ( $line !~ /^\s*\d+\s+\w/ ) { next; }

		$line =~ s/^\s*//g;
		my @e = split /\s+/, $line;
		my $resi = shift @e;
		my $res  = shift @e;

		$sequence .= $res;
		my @log_probs;
		for ( 1 .. 20 ) {
			my $log_prob = shift @e;
			push @log_probs, sprintf( $format_string, $log_prob );
		}

		my $profile_line = join ' ', @log_probs;
		push @profile_lines, $profile_line;
	}

	print "ENTRY: $id\n";
	print "SEQUENCE: $sequence\n";
	print join "\n", @profile_lines;
	print "\n--\n";
}
