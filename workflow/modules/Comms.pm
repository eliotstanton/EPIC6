package Comms;

use strict;
use warnings;
use Term::ANSIColor qw(:constants);

# ---------------------------------------------------------------------------- #

# Module name:	Comms
# Created by:	Eliot Stanton (estanton@wisc.edu)
# Created on:	08 September, 2023
# Modified:	09 October, 2023
# Description:	Subroutine for printing various messages to stdout and to log
#		files.

# ---------------------------------------------------------------------------- #

# Subroutines:
# - Comms::config_missing	Warns user that config file is missing data
# - Comms::container_info	Prints information about a given container
# - Comms::date_time		Returns the current date and time
# - Comms::dir_empty:		Warns user about empty directory
# - Comms::dir_missing:		Warns user about a missing directory
# - Comms::exit_code:		Warns user about non-zero exit code
# - Comms::file_missing:	Warns user that a file is missing
# - Comms::file_present:	Reports file as already present
# - Comms::file_structure	Prints analysis file structure
# - Comms::help:		Prints help file
# - Comms::num_contigs:		Prints the number of contigs in file
# - Comms::pair_missing:	Warns user about missing sequencing file mate
# - Comms::print_command	Prints command to user
# - Comms::start_finish		Prints messages when starting and finishing jobs
# - Comms::terminate_code	Prints message when job doesn't exit cleanly

# ---------------------------------------------------------------------------- #

# Subroutine:	Comms::config_missing
# Description:	Warns user that there is an issue with the config file missing
#		data.

# --------------------------------------

sub config_missing {

	# Arguments:
	my ( $var_string ) = @_;

	# ------------------------------

	print RED,"    ERROR: $var_string details missing from config file!\n", RESET;

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Comms::container_info
# Description:	Prints information about container

# --------------------------------------

sub container_info {

	# Arguments:
	my ( $hash_config, $var_command, $file_log ) = @_;

	# Data structures:
	my $hash_software	= $hash_config->{software};

	# Variables:
	my $var_tag		= $hash_software->{$var_command}->{tag};
	my $var_repo		= $hash_software->{$var_command}->{repository};
	my $var_cmd		= $hash_software->{$var_command}->{command};
	my $file_def		= $hash_software->{$var_command}->{def_file};
	my $var_string;

	# ------------------------------

	# Add container information to string:
	$var_string = "docker://$var_repo/$var_cmd:$var_tag" unless $file_def;
	$var_string = "$file_def -> $var_cmd:$var_tag" if $file_def;

	# Print to command-line:
	print WHITE, "    $var_string\n";

	# Print to $file_log if defined:
	if ( $file_log ) {

		Edit::append_string ( $file_log, "$var_string\n" );

	}

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Comms::date_time
# Description:	Returns current date and time.

# --------------------------------------

sub date_time {

	# Variables:
	my $var_date;
	my $var_time;

	# ------------------------------

	# Return local date and time:
	my ( $var_sec, $var_min, $var_hour, $var_mday, $var_mon,
	$var_year) = localtime();

	# Move date up by one:
	$var_mon += 1;

	# Add a preceding zero to $variables as needed:
	$var_mon = 0 . $var_mon if $var_mon < 10;
	$var_mday = 0 . $var_mday if $var_mday < 10;
	$var_hour = 0 . $var_hour if $var_hour < 10;
	$var_min = 0 . $var_min if $var_min < 10;
	$var_sec = 0 . $var_sec if $var_sec < 10;

	# Add 1900 to year:
	$var_year += 1900;

	# Format date and time:
	$var_date	= "$var_year-$var_mon-$var_mday";
	$var_time	= "$var_hour:$var_min:$var_sec";

	# ------------------------------

	# Return $var_date and $var_time and end subroutine:
	return ($var_date, $var_time);

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Comms::dir_empty
# Description:	Warns user that directory is empty

# --------------------------------------

sub dir_empty {

	# Arguments:
	my ( $var_string ) = @_;

	# ------------------------------

	print RED,"    ERROR: Directory $var_string is empty!\n", RESET;

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Comms::dir_missing
# Description:	Warns user that a directory is missing.

# --------------------------------------

sub dir_missing {

	# Arguments:
	my ( $var_string ) = @_;

	# ------------------------------

	print RED,"    ERROR: Directory $var_string is missing!\n", RESET;

	# ------------------------------

	# End subroutine:
	return;

}


# ---------------------------------------------------------------------------- #

# Subroutine:	Error::exit_code
# Description:	Prints message warning user that process did not exit cleanly.

# --------------------------------------

sub exit_code {

	# Arguments:
	my ( $file_in ) = @_;

	# -----------------------------

	print RED, "    ERROR: See $file_in for more details.\n", RESET;

	# -----------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Comms::file_missing
# Description:	Warns user that a file is missing.

# --------------------------------------

sub file_missing {

	# Arguments:
	my ( $var_string ) = @_;

	# ------------------------------

	print RED,"    ERROR: $var_string file is missing!\n", RESET;
	
	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Comms::file_present
# Description:	Prints message to user that results are already present.

# --------------------------------------

sub file_present {

	# Arguments:
	my ( $file_in ) = @_;

	# -----------------------------

	if ( $file_in ) {

		print YELLOW, "    Results ($file_in) present, skipping.\n", RESET;

	}

	else { print YELLOW, "    Result files present, skipping.\n", RESET; }

	# -----------------------------

	# End subroutine:
	return;

}


# ---------------------------------------------------------------------------- #

# Subroutine:	Comms::help
# Description:	Prints help message

# --------------------------------------

sub help {

	my ( $hash_congif ) = @_;

	# ------------------------------

	# Define help message:
	my $var_string1	= "EPIC workflow\n"; 
	$var_string1 .= "  Run [OPTIONS]\n";
	$var_string1 .= "	-h, Print this help\n";
	$var_string1 .= "	-i, Print Analysis file structure\n";
	$var_string1 .= "	-c, Run complete check\n";
	$var_string1 .= "	-f, Forked mode (default)\n";
	$var_string1 .= "	-l, Linear mode\n";
	$var_string1 .= "	-s, Sample IDs to run (list or CSV)\n";
	$var_string1 .= "	-r, RGI CARD database (default)\n";
	$var_string1 .= "	-w, RGI WildCARD database\n";
	$var_string1 .= "	-u, Submit Slurm jobs (under development)\n";
	$var_string1 .= "	-a, Run end-stage analysis\n\n";

	$var_string1 .= "    Individual stages:\n";
	$var_string1 .= "	-b, Subsample original reads\n";
	$var_string1 .= "	-p, trim and filter reads\n";
	$var_string1 .= "	-q, perform quality control on reads\n";
	$var_string1 .= "	-y, perform taxonomic and metabolic profiling\n";
	$var_string1 .= "	-e, perform ARG prediction\n";
	$var_string1 .= "	-m, perform assembly and annotation\n";

	print WHITE, "$var_string1\n";


	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Comms::num_contigs
# Description:	Prints message to user stating number of contigs present.

# --------------------------------------

sub num_contigs {

	# Arguments:
	my ( $file_in, $var_min_contigs ) = @_;

	# Variables:
	my $var_num_contigs;

	# ------------------------------

	$var_num_contigs = `grep ">" $file_in | wc -l`;

	chomp ( $var_num_contigs );

	my $var_string = "$var_num_contigs contigs (> $var_min_contigs bp)";
	$var_string .= " present in $file_in.\n";

	print WHITE, "    $var_string", RESET;

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Comms::print_analysis
# Description:	Prints analysis file structure

# --------------------------------------

sub print_analysis {

	# Arguments:
	my ( $hash_config ) = @_;

	# Data structures:
	my $hash_filenames	= $hash_config->{filenames};
	my $hash_analysis	= $hash_config->{filenames}->{analysis};

	# ------------------------------

	# Iterate through %hash_analysis:
	foreach my $var_term ( sort keys %{$hash_analysis} ) {

		Comms::print_command ( "Predicted Humann metabolic gene family profiles:" ) if $var_term eq "humann_gene_families";
		Comms::print_command ( "Predicted Humann metabolic pathway abundance profiles:" ) if $var_term eq "humann_path_abundance";

		Comms::print_command ( "Predicted ARG resistance mechanisms:" ) if $var_term eq "metaphlan_tax";
		Comms::print_command ( "Predicted ARG AMR families:" ) if $var_term eq "rgi_AMR_fam";
		Comms::print_command ( "Predicted ARG ARO terms:" ) if $var_term eq "rgi_ARO_term";

		Comms::print_command ( "Predicted ARG drug classes:" ) if $var_term eq "rgi_drug_class";
		Comms::print_command ( "Predicted ARG reference sequences:" ) if $var_term eq "rgi_ref_seq";
		Comms::print_command ( "Predicted ARG resistance mechanisms:" ) if $var_term eq "rgi_res_mech";

		foreach my $var_tag ( sort keys %{$hash_analysis->{$var_term}} ) {

			my $var_file	= $hash_analysis->{$var_term}->{$var_tag};

			print "      - $var_file\n";

		}

		print "\n";

	}

	Comms::print_command ( "Other files:" );
	print "      - $hash_filenames->{file_reads}\n";
	print "      - $hash_filenames->{file_phylum_abundance}\n";
	print "      - $hash_filenames->{file_metadata}\n";

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Error::pair_missing
# Description:	Warns user that a file is missing sequencing mate.

# --------------------------------------

sub pair_missing {

	# Arguments:
	my ( $var_string ) = @_;

	# ------------------------------

	print RED, "    ERROR: File pair for $var_string not detected!\n", RESET;

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #


# Subroutine:	Comms::print_command
# Description:	Prints command line message being run.

# --------------------------------------

sub print_command {

	# Arguments:
	my ( $var_string, $file_log ) = @_;

	# ------------------------------

	print WHITE, "    $var_string\n", RESET;

	# Print to $file_log if defined:
	if ( $file_log ) {

		Edit::append_string ( $file_log, "$var_string\n" );

	}

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	start_finish
# Description:	Prints start message to user and to log (optional).

# --------------------------------------

sub start_finish {

	# Arguments:
	my ( $var_process, $var_sample, $var_type, $file_log ) = @_;

	# Variables:
	my $var_date;
	my $var_time;
	my $var_string;

	# ------------------------------

	# Define current date and time:
	( $var_date, $var_time ) = Comms::date_time;

	# Print date and time:
	print "  $var_date $var_time";

	# Print start message:
	if ( $var_type eq "start" ) {

		$var_string = "Starting process $var_process for $var_sample:", RESET;

		print BLUE, "  $var_string\n";

	}

	# Print stop message:
	if ( $var_type eq "finish" ) {

		$var_string = "Finished process $var_process for $var_sample:", RESET;

		print GREEN, "  $var_string\n\n";

	}

	# Print to $file_log if defined:
	if ( $file_log ) {

		Edit::append_string ( $file_log, "$var_date $var_time $var_string\n" );

	}

	# ------------------------------

	# Return $var_string and end subroutine:
	return $var_string;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Error::terminate_code
# Description:	Prints message to use that process did not exit cleanly

# --------------------------------------

sub terminate_code {

	# Argments:
	my ( $var_in ) = @_;

	# ------------------------------

	print RED, "    ERROR: $var_in did not exit cleanly, terminating!\n", RESET;

	# ------------------------------

	# End subroutine;
	return;

}

# ---------------------------------------------------------------------------- #

1;
