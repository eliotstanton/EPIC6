package Stats;

use warnings;
use strict;

# ---------------------------------------------------------------------------- #

# Module name:  Stats
# Created by:   Eliot Stanton (estanton@wisc.edu)
# Created on:   13 September, 2023
# Modified:     19 September, 2023
# Description:  Handles end-point statistics primarily by coordinating calling R
#		scrpts.

# ---------------------------------------------------------------------------- #

# Subroutine:	abundance
# Description:	Perform differential abundance analysis for taxonomic, ARG, and
#		metabolic data.

# --------------------------------------

sub abundance {

	# Arguments:
	my ( $hash_config ) = @_;

	# Data structures:
	my $hash_filenames		= $hash_config->{filenames};
	my $hash_analysis		= $hash_config->{filenames}->{analysis};

	# File paths:
	my $dir_analysis		= $hash_config->{directories}->{analysis};
	my $dir_out			= "$dir_analysis/Maaslin2";
	my $file_metadata		= $hash_filenames->{file_metadata};
	my $flag_wildcard		= $hash_config->{opts}->{flag_wildcard};

	# ------------------------------

	# Iterate through terms in $hash_analysis:
	foreach my $var_term ( sort keys %{$hash_analysis} ) {

		next if $var_term =~ /rgi_ref_seq/;
		next if $var_term =~ /rgi_drug_class/;
		next if $var_term =~ /rgi_AMR_fam/;
		next if $var_term =~ /rgi_ARO_term/;
		next if $var_term =~ /rgi_res_mech/;
		next if $var_term =~ /metaphlan_tax/;
		next if $var_term =~ /humann_path/;
		next if $var_term =~ /humann_gene/;

		# Define output directory for Maaslin2:
		my $dir_out	= "$dir_out/$var_term";
		$dir_out	= "$dir_out\_wildcard" if $var_term =~ /rgi/ && $flag_wildcard;

		# Create Maaslin2 output directory:
		system ( "mkdir -p $dir_out" );

		# Define string for calling abundance R script:
		my $var_string1 = "r-base \\\n";
		$var_string1 .= "	workflow/stats/abundance.R \\\n";
		$var_string1 .= "		$file_metadata \\\n";

#		$var_string1 .= "		$hash_analysis->{$var_term}->{file_merged_abs} \\\n" unless $var_term =~ /humann/;
#		$var_string1 .= "		$hash_analysis->{$var_term}->{file_merged_cpm} \\\n" if $var_term =~ /humann/;

		$var_string1 .= "		$hash_analysis->{$var_term}->{file_merged_norm} \\\n" unless $var_term eq "metaphlan_tax";
		$var_string1 .= "		$hash_analysis->{$var_term}->{file_merged_abs} \\\n" if $var_term eq "metaphlan_tax";

		$var_string1 .= "		$hash_analysis->{$var_term}->{file_DAA} \\\n";

		$var_string1 .= "		$hash_analysis->{$var_term}->{file_merged_norm} \\\n" unless $var_term =~ /metaphlan_tax/;
		$var_string1 .= "		$hash_analysis->{$var_term}->{file_merged_rel} \\\n" if $var_term =~ /metaphlan_tax/;
		$var_string1 .= "		$dir_out";

		# Print command for calling R script to user:
		Comms::print_command ( $var_string1 );

		# Run command:
		system ($var_string1);

	}

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	diversity
# Description:	Handles calculating alpha and beta diversity for taxonomic,
#		ARG, and metabolic data

# --------------------------------------

sub diversity {

	# Arguments:
	my ( $hash_config )	= @_;

	# Data structures:
	my $hash_filenames	= $hash_config->{filenames};
	my $hash_analysis	= $hash_config->{filenames}->{analysis};

	# Filepaths:
	my $file_metadata		= $hash_filenames->{file_metadata};
	my $file_phylum_abundance	= $hash_filenames->{file_phylum_abundance};
	my $file_phylo_pdf		= $hash_filenames->{file_phylo_pdf};

	# ------------------------------

	# Loop through calculating alpha and beta diversity:
	foreach my $var_term ( sort keys %{$hash_analysis} ) {

		print "$var_term:\n";

		# Define string to assess alpha diversity:
		my $var_string1 = "r-base \\\n";
		$var_string1 .= "	workflow/stats/alpha_div.R \\\n";
		$var_string1 .= "		$file_metadata \\\n";
		$var_string1 .= "		$hash_analysis->{$var_term}->{file_merged_rel} \\\n";
		$var_string1 .= "		$hash_analysis->{$var_term}->{file_alpha_div} \\\n";
		$var_string1 .= "		$hash_analysis->{$var_term}->{file_alpha_div_pdf}";

		# Define string to assess beta diversity:
		my $var_string2 = "r-base \\\n";
		$var_string2 .= "	workflow/stats/beta_div.R \\\n";
		$var_string2 .= "		$file_metadata \\\n";
		$var_string2 .= "		$hash_analysis->{$var_term}->{file_merged_rel} \\\n";
		$var_string2 .= "		$hash_analysis->{$var_term}->{file_beta_div} \\\n";
		$var_string2 .= "		$hash_analysis->{$var_term}->{file_beta_div_pdf}";

		# Define string for generating abundance plot:
		my $var_string3 = "r-base \\\n";
		$var_string3 .= "	workflow/stats/ARGs_abundance.R \\\n";
		$var_string3 .= "	$hash_analysis->{$var_term}->{file_merged_norm} \\\n";
		$var_string3 .= "	$file_metadata \\\n";
		$var_string3 .= "	$hash_analysis->{$var_term}->{file_abundance_pdf }";

		# Print command to assess alpha diverity:
		Comms::print_command ( $var_string1 );

		# Run command:
		system ($var_string1);

		# Print command to assess beta diversity:
		Comms::print_command ( $var_string2 );

		# Run command:
		system ($var_string2);

		# Move on unless working with RGI data:
		next unless $var_term =~ /rgi/;

		# Print command to generate abundance plot:
		Comms::print_command ( $var_string3 );

		# Run command:
		system ( $var_string3 );

	}

	# ------------------------------

	# Define string for phylum abundance figure:
	my $var_string4 = "r-base \\\n";
	$var_string4 .= "	workflow/stats/phylum_abundance.R \\\n";
	$var_string4 .= "		$hash_analysis->{metaphlan_tax}->{file_merged_rel} \\\n";
	$var_string4 .= "		$file_metadata \\\n";
	$var_string4 .= "		$file_phylum_abundance";

	print "phylum_abundance:\n";

	# Print command for generating phylum abundance figure:
	Comms::print_command ( $var_string4 );

	# Run command for generating phylum abundance figure::
	system ($var_string4);

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	process
# Description:	Processes DAA output and generate LFC figures.

# --------------------------------------

sub process {

	# Arguments:
	my ( $hash_config ) = @_;

	# Data structures:
	my $hash_analysis	= $hash_config->{filenames}->{analysis};

	# File paths:

	# Variables:

	# ------------------------------

	# Iterate through DAA output files:
	foreach my $var_term ( sort keys %{$hash_analysis} ) {

#		next if $var_term =~ /rgi_ref_seq/;
#		next if $var_term =~ /rgi_drug_class/;
#		next if $var_term =~ /rgi_AMR_fam/;
#		next if $var_term =~ /rgi_ARO_term/;
#		next if $var_term =~ /rgi_res_mech/;
		next if $var_term =~ /metaphlan_tax/;
		next if $var_term =~ /humann_path/;
#		next if $var_term =~ /humann_gene/;

		# Define path for file to be processed:
		my $file_in = $hash_analysis->{$var_term}->{file_DAA};

		# Define path for output file:

		print "$var_term - $file_in\n";

		# Define @array_out to hold output:
		my @array_out;

		# Import contents of $file_in to @array_in:
		my @array_in = @{Edit::file_to_array ( $file_in )};

		# Store header as a string:
		my $var_header = $array_in[1];

		# Remove the first two lines from @array_in:
		splice @array_in, 0, 2;

		# Iterate through @array_in and convert each line into an
		# individual array:
		for ( my $i = 0; $i < scalar @array_in; $i++ ) {

			my @array_temp	= ( split " ", $array_in[$i] );

			unless ( $array_temp[1] ) {

				splice @array_in, $i, 1;

				$i--;

				next;

			}

			else {

				$array_in[$i]	= \@array_temp;

			}

		}

#		print "$var_header\n";

		# Sort by the 15th and 16th elements in each line:
		@array_in = sort { $a->[17] cmp $b->[17] || $b->[18] cmp $a->[18] || $b->[19] cmp $a->[19] } @array_in;

		for ( my $i = 0; $i < scalar @array_in; $i++ ) {

			next if $array_in[$i][17] eq "TRUE";

			if ( $array_in[$i][18] eq "TRUE" || $array_in[$i][19] eq "TRUE" ) {

				print "$i: @{$array_in[$i]}\n";

			}

		}		

#		print "@{$array_in[0]}\n";

		# Iterate through @array_in and export appropriate elements to
		# @array_out

		# Print @array_out to $file_out:

		# Call R script to generate figure from $file_out:

	}

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

1;
