package Stats;

use warnings;
use strict;

# ---------------------------------------------------------------------------- #

# Module name:  Stats
# Created by:   Eliot Stanton (estanton@wisc.edu)
# Created on:   13 September, 2023
# Modified:     29 October, 2023
# Description:  Handles end-point statistics primarily by coordinating calling R
#		scrpts.

# ---------------------------------------------------------------------------- #

# Subroutines:
# - abundance:	Calls R script for performing differential abundance analysis
# - divesity:	Calls R scripts for evaluating alpha and beta diversity
# - process:	Processes DAA output to identify significant results for user

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

	# Variables:
	my $flag_slurm			= $hash_config->{opts}->{flag_slurm};
	my $var_process_abundance	= "8";
	my $obj_process_abundance	= Parallel::ForkManager->new($var_process_abundance);

	# ------------------------------

	# Iterate through terms in $hash_analysis:
	foreach my $var_term ( sort keys %{$hash_analysis} ) {

                # Start forking:
                $obj_process_abundance->start and next;

#		next if $var_term =~ /rgi_ref_seq/;
#		next if $var_term =~ /rgi_drug_class/;
#		next if $var_term =~ /rgi_AMR_fam/;
#		next if $var_term =~ /rgi_ARO_term/;
#		next if $var_term =~ /rgi_res_mech/;
#		next if $var_term =~ /metaphlan_tax/;
#		next if $var_term =~ /humann_path/;
#		next if $var_term =~ /humann_gene/;

		# Define output directory for Maaslin2:
		my $dir_out	= "$dir_out/$var_term";
		$dir_out	= "$dir_out\_wildcard" if $var_term =~ /rgi/ && $flag_wildcard;

		# Create Maaslin2 output directory:
		system ( "mkdir -p $dir_out" );

		# Define string for calling abundance R script:
		my $var_string1 = "r-base \\\n";
		$var_string1 .= "	workflow/stats/abundance.R \\\n";
		$var_string1 .= "		$file_metadata \\\n";

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

		# End forking:
		$obj_process_abundance->finish;

	}

	# Wait for all children from $var_process_humann:
	$obj_process_abundance->wait_all_children;

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

	# Variables:
	my $flag_slurm			= $hash_config->{opts}->{flag_slurm};
	my $var_process_diversity	= "8";
	my $obj_process_diversity	= Parallel::ForkManager->new($var_process_diversity);

	# ------------------------------

	# Loop through calculating alpha and beta diversity:
	foreach my $var_term ( sort keys %{$hash_analysis} ) {

#		print "$var_term:\n";
                # Start forking:
		$obj_process_diversity->start and next;

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

		# Print command to assess alpha diverity:
		Comms::print_command ( $var_string1 );

		# Run command:
		system ($var_string1);

		# Print command to assess beta diversity:
		Comms::print_command ( $var_string2 );

		# Run command:
		system ($var_string2);

		# Generate abundance plots for ARGs:
		if ( $var_term =~ /rgi/ ) {

			# Define string for generating abundance plot:
			my $var_string3 = "r-base \\\n";
			$var_string3 .= "	workflow/stats/ARGs_abundance.R \\\n";
			$var_string3 .= "	$hash_analysis->{$var_term}->{file_merged_norm} \\\n";
			$var_string3 .= "	$file_metadata \\\n";
			$var_string3 .= "	$hash_analysis->{$var_term}->{file_abundance_pdf }";

			# Print command to generate abundance plot if working wit rgi:
			Comms::print_command ( $var_string3 ) if $var_term =~ /rgi/;

			# Run command for generating rgi plots:
			system ( $var_string3 ) if $var_term =~ /rgi/;

		}

		# End forking:
		$obj_process_diversity->finish;

	}

	# Wait for all children from $var_process_humann:
	$obj_process_diversity->wait_all_children;

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

# Subroutine:	process_DAA
# Description:	Processes DAA output and generate LFC figures.

# --------------------------------------

sub process_DAA {

	# Arguments:
	my ( $hash_config ) = @_;

	# Data structures:
	my $hash_analysis	= $hash_config->{filenames}->{analysis};

	# Variables:
	my $flag_slurm		= $hash_config->{opts}->{flag_slurm};

	# ------------------------------

	# Iterate through DAA output files:
	foreach my $var_term ( sort keys %{$hash_analysis} ) {

#		next if $var_term =~ /rgi_ref_seq/;
#		next if $var_term =~ /rgi_drug_class/;
#		next if $var_term =~ /rgi_AMR_fam/;
#		next if $var_term =~ /rgi_ARO_term/;
#		next if $var_term =~ /rgi_res_mech/;
#		next if $var_term =~ /metaphlan_tax/;
#		next if $var_term =~ /humann_path/;
#		next if $var_term =~ /humann_gene/;

		# Define path for output files:
		my $file_out	= $hash_analysis->{$var_term}->{"file_DAA_processed"};
		my $file_plot	= $hash_analysis->{$var_term}->{"file_DAA_pdf"};

		# Define path for file to be processed:
		my $file_in = $hash_analysis->{$var_term}->{file_DAA};

		# Define @array_out to hold output:
		my @array_out;

		# Import contents of $file_in to @array_in:
		my @array_in = @{Edit::file_to_array ( $file_in )};

		# Remove everything until [[1]]$res is detected:
		for ( my $i = 0; $i < scalar @array_in; $i++ ) {

			#print "$array_in[$i]\n";

			last if $array_in[$i] =~ /\$res/;

			splice @array_in, $i, 1;

			$i--;

		}

		# Store header as a string:
		my $var_header = $array_in[1];

		# Remove the first two lines from @array_in:
		splice @array_in, 0, 2;

		# Iterate through @array_in and convert each line into an
		# individual array:
		for ( my $i = 0; $i < scalar @array_in; $i++ ) {

			# Move on if TRUE isn't detected:
			next unless $array_in[$i] =~ /TRUE/ || $array_in[$i] =~ /FALSE/;

			$array_in[$i] = substring $array_in[$i], -3 if $array_in[$i] =~ /^s__/;

			# Split line into temporary array:
			my @array_temp	= ( split " ", $array_in[$i] );

			# Remove metaphlan formatting:
			if ( $array_temp[1] =~ /^s__/ ) {

				$array_temp[1] = substr $array_temp[1], 3;

			}

			# Move on unless second element in temporary array is defined:
			unless ( $array_temp[1] ) {

				next;

			}

			# Store temporary array in @array_out:
			else {

				push @array_out, \@array_temp;

			}

		}

		# Sort by the 15th and 16th elements in each line:
		@array_out = sort { $a->[17] cmp $b->[17] || $b->[18] cmp $a->[18] || $b->[19] cmp $a->[19] } @array_out;

		# Define string to hold data to be printed to file:
		my $var_string;

		# Iterate through @array_out adding each line to $var_string:
		for ( my $i = 0; $i < scalar @array_out; $i++ ) {

			# Add lines with significant results:
			if ( $array_out[$i][18] eq "TRUE" || $array_out[$i][19] eq "TRUE" ) {

				$var_string .= "@{$array_out[$i]}\n";

			}

		}

		# Add header if $var_string is defined:
		$var_string = "$var_header\n$var_string" if $var_string;

		# Print $var_string in $file_out:
		Edit::write_string ( $file_out, $var_string ) if $var_string;

		# Define string for calling R script that generates waterfall
		# figure:
		my $var_string2 = "r-base \\\n";
		$var_string2 .= "	workflow/stats/lfc.R \\\n";
		$var_string2 .= "	$file_out \\\n";
		$var_string2 .= "	$file_plot \\\n";

		# Call R script to generate figure from $file_out:
		system ( $var_string2 ) if $var_string;

	}

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

1;
