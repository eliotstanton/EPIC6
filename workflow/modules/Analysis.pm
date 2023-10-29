package Analysis;

use strict;
use warnings;
use Cwd;

# ---------------------------------------------------------------------------- #

# Module name:	Analysis
# Created by:	Eliot Stanton (estanton@wisc.edu)
# Created on:	13 September, 2023
# Modified:	29 October, 2023
# Description:	Handles end-point analysis

# ---------------------------------------------------------------------------- #

# Subroutines:
# - condense_humann	Condenses humann output files by removing individual
#			inputs.
# - count_reads		Counds reads present in original, trimmed, filtered
#			and deduplicated FASTQ files.
# - merge_metaphlan:	Merges individual metaphlan output for relative and
#			estimated read coverage.
# - merge_rgi		Merged individual rgi output for relative and absolute
#			read coverage.
# - print_reads		Prints hash containing number of reads in FASTQ files.
# - process_humann	Regroups, renames, and merges individual human output.

# ---------------------------------------------------------------------------- #

# Subroutine name:	Analysis::condense_humann
# Description:		Condenses humann output files by removing individual
#			inputs.

# --------------------------------------

sub condense_humann {

	# Arguments:
	my ( $file_in, $file_out ) = @_;

	# Data structures:
	my @array_in;
	my @array_out;

	# Variables:
	my $dir_in;
	my $var_string;

	# ------------------------------

	# Import $file_in to array:
	@array_in = @{Edit::file_to_array ( $file_in )};

	# Iterate through @array_in storing each gene family to @array_out:
	for ( my $i = 1; $i < scalar @array_in; $i++ ) {

		# Create temporary array:
		my @array_temp = ( split /[\t]/, $array_in[$i] );

		# Add line to @array_out if it doesn't contain bar characters:
		unless ( $array_temp[0] =~ /\|/ ) {

			#print "$i: @array_temp\n"

			push @array_out, $array_in[$i];

		}

	}

	# Write @array_out to $file_out:
	$var_string	= "$array_in[0]\n";
	$var_string	.= join "\n", @array_out;
	Edit::write_string ( $file_out, $var_string );

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Analysis::count_reads
# Description:	Counts the reads present in original, trimmed, filtered, and
#		deduplicated FASTQ files.

# --------------------------------------

sub count_reads {

	# Arguments:
	my ( $hash_config, $var_name_ID, $var_threads ) = @_;

	# Data structure:
	my $hash_files		= $hash_config->{filenames}->{$var_name_ID};
	my $hash_num;
	my $hash_out;

	# Filenames:
	my $file_fastq_R1	= $hash_files->{file_fastq_R1};
	my $file_trimmed_R1	= $hash_files->{file_trimmed_R1};
	my $file_filtered_R1	= $hash_files->{file_filtered_R1};
#	my $file_dedupe_R1	= $hash_files->{file_dedupe_R1};

	# Variables:
	$hash_num->{$file_fastq_R1}	= "1";
	$hash_num->{$file_trimmed_R1}	= "1";
	$hash_num->{$file_filtered_R1}	= "1";
#	$hash_num->{$file_dedupe_R1}	= "1";

	# ------------------------------

	# Iterate through keys in %hash_out:
	foreach my $file_in ( keys %{$hash_num} ) {

		next unless -e $file_in;

		# Define string for decompressing file and reading the
		# number of lines:
		my $var_string = "pigz \\\n";
		$var_string .= "	-dc -p $var_threads \\\n";
		$var_string .= "	$file_in | wc -l ";

		# Print command to user:
		Comms::print_command ( "$var_string" );

		# Run commands:
		$hash_num->{$file_in}	= `$var_string`/4;

#		print "$file_in: $hash_num->{$file_in}\n";

	}

	# Store $hash_num as a reference using $var_name_ID as the value:
	$hash_out->{$var_name_ID} = $hash_num;

	# ------------------------------

	# Return $hash_out and end subroutine:
	return $hash_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Analysis::merge_metaphlan
# Description:	Gathers and mergers metaphlan output data into files for
#		relative baundnace and estimated read coverage.

# --------------------------------------

sub merge_metaphlan {

	# Arguments:
	my ( $hash_config ) = @_;

	# Data-structures:
	my @array_in;
	my $hash_clade;
	my $hash_analysis		= $hash_config->{filenames}->{analysis};

	# File paths:
	my $dir_metaphlan		= $hash_config->{subdirectories}->{metaphlan};
	my $dir_analysis		= $hash_config->{directories}->{analysis};

	my $file_metaphlan_abs		= $hash_analysis->{metaphlan_tax}->{file_merged_abs};
	my $file_metaphlan_rel		= $hash_analysis->{metaphlan_tax}->{file_merged_rel};

	# Variables:
	my $var_string	= ""; #"Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies";
	my $var_string_rel;
	my $var_string_abs;

	# ------------------------------

	# Import contents of $dir_metaphlan to @array_in:
	@array_in	= @{ Edit::dir_to_array($dir_metaphlan) };

	# Remove non-tsv files from @array_in:
	@array_in = grep {/\.tsv$/} @array_in;

	# Sort @array_in:
	@array_in = sort { $a cmp $b } @array_in;

	# Iterate through @array_in and store relative and absolute abundance
	# in hashes:
	foreach my $file_in ( @array_in ) {

		# Remove extra parts from $file_in:
		my $file_in_edit = ( split "_metaphlan", $file_in )[0];

		# Add first line to $var_string to act as header with file names:
		$var_string .= "\t$file_in_edit";

#		print "$var_string\n";

		# Import $file_in to @array_file:
		my @array_file = @{ Edit::file_to_array ( "$dir_metaphlan/$file_in" ) };

#		print "@array_file\n";

		# Iterate through @array_file to extract relative and absolute
		# abundance and store in $hash_clade:
		foreach my $var_line ( @array_file ) {

#			print "$var_line\n";

			# Skip comment lines:
			next if $var_line =~/^\#/;

			# Skip lines above species level:
			next unless $var_line =~/\|s__/;

			# Skip lines below species level:
			next if $var_line =~/\|t__/;

			#print "\t$var_line\n";

			# Split lines into temperary arrays:
			my @array_temp = split "\t", $var_line;

			# Define clade, relative, and absolute abundance:
			my $var_clade		= $array_temp[0];
			my $var_rel_abund	= $array_temp[2];
			my $var_abs_abund	= $array_temp[4];

			# Split $var_clade to define domain, phylum, class,
			# family, genus, species, and strain:
			@array_temp		= split "\\|", $var_clade;

			# Define the number of elements in @array_temp:
			my $var_rank		= scalar @array_temp;

			# Store relative and absolute abundance in $hash_clade:
			$hash_clade->{$var_rank}->{$var_clade}->{$file_in}->{var_rel} = $var_rel_abund;
			$hash_clade->{$var_rank}->{$var_clade}->{$file_in}->{var_abs} = $var_abs_abund;

		}

	}

	# ------------------------------

	# Add line break to $var_string and copy over to $var_string_rel and 
	# $var_string_abs:
	$var_string	.= "\n";
	$var_string_rel	= $var_string;
	$var_string_abs	= $var_string;
	$var_string	= "";

	# Iterate through sorted $hash_clade adding clades and abundance
	# data to $var_string_rel and $var_string_abs:
	foreach my $var_rank ( sort keys %{$hash_clade} ) {

		# Iterate through clades:
		foreach my $var_clade ( sort keys %{$hash_clade->{$var_rank}} ) {

			# Add $var_clade to both $var_string_rel and
			# $var_string_abs:
			$var_string_rel .= "$var_clade";
			$var_string_abs .= "$var_clade";

			# Iterate through the files in @array_in:
			foreach my $file_in ( @array_in ) {

				# Define $hash_file:
				my $hash_file = $hash_clade->{$var_rank}->{$var_clade}->{$file_in};

				# Define default abundance values:
				my $var_rel_abund	= 0;
				my $var_abs_abund	= 0;

				# Add relative and absolute abundance values if they
				# exist:
				$var_rel_abund = $hash_file->{var_rel} if $hash_file->{var_rel};
				$var_abs_abund = $hash_file->{var_abs} if $hash_file->{var_abs};

				# Add abundances to $var_string_rel and $var_string_abs:
				$var_string_rel .= "\t$var_rel_abund";
				$var_string_abs .= "\t$var_abs_abund";

			}

			# Add line breaks to strings:
			$var_string_rel .= "\n";
			$var_string_abs .= "\n";

		}

	}

	# Print $var_string to $file_out:
	Edit::write_string ( $file_metaphlan_rel, $var_string_rel );
	Edit::write_string ( $file_metaphlan_abs, $var_string_abs );

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	merge_rgi
# Description:	Handles analyzing rgi output.

# --------------------------------------

sub merge_rgi {

	# Arguments:
	my ( $hash_config, $array_fastq ) = @_;

	# Data structures:
	my @array_fastq = @{$array_fastq};
	my $hash_files	= $hash_config->{filenames};
	my $hash_info;
	my $hash_samples;
	my $hash_total;

	my $hash_analysis		= $hash_config->{filenames}->{analysis};

	# File paths:
	my $dir_rgi			= $hash_config->{subdirectories}->{rgi};

	my $file_ref_seq_rel		= $hash_analysis->{rgi_ref_seq}->{file_merged_rel};
	my $file_ref_seq_abs		= $hash_analysis->{rgi_ref_seq}->{file_merged_abs};
	my $file_ref_seq_norm		= $hash_analysis->{rgi_ref_seq}->{file_merged_norm};

	my $file_ARO_term_rel		= $hash_analysis->{rgi_ARO_term}->{file_merged_rel};
	my $file_ARO_term_abs		= $hash_analysis->{rgi_ARO_term}->{file_merged_abs};
	my $file_ARO_term_norm		= $hash_analysis->{rgi_ARO_term}->{file_merged_norm};

	my $file_AMR_fam_rel		= $hash_analysis->{rgi_AMR_fam}->{file_merged_rel};
	my $file_AMR_fam_abs		= $hash_analysis->{rgi_AMR_fam}->{file_merged_abs};
	my $file_AMR_fam_norm		= $hash_analysis->{rgi_AMR_fam}->{file_merged_norm};

	my $file_drug_class_rel		= $hash_analysis->{rgi_drug_class}->{file_merged_rel};
	my $file_drug_class_abs		= $hash_analysis->{rgi_drug_class}->{file_merged_abs};
	my $file_drug_class_norm	= $hash_analysis->{rgi_drug_class}->{file_merged_norm};

	my $file_res_mech_rel		= $hash_analysis->{rgi_res_mech}->{file_merged_rel};
	my $file_res_mech_abs		= $hash_analysis->{rgi_res_mech}->{file_merged_abs};
	my $file_res_mech_norm		= $hash_analysis->{rgi_res_mech}->{file_merged_norm};

	# Variables:
	my $flag_wildcard		= $hash_config->{opts}->{flag_wildcard};
	my $var_min_reads		= $hash_config->{variables}->{min_reads};
	my $var_sample_ID		= join "\t", @array_fastq;

	my $var_ref_seq_fpkm		= "\t$var_sample_ID";
	my $var_ref_seq_reads		= "\t$var_sample_ID";
	my $var_ref_seq_relative	= "\t$var_sample_ID";

	my $var_ARO_term_fpkm		= "\t$var_sample_ID";
	my $var_ARO_term_reads		= "\t$var_sample_ID";
	my $var_ARO_term_relative	= "\t$var_sample_ID";

	my $var_ARM_fam_fpkm		= "\t$var_sample_ID";
	my $var_ARM_fam_reads		= "\t$var_sample_ID";
	my $var_ARM_fam_relative	= "\t$var_sample_ID";

	my $var_drug_class_fpkm		= "\t$var_sample_ID";
	my $var_drug_class_reads	= "\t$var_sample_ID";
	my $var_drug_class_relative	= "\t$var_sample_ID";

	my $var_res_mech_fpkm		= "\t$var_sample_ID";
	my $var_res_mech_reads		= "\t$var_sample_ID";
	my $var_res_mech_relative	= "\t$var_sample_ID";

	# ------------------------------

	# Import data from rgi output files:
	foreach my $var_name_ID ( @array_fastq ) {

		# Define RGI output files:
		my $file_allele = $hash_files->{$var_name_ID}->{file_allele};
		my $file_stats	= $hash_files->{$var_name_ID}->{file_stats};

		# Check that files hold multiple lines of data:
		unless ( -e $file_allele && -e $file_stats ) {

			Comms::file_missing ( $file_allele ) unless -e $file_allele;
			Comms::file_missing ( $file_stats ) unless -e $file_stats;

			exit;

		}

		# Skip if files aren't present:
		next unless -e "$file_allele" && -e "$file_stats";

		# Import $file_allele to @array_file and @array_stats:
		my @array_file		= @{ Edit::file_to_array ( "$file_allele" ) };
		my @array_stats		= @{ Edit::file_to_array ( "$file_stats" ) };

		# Determine number of reads analyzed from @array_stats:
		my $var_total_reads	= ( split " ", $array_stats[5])[2];
		my $var_mapped_reads	= ( split " ", $array_stats[6])[2];

		# Iterate through data in @array_file:
		for ( my $i = 1; $i < scalar @array_file; $i ++ ) {

			my @array_data	= split /[\t]/, $array_file[$i];

			# Define each tab as a variable:
			my $var_ref_seq		 	= $array_data[0];
			my $var_ARO_term		= $array_data[1];
			my $var_ARO_accession		= $array_data[2];
			my $var_reference_model_type	= $array_data[3];
			my $var_reference_db		= $array_data[4];
			my $var_reference_allele_source = $array_data[5];
			my $var_observed_in_genomes	= $array_data[6];
			my $var_observed_in_plasmids	= $array_data[7];
			my $var_observed_in_pathogens	= $array_data[8];
			my $var_completely_mapped_reads = $array_data[9];
			my $var_mapped_reads_flank_seq	= $array_data[10];
			my $var_all_mapped_reads	= $array_data[11];
			my $var_percent_coverage	= $array_data[12];
			my $var_length_coverage	 	= $array_data[13];
			my $var_average_MAPQ		= $array_data[14];
			my $var_mate_pair_linkage	= $array_data[17];
			my $var_reference_length	= $array_data[16];
			my $var_AMR_fam		 	= $array_data[17];
			my $var_drug_class		= $array_data[18];
			my $var_res_mech		= $array_data[19];
			my $var_depth			= $array_data[20];
			my $var_SNPs			= $array_data[21];
			my $var_consensus_sequence_DNA	= $array_data[22];
			my $var_consensus_sequence_prot = $array_data[23];

			# Skip if number of reads is below threshold set by $var_min_reads:
			next if $var_all_mapped_reads < $var_min_reads;

			# Calculate normalized data from read counts:
			my $var_fpkm = ($var_all_mapped_reads*10**9) / ($var_reference_length*$var_total_reads);

			# Round $var_FPKM to six decimal places:
			$var_fpkm = sprintf("%.6f", $var_fpkm);

			# Store needed tags in $hash_info under ref_seq term:
			$hash_info->{ref_seq}->{$var_ref_seq}->{var_ARO_term}	= $var_ARO_term;
			$hash_info->{ref_seq}->{$var_ref_seq}->{var_AMR_fam}	= $var_AMR_fam;
			$hash_info->{ref_seq}->{$var_ref_seq}->{var_drug_class} = $var_drug_class;
			$hash_info->{ref_seq}->{$var_ref_seq}->{var_res_mech}	= $var_res_mech;

			# Store other tags in $hash_info:
			$hash_info->{ARO_term}->{$var_ARO_term}		 = $var_ARO_term;

			my @array_AMR_fam = split /; /, $var_AMR_fam;

			foreach my $var_AMR_fam ( @array_AMR_fam ) {

				$hash_info->{AMR_fam}->{$var_AMR_fam}			= $var_AMR_fam;

			}

			my @array_drug_class = split /; /, $var_drug_class;

			foreach my $var_drug_class ( @array_drug_class ) {

				$hash_info->{drug_class}->{$var_drug_class}	= $var_drug_class;

			}

			my @array_res_mech = split /; /, $var_res_mech;

			foreach my $var_res_mech ( @array_res_mech ) {

				$hash_info->{res_mech}->{$var_res_mech}	 = $var_res_mech;

			}

			# Define temporary array:
			my $hash_temp;

			# Store data for sample in $hash_temp:
			$hash_temp->{ref_seq}->{$var_ref_seq}->{var_fpkm}	= $var_fpkm;
			$hash_temp->{ref_seq}->{$var_ref_seq}->{var_reads}	= $var_all_mapped_reads;

			$hash_temp->{ARO_term}->{$var_ARO_term}->{var_fpkm}	= $var_fpkm;
			$hash_temp->{ARO_term}->{$var_ARO_term}->{var_reads}	= $var_all_mapped_reads;

			foreach my $var_AMR_fam ( @array_AMR_fam ) {

				$hash_temp->{AMR_fam}->{$var_AMR_fam}->{var_fpkm}	= $var_fpkm;
				$hash_temp->{AMR_fam}->{$var_AMR_fam}->{var_reads}	= $var_all_mapped_reads;

			}

			foreach my $var_drug_class ( @array_drug_class ) {

				$hash_temp->{drug_class}->{$var_drug_class}->{var_fpkm}	 = $var_fpkm;
				$hash_temp->{drug_class}->{$var_drug_class}->{var_reads}	= $var_all_mapped_reads;

			}

			foreach my $var_res_mech ( @array_res_mech ) {

				$hash_temp->{res_mech}->{$var_res_mech}->{var_fpkm}	= $var_fpkm;
				$hash_temp->{res_mech}->{$var_res_mech}->{var_reads}	= $var_all_mapped_reads;

			}

			# Iteratively store $hash_temp in $hash_sample:
			foreach my $var_term ( sort keys %{$hash_temp} ) {

				foreach my $var_tag ( sort keys %{$hash_temp->{$var_term}} ) {

					my $var_fpkm	= $hash_temp->{$var_term}->{$var_tag}->{var_fpkm};
					my $var_reads	= $hash_temp->{$var_term}->{$var_tag}->{var_reads};

					$hash_samples->{$var_name_ID}->{$var_term}->{$var_tag}->{var_fpkm}	+= $var_fpkm;
					$hash_samples->{$var_name_ID}->{$var_term}->{$var_tag}->{var_reads}	+= $var_reads;

				}

			}

		}

	}

	# ------------------------------

	# Iterate through each sample and record the total number of reads to calculate relative abundance:
	foreach my $var_term ( sort keys %{$hash_info} ) {

		foreach my $var_name_ID ( @array_fastq ) {

			# Define variable to hold total number of reads:
			my $var_total	= 0;


			# Iterate through list of tags:
			foreach my $var_tag ( sort keys %{$hash_info->{$var_term}} ) {

				my $var_reads = 0;

				if ( $hash_samples->{$var_name_ID}->{$var_term}->{$var_tag}->{var_reads} ) {

					$var_reads = $hash_samples->{$var_name_ID}->{$var_term}->{$var_tag}->{var_reads};

				}

				$var_total += $var_reads;

			}

			# Store total in $hash_total:
			$hash_total->{$var_name_ID}->{$var_term} = $var_total;

		}

	}

	# ------------------------------

	# Iterate through $hash_info creating tables:
	foreach my $var_term ( sort keys %{$hash_info} ) {

		# Deifne strings to hold data:
		my $var_string_fpkm;
		my $var_string_reads;
		my $var_string_relative;

		# Iterate through list of tags:
		foreach my $var_tag ( sort keys %{$hash_info->{$var_term}} ) {

			# Replace spaces with underscores:
			my @array_temp	= ( split /[\s,\|,\:]/, $var_tag );
			my $var_tag2	= join "_", @array_temp;

			# Add initial tags to strings:
			$var_string_fpkm .= $var_tag2;
			$var_string_reads .= $var_tag2;
			$var_string_relative .= $var_tag2;

			# Handle additional information under ref_seq term:
			if ( $var_term eq "ref_seq" ) {

				my $var_AMR_fam	 	= $hash_info->{$var_term}->{$var_tag}->{var_AMR_fam};
				my $var_ARO_term	= $hash_info->{$var_term}->{$var_tag}->{var_ARO_term};
				my $var_drug_class	= $hash_info->{$var_term}->{$var_tag}->{var_drug_class};
				my $var_res_mech	= $hash_info->{$var_term}->{$var_tag}->{var_res_mech};

			}

			# Iterate through each sample:
			foreach my $var_name_ID ( @array_fastq ) {

				# Define variables to hold data:
				my $var_fpkm		= "0";
				my $var_reads		= "0";
				my $var_relative	= "0";

				# Define total number of reads:
				my $var_total	= $hash_total->{$var_name_ID}->{$var_term};

				# If data is available redefine $var_fpkm and $var_reads:
				if ( $hash_samples->{$var_name_ID}->{$var_term}->{$var_tag}->{var_fpkm} ) {

					$var_fpkm = $hash_samples->{$var_name_ID}->{$var_term}->{$var_tag}->{var_fpkm};

				}

				if ( $hash_samples->{$var_name_ID}->{$var_term}->{$var_tag}->{var_reads} ) {

					$var_reads = $hash_samples->{$var_name_ID}->{$var_term}->{$var_tag}->{var_reads};

					$var_relative = sprintf("%.6f", $var_reads/$var_total );

				}

				# Add data to strings:
				$var_string_fpkm .= "\t$var_fpkm";
				$var_string_reads .= "\t$var_reads";
				$var_string_relative .= "\t$var_relative";

			}

			# Add line breaks to strings:
			$var_string_fpkm .= "\n";
			$var_string_reads .= "\n";
			$var_string_relative .= "\n";

		}

		# Store res_seq strings:
		if ( $var_term eq "ref_seq" ) {

			$var_ref_seq_fpkm = "# Normalized FPKM values\n$var_ref_seq_fpkm\n$var_string_fpkm";
			$var_ref_seq_reads = "# Read counts\n$var_ref_seq_reads\n$var_string_reads";

			$var_ref_seq_relative = "# Relative abundance\n$var_ref_seq_relative\n$var_string_relative";			

		}

		# Store ARO_term strings:
		elsif ( $var_term eq "ARO_term" ) {

			$var_ARO_term_fpkm = "# Normalized FPKM values\n$var_ARO_term_fpkm\n$var_string_fpkm";
			$var_ARO_term_reads = "# Read counts\n$var_ARO_term_reads\n$var_string_reads";
			$var_ARO_term_relative = "# Relative abundance\n$var_ARO_term_relative\n$var_string_relative";

		}

		# Store ARM_fam strings:
		elsif ( $var_term eq "AMR_fam" ) {

			$var_ARM_fam_fpkm = "# Normalized FPKM values\n$var_ARM_fam_fpkm\n$var_string_fpkm";
			$var_ARM_fam_reads = "# Read counts\n$var_ARM_fam_reads\n$var_string_reads";
			$var_ARM_fam_relative = "# Relative abundance\n$var_ARM_fam_relative\n$var_string_relative";

		}

		# Store drug_class strings:
		elsif ( $var_term eq "drug_class" ) {

			$var_drug_class_fpkm = "# Normalized FPKM values\n$var_drug_class_fpkm\n$var_string_fpkm";
			$var_drug_class_reads = "# Read counts\n$var_drug_class_reads\n$var_string_reads";
			$var_drug_class_relative = "# Relative abundance\n$var_drug_class_relative\n$var_string_relative";

		}

		# Store res_mech strings:
		elsif ( $var_term eq "res_mech" ) {

			$var_res_mech_fpkm = "# Normalized FPKM values\n$var_res_mech_fpkm\n$var_string_fpkm";
			$var_res_mech_reads = "# Read counts\n$var_res_mech_reads\n$var_string_reads";
			$var_res_mech_relative = "# Read counts\n$var_res_mech_relative\n$var_string_relative";

		}

	}

	# ------------------------------

	# Print strings to files:
	Edit::write_string ( $file_ref_seq_norm, $var_ref_seq_fpkm );
	Edit::write_string ( $file_ref_seq_abs, $var_ref_seq_reads );
	Edit::write_string ( $file_ref_seq_rel, $var_ref_seq_relative );

	Edit::write_string ( $file_ARO_term_norm, $var_ARO_term_fpkm );
	Edit::write_string ( $file_ARO_term_abs, $var_ARO_term_reads );
	Edit::write_string ( $file_ARO_term_rel, $var_ARO_term_relative );

	Edit::write_string ( $file_AMR_fam_norm, $var_ARM_fam_fpkm );
	Edit::write_string ( $file_AMR_fam_abs, $var_ARM_fam_reads );
	Edit::write_string ( $file_AMR_fam_rel, $var_ARM_fam_relative );

	Edit::write_string ( $file_drug_class_norm, $var_drug_class_fpkm );
	Edit::write_string ( $file_drug_class_abs, $var_drug_class_reads );
	Edit::write_string ( $file_drug_class_rel, $var_drug_class_relative );

	Edit::write_string ( $file_res_mech_norm, $var_res_mech_fpkm );
	Edit::write_string ( $file_res_mech_abs, $var_res_mech_reads );
	Edit::write_string ( $file_res_mech_rel, $var_res_mech_relative );

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Analysis::print_reads
# Description:	Handles printing hash containing number of reads for FASTQ files
#		to file:

# --------------------------------------

sub print_reads {

	# Arguments:
	my ( $hash_config, $hash_in, $file_reads ) = @_;

	# Data structures:
	my $hash_fastq	= $hash_config->{fastq};
	my @array_in;
	my @array_out;

	# Filepaths:
	my $dir_reads	= $hash_config->{subdirectories}->{reads};

	# ------------------------------

	# Transfer sorted data from $hash_in to @array_in:
	foreach my $var_name_ID ( sort keys %{$hash_in} ) {

		print "$var_name_ID\n";

		my $hash_files		= $hash_config->{filenames}->{$var_name_ID};
		my $file_fastq_R1	= $hash_files->{file_fastq_R1};
		my $file_trimmed_R1	= $hash_files->{file_trimmed_R1};
		my $file_filtered_R1	= $hash_files->{file_filtered_R1};

		my $hash_temp		= $hash_in->{$var_name_ID};
		my $var_num_fastq	= $hash_temp->{$file_fastq_R1};
		my $var_num_trimmed	= $hash_temp->{$file_trimmed_R1};
		my $var_num_filtered	= $hash_temp->{$file_filtered_R1};

		my $var_per_trimmed	= sprintf( "%.2f", ( $var_num_trimmed / $var_num_fastq ) * 100 ) if $var_num_trimmed;
		my $var_per_filtered	= sprintf( "%.2f", ( $var_num_filtered / $var_num_fastq ) * 100 ) if $var_num_filtered;

		# Define temporary arrays:
		my @array_temp1		 = ($file_fastq_R1, $var_num_fastq, "");
		my @array_temp2		 = ($file_trimmed_R1, $var_num_trimmed, $var_per_trimmed);
		my @array_temp3		 = ($file_filtered_R1, $var_num_filtered, $var_per_filtered);
		my @array_temp5		 = ("","","");

		print "@array_temp1\n@array_temp2\n@array_temp3\n\n"; #@array_temp4\n\n";

		# Store references to temporary arrays in @array_out:
		push @array_out, \@array_temp1;
		push @array_out, \@array_temp2;
		push @array_out, \@array_temp3;
		push @array_out, \@array_temp5;

		# Create a string with output and store individually in $dir_analysis:
		my $var_string = "@array_temp1\t@array_temp2\t";
		$var_string .= "@array_temp3\n";

		Edit::write_string ( "$var_name_ID\_reads.txt", $var_string );

	}

	# ------------------------------

	my $file_in; my $var_num_reads, my $var_surviving;

# Define body of printout:
format DATA =
    @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   @||||||||||||    @|||||||
$file_in, $var_num_reads, $var_surviving
.

my $var_header  = "    Filename					    Number reads    Survival (%)
    --------------------------------------------------------------------------------";

	# Open $file_surviving_log:
	open ( my $file_write, '>', "$file_reads" ) or die $!;

	select($file_write);
	$~ = "DATA";

	print $file_write "$var_header\n";

	# Iterate through and print out to $file_write:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		$file_in	= $array_out[$i][0];
		$var_num_reads	= $array_out[$i][1];
		$var_surviving	= $array_out[$i][2];

		write;

	}

	# Close $file_write:
	close $file_write;

	# ------------------------------

	select(STDOUT);
	$~ = "DATA";

	print "$var_header\n";

	# Iterate through and print out to stdout:
	for ( my $i = 0; $i < scalar @array_out; $i++ ) {

		$file_in	= $array_out[$i][0];
		$var_num_reads	= $array_out[$i][1];
		$var_surviving	= $array_out[$i][2];
		
		write;

	}

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:   process_humann
# Description:  Handles regrouping, renaming, and merging humann output.

# --------------------------------------

sub process_humann {

	# Arguments:
	my ( $hash_config ) = @_;

	# Data structures:
	my $hash_files		= $hash_config->{filenames};
	my $hash_fastq		= $hash_config->{fastq};
	my $hash_analysis	= $hash_config->{filenames}->{analysis};
	my @array_fastq		= @{Edit::keys_to_array ( $hash_config->{fastq} )};

	# File paths:
	my $dir_container_home	= $hash_config->{directories}->{container_home};
	my $dir_humann		= $hash_config->{subdirectories}->{humann};
	my $file_gene_rpk	= $hash_analysis->{humann_gene_families}->{file_merged_norm};
	my $file_gene_relab	= $hash_analysis->{humann_gene_families}->{file_merged_rel};
	my $file_gene_cpm	= $hash_analysis->{humann_gene_families}->{file_merged_cpm};
	my $file_path_rpk	= $hash_analysis->{humann_path_abundance}->{file_merged_norm};
	my $file_path_relab	= $hash_analysis->{humann_path_abundance}->{file_merged_rel};
	my $file_path_cpm	= $hash_analysis->{humann_path_abundance}->{file_merged_cpm};
	my $dir_cwd		= getcwd

	# Variables:
	my $var_humann_version  = $hash_config->{software}->{humann}->{tag};
	my $var_singularity     = "singularity exec --bind $dir_cwd ";
	$var_singularity        .= "$dir_container_home/humann_$var_humann_version\.sif \\\n";

	my $var_process_humann	= "48";
	my $obj_process_humann	= Parallel::ForkManager->new($var_process_humann);
	my $var_header		= join ( "\t", @array_fastq );
	$var_header		= "\t$var_header\n";

	# ------------------------------

	# Loop through humann output files:
	foreach my $var_name_ID ( sort keys %{$hash_fastq} ) {

		# Start forking:
		$obj_process_humann->start and next;

		# Define individual filenames for input and output:
		my $file_gene		= $hash_files->{$var_name_ID}->{file_humann_gene};
		my $file_path		= $hash_files->{$var_name_ID}->{file_humann_path};
		my $file_gene_rename	= $hash_files->{$var_name_ID}->{file_humann_gene_rename};
		my $file_gene_cond	= $hash_files->{$var_name_ID}->{file_humann_gene_cond};
		my $file_path_cond	= $hash_files->{$var_name_ID}->{file_humann_path_cond};
		my $file_regroup	= $hash_files->{$var_name_ID}->{file_humann_regroup};
		my $file_gene_relab	= $hash_files->{$var_name_ID}->{file_humann_gene_relab};
		my $file_gene_cpm	= $hash_files->{$var_name_ID}->{file_humann_gene_cpm};
		my $file_path_relab	= $hash_files->{$var_name_ID}->{file_humann_path_relab};
		my $file_path_cpm	= $hash_files->{$var_name_ID}->{file_humann_path_cpm};

		# Check that all output files are present:
		unless ( -e $file_gene && -e $file_path ) {

			Comms::file_missing ( $file_gene ) unless -e $file_gene;
			Comms::file_missing ( $file_path ) unless -e $file_path;

			exit;

		}

		# Define string to regroup gene family file to functional
		# categories:
		my $var_string1 = $var_singularity;
		$var_string1 .= "	humann_regroup_table \\\n";
		$var_string1 .= "		--input $file_gene \\\n";
		$var_string1 .= "		--output $file_regroup \\\n";
		$var_string1 .= "		--groups uniref90_rxn \\\n";

                # Define string to rename features in gene family file:
		my $var_string2 = $var_singularity;
		$var_string2	.= "	humann_rename_table \\\n";
		$var_string2	.= "		--input $file_regroup \\\n";
		$var_string2	.= "		--output $file_gene_rename \\\n";
		$var_string2	.= "		--names metacyc-rxn \\\n";

		# Define primary string to renorm files:
		my $var_string3 = $var_singularity;
		$var_string3 .= "	humann_renorm_table \\\n";
		$var_string3 .= "		--special n \\\n";
		$var_string3 .= "		--update-snames \\\n";

		# Define secondary strings for renorming gene families and paths
		# abundance files:
		my $var_string4 = $var_string3;
		$var_string4 .= "		--input $file_gene_cond \\\n";

		my $var_string5 = $var_string3;
		$var_string5 .= "		--input $file_path_cond \\\n";

		# Define tertiary strings for renorming files to either relative
		# abundance or cpm:
		my $var_string6 = $var_string4;
		$var_string6 .= "		--output $file_gene_relab\\\n";
		$var_string6 .= "		--units relab";

		my $var_string7 = $var_string4;
		$var_string7 .= "		--output $file_gene_cpm \\\n";
		$var_string7 .= "		--units cpm";

		my $var_string8 = $var_string5;
		$var_string8 .= "		--output $file_path_relab \\\n";
		$var_string8 .= "		--units relab";

		my $var_string9 = $var_string5;
		$var_string9 .= "		--output $file_path_cpm \\\n";
		$var_string9 .= "		--units cpm";

		# ----------------------

		# Print command to regroup gene family file to user:
		Comms::print_command ( $var_string1 );

		# Run command to regroup gene family file:
		system ( $var_string1 );

		# Print command to rename regrouped gene family files to user:
		Comms::print_command ( $var_string2 );

		# Run command to rename regrouped gene family file:
		system ( $var_string2 );

		# Condense $file_rename and $file_path by removing the individual
		# microbial community inputs:
		condense_humann ( $file_gene_rename, $file_gene_cond );
		condense_humann ( $file_path, $file_path_cond );

		# Print commands for renorming to the user:
		Comms::print_command ( $var_string6 );
		Comms::print_command ( $var_string7 );
		Comms::print_command ( $var_string8 );
		Comms::print_command ( $var_string9 );

		# Run the commands for regrouping:
		system ( $var_string6 );
		system ( $var_string7 );
		system ( $var_string8 );
		system ( $var_string9 );

		# End forking:
		$obj_process_humann->finish;

	}

	# Wait for all children from $var_process_humann:
	$obj_process_humann->wait_all_children;

	# ------------------------------

	# Define primary string for calling join tool:
	my $var_string10 = $var_singularity;
	$var_string10	.= "	humann_join_tables \\\n";
	$var_string10	.= "		--input $dir_humann \\\n";

	# Define secondary strings for joining regrouped files:
	my $var_string11 = $var_string10;
	$var_string11	.= "		--file_name genefamilies_cond \\\n";
	$var_string11	.= "		--output $file_gene_rpk";

	my $var_string12 = $var_string10;
	$var_string12	.= "		--file_name pathabundance_cond \\\n";
	$var_string12	.= "		--output $file_path_rpk";

	# Define secondary strings for joining cpm files:
	my $var_string13 = $var_string10;
	$var_string13	.= "		--file_name genefamilies_cpm \\\n";
	$var_string13	.= "		--output $file_gene_cpm";

	my $var_string14 = $var_string10;
	$var_string14	.= "		--file_name pathabundance_cpm \\\n";
	$var_string14	.= "		--output $file_path_cpm";

	# Define secondary strings for joining abundance files:
	my $var_string15 = $var_string10;
	$var_string15	.= "		--file_name genefamilies_relab \\\n";
	$var_string15	.= "		--output $file_gene_relab";

	my $var_string16 = $var_string10;
	$var_string16	.= "		--file_name pathabundance_relab \\\n";
	$var_string16	.= "		--output $file_path_relab";

	# Print commands for joining tables to user:
	Comms::print_command ( $var_string11 );
	Comms::print_command ( $var_string12 );

	Comms::print_command ( $var_string13 );
	Comms::print_command ( $var_string14 );

	Comms::print_command ( $var_string15 );
	Comms::print_command ( $var_string16 );

	# Run commands for joining tables:
	system ( $var_string11 );
	system ( $var_string12 );

	system ( $var_string13 );
	system ( $var_string14 );

	system ( $var_string15 );
	system ( $var_string16 );

	# ------------------------------

	# Loop through $hash_analysis modifying joined tables to be compliant
	# with importing into R:
	foreach my $var_term ( keys %{$hash_analysis} ) {

		# Move on if the term doesn't relate to humann:
		next unless $var_term =~ /humann/;

		# Iterate through files:
		foreach my $var_tag ( keys %{$hash_analysis->{$var_term}} ) {

			# Move on if $var_tag doesn't relate to merged files:
			next unless $var_tag =~ /file_merged/;

			# Define string to hold outbound data:
			my $var_string;

			# Get file path:
			my $file_in	= $hash_analysis->{$var_term}->{$var_tag};

			# Convert file to array:
			my @array_file = @{Edit::file_to_array( $file_in )};

			# Replace first line of @array_file with $var_header:
			$array_file[0]	= $var_header;

			# Remove third line if beginning of line mateches UNGROUPED:
			splice @array_file, 2, 1 if $array_file[2] =~ /^UNGROUPED/;

			# Remove third line if beginning of line matches UNINTEGRATED:
			splice @array_file, 2, 1 if $array_file[2] =~ /^UNINTEGRATED/;

			# Remove third line if beginning of line matches UNMAPPED:
			splice @array_file, 1, 1 if $array_file[1] =~ /^UNMAPPED/;

			# Iterate through @array_file removing white spacess but
			# leaving tabs:
			foreach my $var_line ( @array_file ) {


				# Replace white spaces with underscores:
				$var_line =~ tr / /\_/;

				# Add current string to $var_string:
				$var_string .= "$var_line\n";

			}

			# Write string to file:
			Edit::write_string ( $file_in, $var_string );

		}

	}

}

# ---------------------------------------------------------------------------- #

1;
