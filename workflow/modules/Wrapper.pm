package Wrapper;

use strict;
use warnings;

# ---------------------------------------------------------------------------- #

# Module name:	Wrapper
# Created by:	Eliot Stanton (estanton@wisc.edu)
# Created on:	11 September, 2023
# Modified:	13 September, 2023
# Description:	Subroutines functioning as wrappers for calling external
#		programs.

# ---------------------------------------------------------------------------- #

# Subroutines:
# - run:	Handles checking that output dir is present, printing
#		information and running command with logging and stats output.
# - bwa		Read alignment
# - fastqc	Quality control
# - dedupe	Read deduplication
# - humann	Metabolic profiling
# - kraken2	Taxonomic classification
# - metaphlan	Taxonomic classification
# - quast	Assembly quality control
# - prokka	Assembly annotation
# - rgi_load	RGI database loading
# - rgi_bwt	ARG prediction
# - seqtk	Read subsampling
# - spades	Genome assembly
# - trimmomatic	Read trimming

# ---------------------------------------------------------------------------- #

# Subroutine: 	Wrapper::run
# Description:	Handles checking that output directory is present, printing
#		information to user and logs, and running the command.

# --------------------------------------

sub run {

	# Arguments:
	my ( $hash_config, $var_name_ID, $var_command, $var_string, $dir_out,
	$file_log ) = @_;

	# File paths:
	my $dir_logs	= $hash_config->{directories}->{logs};

	# Variables:
	my $var_stats	= $hash_config->{variables}->{stats};
	my $var_return;

	# ------------------------------

	# Check that output directory is present:
	system ( "mkdir -p $dir_out" ) unless -d "$dir_out";

	# Check that logs directory is present:
	system ( "mkdir -p $dir_logs/$var_name_ID" ) unless -d "$dir_logs/$var_name_ID";

	# Print container information:
	Comms::container_info ( $hash_config, $var_command, $file_log );

	# Print command:
	Comms::print_command ( $var_string, $file_log );

	# Run command:
	$var_return = system ( "$var_stats $var_string &>> $file_log" );

	# ------------------------------

	# TODO: Handle error messages

	# ------------------------------
	
	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	bwa
# Description:	Read filtering using aligner.

# --------------------------------------

sub bwa {

	# Arguments:
	my ( $hash_config, $var_name_ID, $file_trimmed_R1, $file_trimmed_R2,
	$file_concatenated, $file_bam, $file_filtered_R1, $file_filtered_R2,
	$file_log, $var_threads ) = @_;

	# Data structures:
	my $dir_logs	= $hash_config->{directories}->{logs};

	# File paths:
	my $dir_out	= $hash_config->{subdirectories}->{bwa};

	# Variables:
	my $var_stats	= $hash_config->{variables}->{stats};
	my $var_return;

	# ------------------------------

	# Define command for running bwa-mem:
	my $var_string1 = "bwa mem \\\n";
	$var_string1 .= "	-t $var_threads \\\n";
	$var_string1 .= "	$file_concatenated \\\n";
	$var_string1 .= "	$file_trimmed_R1 \\\n";
	$var_string1 .= "	$file_trimmed_R2 ";

	# Define command for running SAMtools sort:
	my $var_string2 = "samtools sort \\\n";
	$var_string2 .= "	-@ $var_threads \\\n";
	$var_string2 .= "	-o $file_bam";

	# Define command for running SAMtools fastq:
	my $var_string3 = "samtools fastq \\\n";
	$var_string3 .= "	-@ $var_threads \\\n";
	$var_string3 .= "	-f 12 \\\n";
	$var_string3 .= "	-1 $file_filtered_R1 \\\n";
	$var_string3 .= "	-2 $file_filtered_R2 \\\n";
	$var_string3 .= "	$file_bam";

	# Check if $file_filtered_R1 $file_filtered_R2 is present and has
	# non-zero size:
	if ( -s $file_filtered_R1 && -s $file_filtered_R2 ) {

		Comms::file_present ( $file_filtered_R1 );
		Comms::file_present ( $file_filtered_R2 );

		$var_return = 0;

	}

	else {

		# Check that $file_trimmed_R1 and $file_trimmed_R2 are present and
		# have non-zero sizes:
		unless ( -s $file_trimmed_R1 && -s $file_trimmed_R2 ) {

			Comms::file_missing ( $file_trimmed_R1 ) unless -e $file_trimmed_R1;
			Comms::file_missing ( $file_trimmed_R2 ) unless -e $file_trimmed_R2;

			$var_return = 1;

		}

		else {

			# Check that output directory is present:
			system ( "mkdir -p $dir_out" ) unless -d "$dir_out";

			# Check that logs directory is present:
			system ( "mkdir -p $dir_logs/$var_name_ID" ) unless -d "$dir_logs/$var_name_ID";

			# Print container information:
			Comms::container_info ( $hash_config, "bwa", $file_log );

			# Print command:
			Comms::print_command ( "$var_string1 2>> $file_log | $var_string2", $file_log );

			# Run command:
			$var_return = system ( "$var_stats $var_string1 2>> $file_log | $var_string2 &>> $file_log" );

			# Print container information:
			Comms::container_info ( $hash_config, "samtools", $file_log );

			# Print command:
			Comms::print_command ( $var_string3, $file_log );

			# Run command:
			$var_return = system ( "$var_stats $var_string3 &>> $file_log" )

			# Remove $file_bam:
#			system ( "rm $dir_out/$file_bam" ) if -e "$dir_out/$file_bam";			

		}

	}

	# ------------------------------

	# Handle error code from bwa:
	if ( $var_return != 0 ) {

		Comms::exit_code ( $file_log );

		system ( "rm $file_filtered_R1" ) if -e $file_filtered_R1;
		system ( "rm $file_filtered_R2" ) if -e $file_filtered_R2;

	}

	# ------------------------------

	# End subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	fastqc
# Description:	Handles qulity control of sequencing reads.

# --------------------------------------

sub fastqc {

	# Arguments:
	my ( $hash_config, $var_name_ID, $file_fastq_R1, $file_fastq_R2,
	$file_log, $var_threads ) = @_;

	# Data structures:

	# File paths:
	my $dir_out	= $hash_config->{subdirectories}->{fastqc};

	# Variables:
	my $file_out_R1	= ( split /\//, $file_fastq_R1 )[-1];
	my $file_out_R2 = ( split /\//, $file_fastq_R2 )[-1];
	$file_out_R1	= ( split /\./, $file_out_R1 )[0];
	$file_out_R2	= ( split /\./, $file_out_R2 )[0];
	$file_out_R1	= "$dir_out/$file_out_R1\_fastqc.zip";
	$file_out_R2	= "$dir_out/$file_out_R2\_fastqc.zip";
	my $var_return;


	# ------------------------------

	# Define command for fastqc:
	my $var_string1 = "fastqc \\\n";
	$var_string1 .= "	--outdir $dir_out \\\n";
	$var_string1 .= "	--threads $var_threads \\\n";
	$var_string1 .= "	$file_fastq_R1 $file_fastq_R2";

	# Check if $file_out is present and has non-zero size:
	if ( -s $file_out_R1 && -s $file_out_R2 ) {

		Comms::file_present ( $file_out_R1 );
		Comms::file_present ( $file_out_R2 );

		$var_return = 0;

	}

	else {

		# Check that $file_fastq_R1 and $file_fastq_R2 are present and
		# have non-zero sizes:
		unless ( -s $file_fastq_R1 && -s $file_fastq_R2 ) {

			Comms::file_missing ( $file_fastq_R1 ) unless -e $file_fastq_R1;
			Comms::file_missing ( $file_fastq_R2 ) unless -e $file_fastq_R2;

			$var_return = 1;

		}

		else {

			# Run command:
			$var_return = Wrapper::run ( $hash_config, $var_name_ID, "fastqc",
			"$var_string1", $dir_out, $file_log );

		}

	}

	# ------------------------------
	
	# Handle error exits codes by printing message to user:
	if ( $var_return != 0 ) { 

		Comms::exit_code ( $file_log );
		system ( "rm $file_out_R1" ) if -e $file_out_R1;
		system ( "rm $file_out_R2" ) if -e $file_out_R2;

	}

	# ------------------------------

	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:
# Description:

# --------------------------------------

sub dedupe {

	# Arguments:
	my () = @_;

	# Data structures:

	# File paths:

	# Variables:
	my $var_return;

	# ------------------------------

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	humann
# Description:	Performs metabolic profiling.

# --------------------------------------

sub humann {

	# Arguments:
	my ( $hash_config, $var_name_ID, $file_merged_fastq, $file_humann,
	$file_metaphlan, $file_log, $var_threads, $var_memory ) = @_;

	# Data structures:

	# File paths:
#	my $dir_logs		= $hash_config->{directories}->{logs};
	my $dir_out		= $hash_config->{subdirectories}->{humann};

	my $dir_metaphlan_db	= $hash_config->{databases}->{metaphlan};
	my $dir_humann_db	= $hash_config->{databases}->{humann};
	my $dir_chocophlan_db	= $hash_config->{databases}->{chocophlan};
	my $dir_uniref_db	= $hash_config->{databases}->{uniref};

	# Variables:
	my $var_humann		= $hash_config->{versions}->{humann};
	my $var_metaphlan	= $hash_config->{versions}->{metaphlan};
	my $var_return;

	# ------------------------------

	# Define command for running humann
	my $var_string1 = "humann \\\n";
	$var_string1 .= "	--input $file_merged_fastq \\\n";
	$var_string1 .= "	--output $dir_out \\\n";
	#$var_string1 .= "	--resume \\\n";
	$var_string1 .= "	--nucleotide-database $dir_chocophlan_db \\\n";
	$var_string1 .= "	--protein-database $dir_uniref_db \\\n";
	$var_string1 .= "	--metaphlan-options \\\"--bowtie2db $dir_metaphlan_db ";
	$var_string1 .= "--index $var_metaphlan \\\" \\\n";
	$var_string1 .= "	--remove-temp-output \\\n";
	$var_string1 .= "	--threads $var_threads \\\n";
	$var_string1 .= "	--memory-use $var_memory";

	# Define new command for running humann
	my $var_string2 = "humann \\\n";
	$var_string2 .= "	--input $file_merged_fastq \\\n";
	$var_string2 .= "	--output $dir_out \\\n";
	$var_string2 .= "	--threads $var_threads \\\n";
	$var_string2 .= "	--taxonomic-profile $file_metaphlan \\\n";
	$var_string2 .= "	--memory-use $var_memory \\\n";
	$var_string2 .= "	--search-mode uniref90 \\\n";
	$var_string2 .= "	--nucleotide-database $dir_chocophlan_db \\\n";
	$var_string2 .= "	--protein-database $dir_uniref_db \\\n";
	$var_string2 .= "	--taxonomic-profile $file_metaphlan \\\n";
	$var_string2 .= "	--remove-temp-output";

	# Check if $file_humann is present and has non-zero size:
	if ( -s $file_humann ) {

		Comms::file_present ( $file_humann );

		$var_return = 0;

	}

	else {

		# Check that $file_merged_fastq and $file_metaphlan are
		# present and has non-zero size:
		unless ( -s $file_merged_fastq & -s $file_metaphlan ) {

			Comms::file_missing ( $file_merged_fastq );
			Comms::file_missing ( $file_metaphlan );

			$var_return = 1;

		}

		else {

			# Run command:
			$var_return = Wrapper::run ( $hash_config, $var_name_ID, "humann",
			"$var_string2", $dir_out, $file_log );

		}

	}

	# ------------------------------
	
	# Handle error code from humann:
	if ( $var_return != 0 ) {

		Comms::exit_code ( $file_log );

		system ( "rm $file_humann" ) if -e $file_humann;

	}

	# ------------------------------

	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	kraken2
# Description:	Performs taxonomic profiling.

# --------------------------------------

sub kraken2 {

	# Arguments:
	my ( $hash_config, $var_name_ID, $file_filtered_R1, $file_filtered_R2,
	$file_kraken2_out, $file_kraken2_report, $file_log, $var_threads ) = @_;

	# Data structures:

	# File paths:
#	my $dir_logs		= $hash_config->{directories}->{logs};
	my $dir_out		= $hash_config->{subdirectories}->{kraken2};
	my $dir_kraken2_db	= $hash_config->{databases}->{kraken2};

	# Variables:
	my $var_return;

	# ------------------------------

	# Define command for running kraken2:
	my $var_string1 = "kraken2 \\\n";
	$var_string1 .= "	--db $dir_kraken2_db \\\n";
	$var_string1 .= "	--threads $var_threads \\\n";
	$var_string1 .= "	--output $file_kraken2_out \\\n";
	$var_string1 .= "	--report $file_kraken2_report \\\n";
	$var_string1 .= "	--paired \\\n";
	$var_string1 .= "	--use-names \\\n";
	$var_string1 .= "	--gzip-compressed \\\n";
	$var_string1 .= "	$file_filtered_R1 \\\n";
	$var_string1 .= "	$file_filtered_R2";

	# Check if $file_kraken2_out and $file_kraken2_report are present and
	# have non-zero size:
	if ( -s $file_kraken2_out && -s $file_kraken2_report ) {

		Comms::file_present ( $file_kraken2_out );
		Comms::file_present ( $file_kraken2_report );

		$var_return = 0;

	}

	else {

		# Check that $file_filtered_R1 and $file_filtered_R2 are present and
		# have non-zero sizes:
		unless ( -s $file_filtered_R1 && -s $file_filtered_R2 ) {

			Comms::file_missing ( $file_filtered_R1 ) unless -e $file_filtered_R1;
			Comms::file_missing ( $file_filtered_R2 ) unless -e $file_filtered_R2;

			$var_return = 1;

		}

		else {

			# Run command:
			$var_return = Wrapper::run ( $hash_config, $var_name_ID, "kraken2",
			"$var_string1", $dir_out, $file_log );

		}

	}

	# ------------------------------
	
	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	metaphlan
# Description:	Performs taxonomic profiling.

# --------------------------------------

sub metaphlan {

	# Arguments:
	my ( $hash_config, $var_name_ID, $file_filtered_R1, $file_filtered_R2,
	$file_bowtie, $file_metaphlan, $file_log, $var_threads ) = @_;

	# File paths:
	my $dir_logs		= $hash_config->{directories}->{logs};
	my $dir_out		= $hash_config->{subdirectories}->{metaphlan};
	my $dir_metaphlan_db	= $hash_config->{databases}->{metaphlan};

	# Variables:
	my $var_metaphlan	= $hash_config->{versions}->{metaphlan};
	my $var_return;

	# ------------------------------

	# Define command to run metaphlan:
	my $var_string1 = "metaphlan \\\n";

	# Handle novel run:
	unless ( -s "$file_bowtie" ) {

		$var_string1 .= "	$file_filtered_R1,$file_filtered_R2 \\\n";
		$var_string1 .= "	--input_type fastq \\\n";
		$var_string1 .= "	--bowtie2out $file_bowtie \\\n";

	}

	# Handle repeat run using bowte2 alignment:
	if ( -s "$file_bowtie" ) {

		$var_string1 .= "	$file_bowtie \\\n";
		$var_string1 .= "	--input_type bowtie2out \\\n";

	}

	# Add remaining flags:
	$var_string1 .= "	-t rel_ab_w_read_stats \\\n";
	$var_string1 .= "	--biom $file_metaphlan\.biom \\\n";
#	$var_string1 .= "	--tax_lev p \\\n";
	$var_string1 .= "	--index $var_metaphlan \\\n";
	$var_string1 .= "	--bowtie2db $dir_metaphlan_db \\\n";
	$var_string1 .= "	--nproc $var_threads \\\n";
	$var_string1 .= "	-o $file_metaphlan";

	# Check if $file_metaphlan is present and has non-zero size:
	if ( -s $file_metaphlan ) {

		Comms::file_present ( $file_metaphlan );

		$var_return = 0;

	}

	else {

		# Check that $file_filtered_R1 and $file_filtered_R2 are present and
		# have non-zero sizes:
		unless ( -s $file_filtered_R1 && -s $file_filtered_R2 ) {

			Comms::file_missing ( $file_filtered_R1 ) unless -e $file_filtered_R1;
			Comms::file_missing ( $file_filtered_R2 ) unless -e $file_filtered_R2;

			$var_return = 1;

		}

		else {

			# Run command:
			$var_return = Wrapper::run ( $hash_config, $var_name_ID, "metaphlan",
			"$var_string1", $dir_out, $file_log );

		}

	}

	# ------------------------------

	# Handle error code from metaphlan:
	if ( $var_return != 0 ) {

		Comms::exit_code ( $file_log );

		system ( "rm $file_metaphlan" ) if -e $file_metaphlan;

	}

	# ------------------------------
	
	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	quast
# Description:	Evaluate metagenomic assembly.

# --------------------------------------

sub quast {

	# Arguments:
	my ( $hash_config, $var_name_ID, $file_in, $file_quast, 
	$file_log, $var_threads ) = @_;

	# File paths:
	my $dir_out	= $hash_config->{subdirectories}->{quast};
	$dir_out	= "$dir_out/$var_name_ID";

	# Variables:
	my $var_return;

	# ------------------------------

	# Define command for running quast:
	my $var_string1 = "quast \\\n";
	$var_string1 .= "	--output-dir $dir_out \\\n";
	$var_string1 .= "	--threads $var_threads \\\n";
	$var_string1 .= "	--circos \\\n";
#	$var_string1 .= "	--max-ref-number 0 \\\n";
	$var_string1 .= "	$file_in";

	# Check if $file_quast is present and has non-zero size:
	if ( -s $file_quast ) {

		Comms::file_present ( $file_quast );

		$var_return = 0;

	}

	else {

		# Check that $file_in is present and has non-zero sizes:
		unless ( -s $file_in ) {

			Comms::file_missing ( $file_in );

			$var_return = 1;

		}

		else {

			# Run command:
			$var_return = Wrapper::run ( $hash_config, $var_name_ID, "quast",
			"$var_string1", $dir_out, $file_log );

		}

	}

	# ------------------------------
	
	# Return $var_eturn and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	prokka
# Description:	Performs annotation on metagenomic assemblies

# --------------------------------------

sub prokka {

	# Arguments:
	my ( $hash_config, $var_name_ID, $file_spades, $file_prokka, $file_log,
	$var_threads ) = @_;

	# Data structures:

	# File paths:
	my $dir_out	= $hash_config->{subdirectories}->{prokka};
	$dir_out	= "$dir_out/$var_name_ID";

	# Variables:
	my $var_return;

	# ------------------------------

	# Define command for running prokka:
	my $var_string1 = "prokka \\\n";
	$var_string1 .= "	--metagenome \\\n";
	$var_string1 .= "	--force \\\n";
	$var_string1 .= "	--compliant \\\n";
	$var_string1 .= "	--prefix $var_name_ID \\\n";
	$var_string1 .= "	--outdir $dir_out \\\n";
	$var_string1 .= "	--cpus $var_threads \\\n";
	$var_string1 .= "	$file_spades";

	# Check if $file_out is present and has non-zero size:
	if ( -s $file_prokka ) {

		Comms::file_present ( $file_prokka );

		$var_return = 0;

	}

	else {

		# Check that $file_spades is present and has non-zero sizes:
		unless ( -s $file_spades ) {

			Comms::file_missing ( $file_spades );

			$var_return = 1;

		}

		else {

			# Run command:
			$var_return = Wrapper::run ( $hash_config, $var_name_ID, "prokka",
			"$var_string1", $dir_out, $file_log );

		}

	}

	# ------------------------------
	
	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	rgi_load
# Description:	Loads RGI databases.

# --------------------------------------

sub rgi_load {

	# Arguments:
	my ( $hash_config, $var_database ) = @_;

	# Data structures:

	# File paths:
	my $dir_rgi_card_db	= $hash_config->{databases}->{rgi_card};
	my $dir_rgi_wild_db	= $hash_config->{databases}->{rgi_wildcard};

	# Variables:
	my $var_rgi_card	= $hash_config->{versions}->{rgi_card};
	my $var_rgi_wildcard	= $hash_config->{versions}->{rgi_wildcard};
	
	my $var_return;

	# ------------------------------

	# Define string to clean old databases:
	my $var_string1 = "rgi clean --local";

	# Define string for loading standard rgi database:
	my $var_string2 = "rgi load \\\n";
	$var_string2 .= "	-i $dir_rgi_card_db/card.json \\\n";
	$var_string2 .= "	--card_annotation $dir_rgi_card_db/card_database_v$var_rgi_card\.fasta \\\n";
	$var_string2 .= "	--local";

	# Define string for loading WildCARD rgi database:
	my $var_string3 = "rgi load \\\n";
	$var_string3 .= "	--wildcard_annotation $dir_rgi_wild_db";
	$var_string3 .= "/wildcard_database_v$var_rgi_wildcard\.fasta \\\n";
	$var_string3 .= "	--card_json $dir_rgi_card_db/card.json \\\n";
	$var_string3 .= "	--wildcard_index $dir_rgi_wild_db";
	$var_string3 .= "/index-for-model-sequences.txt \\\n";
	$var_string3 .= "	--card_annotation $dir_rgi_card_db";
	$var_string3 .= "/card_database_v$var_rgi_card\.fasta \\\n";
	$var_string3 .= "	--local";

	# ------------------------------

	# Remove old databases:
	Comms::print_command ( $var_string1 );
	system ( "$var_string1" );

	if ( $var_database eq "CARD" ) {

		Comms::print_command ( $var_string2 );

		# Load rgi database:
		$var_return = system ( "$var_string2" );

	}

	else {

		Comms::print_command ( $var_string3 );

		# Load rgi database:
		$var_return = system ( "$var_string3" );

	}

	# ------------------------------
	
	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	rgi_bwt
# Description:	Runs ARG prediction using RGI.

# --------------------------------------

sub rgi_bwt {

	# Arguments:
	my ( $hash_config, $var_name_ID, $file_fastq_R1, $file_fastq_R2,
	$file_rgi, $file_log, $var_threads ) = @_;

	# Data structures:

	# File paths:
	my $dir_out		= $hash_config->{subdirectories}->{rgi};
	$dir_out		= "$dir_out/$var_name_ID";

	# Variables:
	my $flag_wildcard	= $hash_config->{opts}->{flag_wildcard};
	my $var_return;

	# ------------------------------

	# Redefine $dir_out and $file_rgi for wildcard:
	$dir_out = "$dir_out\_wildcard" if $flag_wildcard;
	$file_rgi = "$file_rgi\_wildcard" if $flag_wildcard;

	# Define string for running rgi-bwt:
	my $var_string1 = "rgi bwt \\\n";
	$var_string1 .= "	-1 $file_fastq_R1 \\\n";
	$var_string1 .= "	-2 $file_fastq_R2 \\\n";
	$var_string1 .= "	--include_wildcard \\\n" if $flag_wildcard;
#	$var_string1 .= "	--include_other_models \\\n";
	$var_string1 .= "	-n $var_threads \\\n";
	$var_string1 .= "	-o $file_rgi \\\n";
	$var_string1 .= "	--local";

	# Add file ending to $file_rgi for skipping completed samples:
	$file_rgi .= ".allele_mapping_data.json";

	# Check if $file_out is present and has non-zero size:
	if ( -s $file_rgi ) {

		Comms::file_present ( $file_rgi );

		$var_return = 0;

	}

	else {

		# Check that $file_fastq_R1 and $file_fastq_R2 are present and
		# have non-zero sizes:
		unless ( -s $file_fastq_R1 && -s $file_fastq_R2 ) {

			Comms::file_missing ( $file_fastq_R1 ) unless -e $file_fastq_R1;
			Comms::file_missing ( $file_fastq_R2 ) unless -e $file_fastq_R2;

			$var_return = 1;

		}

		else {

			# Run command:
			$var_return = Wrapper::run ( $hash_config, $var_name_ID, "rgi",
			"$var_string1", $dir_out, $file_log );

		}

	}

	# Remove extra lines from log file:
	system ( "sed -i '/mapped query cannot have zero coordinate/d' $file_log" );

	# ------------------------------
	
	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Wrapper::seqtk
# Description:	Subsamples FASTQ reads

# --------------------------------------

sub seqtk {

	# Arguments:
	my ( $hash_config, $var_name_ID, $file_in, $file_out, $file_log,
	$var_threads ) = @_;

	# File paths:
	my $dir_logs	= $hash_config->{directories}->{logs};
	my $dir_out	= $hash_config->{subdirectories}->{seqtk};

	# Variables:
	my $var_sub	= $hash_config->{variables}->{sub};
	my $var_stats	= $hash_config->{variables}->{stats};
	my $var_return;

	# ------------------------------

	my $var_string1 = "seqtk sample \\\n";
	$var_string1 .= "	-s100 \\\n";
	$var_string1 .= "	$file_in \\\n";
	$var_string1 .= "	$var_sub";

	my $var_string2 = "pigz \\\n";
	$var_string2 .= "		--processes $var_threads \\\n";
	$var_string2 .= "		--fast \> $file_out";

	# Check if $file_out is present and has non-zero size:
	if ( -s $file_out ) {

		Comms::file_present ( $file_out );

		$var_return = 0;

	}

	else {

		# Check that $file_in is present and has non-zero size:
		unless ( -s $file_in ) {

			Comms::file_missing ( $file_in );

			$var_return = 1;

		}

		else {

			# Check that output directory is present:
			system ( "mkdir -p $dir_out" ) unless -d "$dir_out";

			# Check that logs directory is present:
			system ( "mkdir -p $dir_logs/$var_name_ID" ) unless -d "$dir_logs/$var_name_ID";

			# Print container information:
			Comms::container_info ( $hash_config, "seqtk", $file_log );

			# Print command:
			Comms::print_command ( "$var_string1 2>> $file_log | $var_string2", $file_log );

			# Run command:
			$var_return = system ( "$var_stats $var_string1 2>> $file_log | $var_string2" );

		}

	}

	# ------------------------------
	
	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	spades
# Description:	Perform metagenomic assembly.

# --------------------------------------

sub spades {

	# Arguments:
	my ( $hash_config, $var_name_ID, $file_filtered_R1, $file_filtered_R2,
	$file_spades, $file_log, $var_threads, $var_memory ) = @_;

	# Data structures:

	# File paths:
	my $dir_out	= $hash_config->{subdirectories}->{spades};

	# Variables:
	my $var_return;

	# ------------------------------

	# Define command for running spades:
	my $var_string1 = "spades \\\n";
	$var_string1 .= "	-o $dir_out/$var_name_ID \\\n";
	$var_string1 .= "	--meta \\\n";
	$var_string1 .= "	-1 $file_filtered_R1 \\\n";
	$var_string1 .= "	-2 $file_filtered_R2 \\\n";
	$var_string1 .= "	--threads $var_threads \\\n";
	$var_string1 .= "	--memory $var_memory";

	# Check if $file_out is present and has non-zero size:
	if ( -s $file_spades ) {

		Comms::file_present ( $file_spades );

		$var_return = 0;

	}

	else {

		# Check that $file_filtered_R1 and $file_filtered_R2 are present and
		# have non-zero sizes:
		unless ( -s $file_filtered_R1 && -s $file_filtered_R2 ) {

			Comms::file_missing ( $file_filtered_R1 ) unless -e $file_filtered_R1;
			Comms::file_missing ( $file_filtered_R2 ) unless -e $file_filtered_R2;

			$var_return = 1;

		}

		else {

			# Run command:
			$var_return = Wrapper::run ( $hash_config, $var_name_ID, "spades",
			"$var_string1", $dir_out, $file_log );

		}

	}

	# ------------------------------
	
	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	trimmomatic
# Description:	Trims reads below QC threshold or minimum length.

# --------------------------------------

sub trimmomatic {

	# Arguments:
	my ( $hash_config, $var_name_ID, $file_fastq_R1, $file_fastq_R2, 
	$file_trimmed_R1, $file_trimmed_R2, $file_unpaired_R1, $file_unpaired_R2,
	$file_log, $var_threads ) = @_;

	# Data structures:

	# File paths:
	my $dir_out		= $hash_config->{subdirectories}->{trimmomatic};	

	# Variables:
	my $var_adapter		= $hash_config->{variables}->{adapter};
	my $var_leading		= $hash_config->{variables}->{leading};
	my $var_trailing	= $hash_config->{variables}->{trailing};
	my $var_minlen		= $hash_config->{variables}->{minlen};
	my $var_avgqual		= $hash_config->{variables}->{avgqual};	
	my $var_return;

	# ------------------------------

	my $var_string1 = "trimmomatic PE \\\n";
	$var_string1 .= "	-threads $var_threads \\\n";
	$var_string1 .= "	-phred33 \\\n";
	$var_string1 .= "	$file_fastq_R1 \\\n";
	$var_string1 .= "	$file_fastq_R2 \\\n";
	$var_string1 .= "	$file_trimmed_R1 \\\n";
	$var_string1 .= "	$file_unpaired_R1 \\\n";
	$var_string1 .= "	$file_trimmed_R2 \\\n";
	$var_string1 .= "	$file_unpaired_R2 \\\n";
	$var_string1 .= "	ILLUMINACLIP:/Trimmomatic-0.39/\\";
	$var_string1 .= "adapters/$var_adapter\.fa:2:30:10 \\\n";
	$var_string1 .= "	LEADING:$var_leading \\\n";
	$var_string1 .= "	TRAILING:$var_trailing \\\n";
	$var_string1 .= "	MINLEN:$var_minlen \\\n";
	$var_string1 .= "	AVGQUAL:$var_avgqual";

	# Check if $file_out is present and has non-zero size:
	if ( -s $file_trimmed_R1 && -s $file_trimmed_R2 ) {

		Comms::file_present ( $file_trimmed_R1 );
		Comms::file_present ( $file_trimmed_R2 );

		$var_return = 0;

	}

	else {

		# Check that $file_fastq_R1 and $file_fastq_R2 are present and
		# have non-zero sizes:
		unless ( -s $file_fastq_R1 && -s $file_fastq_R2 ) {

			Comms::file_missing ( $file_fastq_R1 ) unless -e $file_fastq_R1;
			Comms::file_missing ( $file_fastq_R2 ) unless -e $file_fastq_R2;

			$var_return = 1;

		}

		else {

			# Run command:
			$var_return = Wrapper::run ( $hash_config, $var_name_ID, "trimmomatic",
			"$var_string1", $dir_out, $file_log );

		}

	}

	# ------------------------------

	# Handle error code from trimmomatic:
	if ( $var_return != 0 ) {

		Comms::exit_code ( $file_log );

		system ( "rm $file_trimmed_R1" ) if -e $file_trimmed_R1;
		system ( "rm $file_trimmed_R2" ) if -e $file_trimmed_R2;

	}

	# ------------------------------

	# End subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

1;
