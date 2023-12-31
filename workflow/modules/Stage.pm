package Stage;

use strict;
use warnings;

# ---------------------------------------------------------------------------- #

# Module name:	Stage
# Created by:	Eliot Stanton (estanton@wisc.edu)
# Created on:	12 September, 2023
# Modified:	13 September, 2023
# Description:	Subroutines acting between wrappers and master control script.

# ---------------------------------------------------------------------------- #

# Subroutines:
# - assemble:		Perform assembly and annotation
# - classify:		Running taxonomic and metabolic profiling
# - predict:		Perform ARG prediction
# - process:		Trimming and filtering reads
# - postprocess:	Running QC on reads
# - subsample:		Subsampling reads using seqtk

# ---------------------------------------------------------------------------- #

# Subroutine:	assemble
# Description:	Runs SPAdes, Quast, and Prokka to assemble and annotate.

# --------------------------------------

sub assemble {

	# Arguments:
	my ( $hash_config, $var_name_ID, $var_threads, $var_memory, $var_jobID ) = @_;

	# Data structures:
	my $hash_filenames	= $hash_config->{filenames}->{$var_name_ID};

	# File paths:
	my $file_filtered_R1	= $hash_filenames->{file_filtered_R1};
	my $file_filtered_R2	= $hash_filenames->{file_filtered_R2};
	my $file_contigs	= $hash_filenames->{file_contigs};
	my $file_spades		= $hash_filenames->{file_spades};
	my $file_quast		= $hash_filenames->{file_quast};
	my $file_prokka		= $hash_filenames->{file_prokka};
	my $file_spades_log	= $hash_filenames->{file_spades_log};
	my $file_quast_log	= $hash_filenames->{file_spades_log};
	my $file_prokka_log	= $hash_filenames->{file_prokka_log};

	# Variables:
	my $flag_slurm		= $hash_config->{opts}->{flag_slurm};
	my $var_min_contig	= $hash_config->{variables}->{min_contig};
	my $var_announce;
	my $var_return;

	# ------------------------------

	# Check relevant commands:
	Check::command ( $hash_config, "spades" );
	Check::command ( $hash_config, "quast" );
	Check::command ( $hash_config, "prokka" );

	# ------------------------------

	Comms::start_finish ( "spades", $var_name_ID, "start", $file_spades_log );

	# Redefine resources if $flag_slurm is defined:
	$var_threads	= $hash_config->{resources}->{slurm}->{spades}->{cpus} if $flag_slurm;
	$var_memory	= $hash_config->{resources}->{slurm}->{spades}->{memory} if $flag_slurm;

	$var_return = Wrapper::spades ( $hash_config, $var_name_ID,
	$file_filtered_R1, $file_filtered_R2, $file_spades, $file_spades_log, 
	$var_threads, $var_memory, $var_jobID );

	# Print error message and die if $var_return doesn't equal zero:
	if ( $var_return != 0 && !$flag_slurm) { Comms::terminate_code ( "spades" ); die; }

	# Trim contigs file:
	Edit::trim_contigs ( $hash_config, $file_contigs, $file_spades,
	$var_min_contig, $var_name_ID ) unless $flag_slurm;

	# Define string to trim contigs file:
	my $var_string = "sbatch \\\n";
	$var_string .= "	--job-name=$var_name_ID\_trim_contigs \\\n";
	$var_string .= "	--dependency=afterok:$var_return \\\n";
	$var_string .= "	--wrap \"workflow/trim_contigs.pl \\\n";
	$var_string .= "		$file_contigs \\\n";
	$var_string .= "		$file_spades \\\n";
	$var_string .= "		$var_min_contig \\\n";
	$var_string .= "		$var_name_ID\"";

	Comms::print_command ( $var_string );

	# Trim contigs file submitting to Slurm:
	$var_return = `$var_string` if $flag_slurm;

	# Extract jobID:
	$var_jobID = (split " ", $var_return)[-1] if $flag_slurm;

	# Calculate and print the number of contigs:
	$var_announce = Comms::num_contigs ( $file_spades, $var_min_contig ) unless $flag_slurm;

	Comms::start_finish ( "spades", $var_name_ID, "finish", $file_spades_log );

	# ------------------------------

	Comms::start_finish ( "quast", $var_name_ID, "start", $file_quast_log );

	# Redefine resources if $flag_slurm is defined:
	$var_threads	= $hash_config->{resources}->{slurm}->{quast}->{cpus} if $flag_slurm;
	$var_memory	= $hash_config->{resources}->{slurm}->{quast}->{memory} if $flag_slurm;	

	$var_return	= Wrapper::quast ( $hash_config, $var_name_ID, 
	$file_spades, $file_quast, $file_quast_log, $var_threads, $var_jobID );

	# Print error message and die if $var_return doesn't equal zero:
	if ( $var_return != 0 && !$flag_slurm ) { Comms::terminate_code ( "quast" ); die; }

	Comms::start_finish ( "quast", $var_name_ID, "finish", $file_quast_log );

	# ------------------------------

	Comms::start_finish ( "prokka", $var_name_ID, "start", $file_prokka_log );

	# Redefine resources if $flag_slurm is defined:
	$var_threads	= $hash_config->{resources}->{slurm}->{prokka}->{cpus} if $flag_slurm;
	$var_memory	= $hash_config->{resources}->{slurm}->{prokka}->{memory} if $flag_slurm;

	$var_return	= Wrapper::prokka ( $hash_config, $var_name_ID, $file_spades,
	$file_prokka, $file_prokka_log, $var_threads, $var_jobID);

	# Print error message and die if $var_return doesn't equal zero:

	if ( $var_return == 0 && !$flag_slurm) { Comms::terminate_code ( "prokka" ); die; }

	Comms::start_finish ( "prokka", $var_name_ID, "finish", $file_prokka_log );

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	classify
# Description:	Coordinates taxonomic and metabolic profiling.

# --------------------------------------

sub classify {

	# Arguments:
	my ( $hash_config, $var_name_ID, $var_threads, $var_memory, $var_jobID ) = @_;

	# Data structures:
	my $hash_filenames	= $hash_config->{filenames}->{$var_name_ID};

	# File paths:
	my $file_filtered_R1	= $hash_filenames->{file_filtered_R1};
	my $file_filtered_R2	= $hash_filenames->{file_filtered_R2};

	my $file_kraken2_out	= $hash_filenames->{file_kraken2_out};
	my $file_kraken2_report	= $hash_filenames->{file_kraken2_report};

	my $file_bowtie		= $hash_filenames->{file_bowtie};
	my $file_metaphlan	= $hash_filenames->{file_metaphlan};

	my $file_merged_fastq	= $hash_filenames->{file_merged_fastq};
	my $file_humann		= $hash_filenames->{file_humann_path};

	my $file_kraken2_log	= $hash_filenames->{file_kraken2_log};
	my $file_metaphlan_log	= $hash_filenames->{file_metaphlan_log};
	my $file_humann_log	= $hash_filenames->{file_humann_log};

	# Variables:
	my $flag_slurm		= $hash_config->{opts}->{flag_slurm};
	my $var_return;

	# ------------------------------

	# Check relevant commands:
	Check::command ( $hash_config, "kraken2" );
	Check::command ( $hash_config, "metaphlan" );
	Check::command ( $hash_config, "humann" );

	# Check relevant databases:
	Check::database ( $hash_config, "kraken2" );
	Check::database ( $hash_config, "metaphlan" );
	Check::database ( $hash_config, "chocophlan" );
	Check::database ( $hash_config, "uniref" );

	# ------------------------------

	Comms::start_finish ( "kraken2", $var_name_ID, "start", $file_kraken2_log );

	# Redefine resources if $flag_slurm is defined:
	$var_threads	= $hash_config->{resources}->{slurm}->{kraken2}->{cpus} if $flag_slurm;

	$var_return = Wrapper::kraken2 ( $hash_config, $var_name_ID, $file_filtered_R1, $file_filtered_R2,
	$file_kraken2_out, $file_kraken2_report, $file_kraken2_log, $var_threads, $var_jobID );

	# Print error message and die if $var_return doesn't equal zero:
	if ( $var_return != 0 && !$flag_slurm) { Comms::terminate_code ( "kraken2" ); die; }

	Comms::start_finish ( "kraken2", $var_name_ID, "finish", $file_kraken2_log );

	# ------------------------------

	Comms::start_finish ( "metaphlan", $var_name_ID, "start", $file_metaphlan_log );

	# Redefine resources if $flag_slurm is defined:
	$var_threads	= $hash_config->{resources}->{slurm}->{metaphlan}->{cpus} if $flag_slurm;

	$var_return = Wrapper::metaphlan ( $hash_config, $var_name_ID, $file_filtered_R1, $file_filtered_R2,
	$file_bowtie, $file_metaphlan, $file_metaphlan_log, $var_threads, $var_jobID );

	# Print error message and die if $var_return doesn't equal zero:
	if ( $var_return != 0 && !$flag_slurm ) { Comms::terminate_code ( "metaphlan" ); die; }

	Comms::start_finish ( "metaphlan", $var_name_ID, "finish", $file_metaphlan_log );

	# ------------------------------

	# If $file_merged_fastq already exists check that it's the correct
	# size:
	if ( -s $file_merged_fastq ) {

		# Determine sizes of $file_filtered_R1, $file_filtered_R2, and
		# $file_merged_fastq:
		my $file_size_R1	= ( stat $file_filtered_R1)[7];
		my $file_size_R2	= ( stat $file_filtered_R2)[7];
		my $file_size_merged	= ( stat $file_merged_fastq)[7];

		# Remove file if $file_size_merged if abridged:
		unless ( $file_size_R1 + $file_size_R2 == $file_size_merged ) {

			my $var_string1 = "rm $file_merged_fastq";

			Comms::print_command ( $var_string1 );

			system ( $var_string1 );

		}

	}

	# Concatenate $file_filtered_R1 and $file_filtered_R2:
	unless ( -e $file_merged_fastq ) {

		my $var_string2	= "cat \\\n";
		$var_string2	.= "	$file_filtered_R1 \\\n";
		$var_string2	.= "	$file_filtered_R2 > $file_merged_fastq";

		system ( $var_string2 ) unless $flag_slurm;

		my $var_string3 = "sbatch \\\n";
		$var_string3 .= "	--dependency=afterok:$var_jobID \\\n" if $var_jobID;
		$var_string3 .= "	--job-name=$var_name_ID\_cat \\\n";
		$var_string3 .= "	--wrap \"$var_string2\"";

		Comms::print_command ( $var_string3 ) if $flag_slurm;

		$var_return = `$var_string3` if $flag_slurm;

		# Extract jobID:
                $var_return = (split " ", $var_return)[-1] if $flag_slurm;

	}

	Comms::start_finish ( "humann", $var_name_ID, "start", $file_humann_log );

	# Redefine resources if $flag_slurm is defined:
	$var_threads	= $hash_config->{resources}->{slurm}->{humann}->{cpus} if $flag_slurm;

	$var_return = Wrapper::humann ( $hash_config, $var_name_ID,
	$file_merged_fastq, $file_humann, $file_metaphlan, $file_humann_log,
	$var_threads, $var_memory, $var_return );

	# Print error message and die if $var_return doesn't equal zero:
	if ( $var_return != 0 && !$flag_slurm ) { Comms::terminate_code ( "humann" ); die; }

	Comms::start_finish ( "humann", $var_name_ID, "finish", $file_humann_log );

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	predict
# Description:	Perform antibiotic resistance gene (ARG) prediction

# --------------------------------------

sub predict {

	# Arguments:
	my ( $hash_config, $var_name_ID, $var_threads, $var_jobID ) = @_;

	# Data structures:
	my $hash_filenames	= $hash_config->{filenames}->{$var_name_ID};

	# File paths:
	my $file_filtered_R1	= $hash_filenames->{file_filtered_R1};
	my $file_filtered_R2	= $hash_filenames->{file_filtered_R2};
	my $file_rgi		= $hash_filenames->{file_rgi};
	my $file_rgi_log	= $hash_filenames->{file_rgi_log};

	# Variables:
	my $flag_slurm		= $hash_config->{opts}->{flag_slurm};
	my $var_return;

	# ------------------------------

	# Check relevant commands:
	Check::command ( $hash_config, "rgi" );

	# ------------------------------

	Comms::start_finish ( "rgi bwt", $var_name_ID, "start", $file_rgi_log );

	# Redefine resources if $flag_slurm is defined:
	$var_threads	= $hash_config->{resources}->{slurm}->{rgi}->{cpus} if $flag_slurm;

	$var_return = Wrapper::rgi_bwt ( $hash_config, $var_name_ID, $file_filtered_R1,
	$file_filtered_R2, $file_rgi, $file_rgi_log, $var_threads, $var_jobID );

	# Print error message and die if $var_return doesn't equal zero:
	if ( $var_return != 0 && !$flag_slurm ) { Comms::terminate_code ( "rgi bwt" ); die; }

	Comms::start_finish ( "rgi bwt", $var_name_ID, "finish", $file_rgi_log );

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	process
# Description:	Handles trimming and fltered reads using trimmomatic and bwa.

# --------------------------------------

sub process {

	# Arguments:
	my ( $hash_config, $var_name_ID, $var_threads ) = @_;

	# Data structures:
	my $hash_filenames	= $hash_config->{filenames}->{$var_name_ID};

	# File paths:
	my $file_concatenated	= $hash_config->{filenames}->{file_concatenated};
	my $file_fastq_R1	= $hash_filenames->{file_fastq_R1};
	my $file_fastq_R2	= $hash_filenames->{file_fastq_R2};
	my $file_trimmed_R1	= $hash_filenames->{file_trimmed_R1};
	my $file_trimmed_R2	= $hash_filenames->{file_trimmed_R2};
	my $file_unpaired_R1	= $hash_filenames->{file_unpaired_R1};
	my $file_unpaired_R2	= $hash_filenames->{file_unpaired_R2};
	my $file_filtered_R1	= $hash_filenames->{file_filtered_R1};
	my $file_filtered_R2	= $hash_filenames->{file_filtered_R2};
	my $file_bam		= $hash_filenames->{file_bam};
	my $file_trimmed_log	= $hash_filenames->{file_trim_log};
	my $file_bwa_log	= $hash_filenames->{file_bwa_log};

	# Variables:
	my $flag_slurm		= $hash_config->{opts}->{flag_slurm};
	my $var_return;

	# ------------------------------

	# Check relevant commands:
	Check::command ( $hash_config, "trimmomatic" );
	Check::command ( $hash_config, "bwa" );

	# Check genome indexes:
	Check::genome ( $hash_config );

	# ------------------------------

	Comms::start_finish ( "trimmomatic", $var_name_ID, "start", $file_trimmed_log );

	# Redefine resources if $flag_slurm is defined:
	$var_threads	= $hash_config->{resources}->{slurm}->{trimmomatic}->{cpus} if $flag_slurm;

	# Trim reads:
	$var_return = Wrapper::trimmomatic ( $hash_config, $var_name_ID, $file_fastq_R1, $file_fastq_R2, 
		$file_trimmed_R1, $file_trimmed_R2, $file_unpaired_R1, $file_unpaired_R2, 
		$file_trimmed_log, $var_threads );

	# Store job ID for SLURM dependencies in $hash_config:
	$hash_config->{jobID}->{$var_name_ID}->{"trimmomatic"} = $var_return if $flag_slurm;

	# Print error message and die if $var_return doesn't equal zero:
	if ( $var_return != 0 && !$flag_slurm ) { Comms::terminate_code ( "trimmomatic" ); die; }		

	Comms::start_finish ( "trimmomatic", $var_name_ID, "finish", $file_trimmed_log );

	# ------------------------------

	Comms::start_finish ( "bwa", $var_name_ID, "start", $file_bwa_log );

	# Redefine resources if $flag_slurm is defined:
	$var_threads = $hash_config->{resources}->{slurm}->{bwa}->{cpus} if $flag_slurm;

	# Filter reads:
	$var_return = Wrapper::bwa ( $hash_config, $var_name_ID, $file_trimmed_R1, $file_trimmed_R2,
	$file_concatenated, $file_bam, $file_filtered_R1, $file_filtered_R2,
	$file_bwa_log, $var_threads, $var_return );

	# Store job ID for SLURM dependencies in $hash_config:
	$hash_config->{jobID}->{$var_name_ID}->{bwa} = $var_return if $flag_slurm;

	# Print error message and die if $var_return doesn't equal zero:
	if ( $var_return != 0 && !$flag_slurm ) { Comms::terminate_code ( "bwa" ); die; }

	Comms::start_finish ( "bwa", $var_name_ID, "finish", $file_bwa_log );

	# ------------------------------
	
	# End subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	postprocess
# Description:	Manage running quality control on pre and post-processed
#		sequencing reads.

# --------------------------------------

sub postprocess {

	# Arguments:
	my ( $hash_config, $var_name_ID, $var_threads, $var_jobID ) = @_;

	# Data structures:
	my $hash_filenames	= $hash_config->{filenames}->{$var_name_ID};

	# File paths:
	my $dir_logs		= $hash_config->{directories}->{logs};
	my $file_fastq_R1	= $hash_filenames->{file_fastq_R1};
	my $file_fastq_R2	= $hash_filenames->{file_fastq_R2};
	my $file_filtered_R1	= $hash_filenames->{file_filtered_R1};
	my $file_filtered_R2	= $hash_filenames->{file_filtered_R2};
	my $file_log		= $hash_filenames->{file_fastqc_log};

	# Variables:
	my $flag_slurm		= $hash_config->{opts}->{flag_slurm};
	my $var_return;

	# ------------------------------

	# Check relevant commands:
	Check::command ( $hash_config, "fastqc" );

	# ------------------------------

	# Check that logs directory is present:
	system ( "mkdir -p $dir_logs/$var_name_ID" ) unless -d "$dir_logs/$var_name_ID";

	Comms::start_finish ("fastqc", $var_name_ID, "start", $file_log );

	# Redefine resources if $flag_slurm is defined:
	$var_threads = $hash_config->{resources}->{slurm}->{fastqc}->{cpus} if $flag_slurm;

	# Run QC on original reads:
	$var_return = Wrapper::fastqc ( $hash_config, $var_name_ID,
		$file_fastq_R1, $file_fastq_R2, $file_log, $var_threads );

	Comms::start_finish ("fastqc", $var_name_ID, "finish", $file_log );

	# ------------------------------

	Comms::start_finish ("fastqc", $var_name_ID, "start", $file_log );

	# Run QC on trimmed and filtered reads:
	$var_return = Wrapper::fastqc ( $hash_config, $var_name_ID,
	$file_filtered_R1, $file_filtered_R2, $file_log, $var_threads, $var_jobID );

	Comms::start_finish ("fastqc", $var_name_ID, "finish", $file_log );

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	subsample
# Description:	Manages subsampling pairs of FASTQ reads

# --------------------------------------

sub subsample {

	# Arguments:
	my ( $hash_config, $var_name_ID, $var_threads ) = @_;

	# Data structures:
	my $hash_filenames	= $hash_config->{filenames}->{$var_name_ID};

	# File paths:
	my $dir_logs		= $hash_config->{directories}->{logs};
	my $file_fastq_R1	= $hash_filenames->{file_fastq_R1};
	my $file_fastq_R2	= $hash_filenames->{file_fastq_R2};
	my $file_sub_R1		= $hash_filenames->{file_sub_R1};
	my $file_sub_R2		= $hash_filenames->{file_sub_R2};
	my $file_log		= $hash_filenames->{file_seqtk_log};

	# Variables:
	my $flag_slurm		= $hash_config->{opts}->{flag_slurm};
	my $var_return;

	# ------------------------------

	# Check relevant commands:
	Check::command ( $hash_config, "seqtk" );

	# ------------------------------

	# Check that logs directory is present:
	system ( "mkdir -p $dir_logs/$var_name_ID" ) unless -d "$dir_logs/$var_name_ID";

	Comms::start_finish ( "seqtk", $var_name_ID, "start", $file_log );

	# Redefine resources if $flag_slurm is defined:
	$var_threads = $hash_config->{resources}->{slurm}->{seqtk}->{cpus} if $flag_slurm;

	# Subsample R1 read:
	$var_return = Wrapper::seqtk ( $hash_config, $var_name_ID, $file_fastq_R1,
			$file_sub_R1, $file_log, $var_threads );

	# Print error message and die if $var_return doesn't equal zero:
	unless ( $var_return == 0 ) { Comms::terminate_code ( "seqtk" ); die; }

	# Subsample R2 read:
	$var_return = Wrapper::seqtk ( $hash_config, $var_name_ID, $file_fastq_R2,
			$file_sub_R2, $file_log, $var_threads );

	# Print error message and die if $var_return doesn't equal zero:
	unless ( $var_return == 0 ) { Comms::terminate_code ( "seqtk" ); die; }

	Comms::start_finish ( "seqtk", $var_name_ID, "finish", $file_log );

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

1;
