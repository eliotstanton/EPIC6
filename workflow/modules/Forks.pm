package Forks;

use strict;

use warnings;

# ---------------------------------------------------------------------------- #

# Module name:	Forks
# Created by:	Eliot Stanton (estanton@wisc.edu)
# Created on:	13 September, 2023
# Modified:	13 September, 2023
# Description:	Contains subroutines that handle forking various stages.

# ---------------------------------------------------------------------------- #

# Subroutines:
# - fork_assemble:	Forks genome assembly and annotation
# - fork_classify:	Forks taxonomic and metabolic profile prediction	
# - fork_count_reads	Forks counting reads for multiple sequencing files.
# - fork_subsample:	Forks subsampling reads
# - fork_postprocess:	Forks QC of sequencing reads
# - fork_predict:	Forks ARG prediction
# - fork_process:	Forks read trimming and filtering	

# ---------------------------------------------------------------------------- #

# Subroutine:	fork_assemble
# Description:	Forks jobs calling Stage::assemble.

# --------------------------------------

sub fork_assemble {

	# Arguments:
	my ( $hash_config, $array_fastq ) = @_;

	# Data structures:
	my @array_fastq		= @$array_fastq;

	# Variables:
	my $var_assemble	= $hash_config->{forks}->{assemble};
	my $obj_assemble	= Parallel::ForkManager->new($var_assemble);

	# ------------------------------

	# Set $var_classify to number of jobs if larger:
	if ( scalar @array_fastq < $var_assemble ) {

		$var_assemble = scalar @array_fastq;

	}

	# Define threads to be used per job:
	my $var_threads	= $hash_config->{resources}->{threads};
	$var_threads	= int($var_threads/$var_assemble + 0.5);

	# Define memory used per job:
	my $var_memory	= $hash_config->{resources}->{memory};
	$var_memory	= int($var_memory/$var_assemble + 0.5);

	# ------------------------------

	# Iterate through @array_fastq assigning assemble jobs:
	foreach my $var_name_ID ( @array_fastq ) {

		# Start forking:
		$obj_assemble->start and next;

		# Call assemble subroutine:
		Stage::assemble ( $hash_config, $var_name_ID, $var_threads,
		$var_memory );

		# End forking:
		$obj_assemble->finish;

	}

	# Wait for all children from $obj_assemble:
	$obj_assemble->wait_all_children;

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	fork_classify
# Description:	Forks jobs calling Stage::classify.

# --------------------------------------

sub fork_classify {

	# Arguments:
	my ( $hash_config, $array_fastq ) = @_;

	# Data structures:
	my @array_fastq		= @$array_fastq;

	# Variables:
	my $var_classify	= $hash_config->{forks}->{classify};
	my $obj_classify	= Parallel::ForkManager->new($var_classify);

	# ------------------------------

	# Set $var_classify to number of jobs if larger:
	if ( scalar @array_fastq < $var_classify ) {

		$var_classify = scalar @array_fastq;

	}

	# Define threads to be used per job:
	my $var_threads	= $hash_config->{resources}->{threads};
	$var_threads	= int($var_threads/$var_classify + 0.5);

	my $var_memory	= "maximum";

	# ------------------------------

	# Iterate through @array_fastq assigning classify jobs:
	foreach my $var_name_ID ( @array_fastq ) {

		# Start forking:
		$obj_classify->start and next;

		# Call classify subroutine:
		Stage::classify ( $hash_config, $var_name_ID, $var_threads,
		$var_memory );

		# End forking:
		$obj_classify->finish;

	}

	# Wait for all children from $obj_classify:
	$obj_classify->wait_all_children;

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	fork_count_reads
# Description:	Handles counting original, trimmed, filtered, and deduplicated
#		read numbers.

# --------------------------------------

sub fork_count_reads {

	# Arguments:
	my ( $hash_config, $array_fastq ) = @_;

	# Data structures:
	my @array_fastq	= @$array_fastq;
	my $hash_files	= $hash_config->{filenames};
	my $hash_out;

	# File paths:
	my $file_reads		= $hash_config->{filenames}->{file_reads};

	# Variables:
	my $var_count_reads	= $hash_config->{forks}->{count_reads};
	my $obj_count_reads	= Parallel::ForkManager->new($var_count_reads);

	# ------------------------------

	# Set $var_count_reads to number of jobs if larger:
	if ( scalar @array_fastq < $var_count_reads ) {

		$var_count_reads = scalar @array_fastq;

	}

	# Define threads to be used per sample:
	my $var_threads	= $hash_config->{resources}->{threads};
	$var_threads	= int($var_threads/$var_count_reads + 0.5);

	# ------------------------------

	# Check that FASTQ files are in place:
	foreach my $var_name_ID ( @array_fastq ) {

		# Define file paths:
		my $file_fastq_R1	= $hash_files->{$var_name_ID}->{file_fastq_R1};
		my $file_fastq_R2	= $hash_files->{$var_name_ID}->{file_fastq_R2};
		my $file_trimmed_R1	= $hash_files->{$var_name_ID}->{file_trimmed_R1};
		my $file_trimmed_R2	= $hash_files->{$var_name_ID}->{file_trimmed_R2};
		my $file_filtered_R1	= $hash_files->{$var_name_ID}->{file_filtered_R1};
		my $file_filtered_R2	= $hash_files->{$var_name_ID}->{file_filtered_R2};

		unless ( -e $file_fastq_R1 && -e $file_fastq_R2 && -e $file_trimmed_R1 &&
		-e $file_trimmed_R2 && -e $file_filtered_R1 && -e $file_filtered_R2 ) {

			Comms::file_missing ( $file_fastq_R1 ) unless -e $file_fastq_R1;
			Comms::file_missing ( $file_fastq_R2 ) unless -e $file_fastq_R2;
			Comms::file_missing ( $file_trimmed_R1 ) unless -e $file_trimmed_R1;
			Comms::file_missing ( $file_trimmed_R2 ) unless -e $file_trimmed_R2;
			Comms::file_missing ( $file_filtered_R1 ) unless -e $file_filtered_R1;
			Comms::file_missing ( $file_filtered_R2 ) unless -e $file_filtered_R2;

			exit;

		}

	}

	# ------------------------------

	# If $file_reads already exists:
	if ( -s $file_reads ) {

		Comms::file_present ( $file_reads );

	}

	# Assess the number of reads and save in $file_reads:
	else {

		# Define subroutine for handling output from forks below:
		$obj_count_reads->run_on_finish (

			sub {

				my $var_name_ID = (keys %{$_[5]})[0];

				$hash_out->{$var_name_ID} = $_[5]{$var_name_ID};

			}

		);

		# Iterate through @array_fastq assigning postprocess jobs:
		foreach my $var_name_ID ( @array_fastq ) {

			# Start forking:
			$obj_count_reads->start and next;

			my $hash_num = Analysis::count_reads ( $hash_config, $var_name_ID, $var_threads );

			# End forking:
			$obj_count_reads->finish(0, $hash_num);

		}

		# Wait for all children from $obj_predict:
		$obj_count_reads->wait_all_children;

		# Print information to user and file:
		Analysis::print_reads ( $hash_config, $hash_out, $file_reads );

	}

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	fork_subsample
# Description:	Forks jobs calling Stage::subsample.

# --------------------------------------

sub fork_subsample {

	# Arguments:
	my ( $hash_config, $array_fastq ) = @_;

	# Data structures:
	my @array_fastq		= @$array_fastq;

	# Variables:
	my $var_subsample	= $hash_config->{forks}->{subsample};
	my $obj_subsample	= Parallel::ForkManager->new($var_subsample);

	# ------------------------------

	# Set $var_classify to number of jobs if larger:
	if ( scalar @array_fastq < $var_subsample ) {

		$var_subsample = scalar @array_fastq;

	}

	# Define threads to be used per job:
	my $var_threads	= $hash_config->{resources}->{threads};
	$var_threads	= int($var_threads/$var_subsample + 0.5);

	# ------------------------------

	# Iterate through @array_fastq assigning subsample jobs:
	foreach my $var_name_ID ( @array_fastq ) {

		# Start forking:
		$obj_subsample->start and next;

		# Call subsample subroutine:
		Stage::subsample ( $hash_config, $var_name_ID, $var_threads );

		# End forking:
		$obj_subsample->finish;

	}

	# Wait for all children from $obj_subsample:
	$obj_subsample->wait_all_children;

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	fork_postprocess
# Description:	Forks jobs calling Stage::postprocess.

# --------------------------------------

sub fork_postprocess {

	# Arguments:
	my ( $hash_config, $array_fastq ) = @_;

	# Data structures:
	my @array_fastq = @$array_fastq;

	# Variables:
	my $var_postprocess = $hash_config->{forks}->{postprocess};
	my $obj_postprocess = Parallel::ForkManager->new($var_postprocess);

	# ------------------------------

	# Set $var_classify to number of jobs if larger:
	if ( scalar @array_fastq < $var_postprocess ) {

		$var_postprocess = scalar @array_fastq;

	}

	# Define threads to be used per job:
	my $var_threads	= $hash_config->{resources}->{threads};
	$var_threads	= int($var_threads/$var_postprocess + 0.5);

	# ------------------------------

	# Check relevant commands:
	Check::command ( $hash_config, "fastqc" );

	# ------------------------------

	 # Iterate through @array_fastq assigning postprocess jobs:
	 foreach my $var_name_ID ( @array_fastq ) {

		# Start forking:
		$obj_postprocess->start and next;

		# Call postprocess subroutine:
		Stage::postprocess ( $hash_config, $var_name_ID, $var_threads );

		# End forking:
		$obj_postprocess->finish;

	}

	# Wait for all children from $obj_postprocess:
	$obj_postprocess->wait_all_children;

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	fork_predict
# Description:	Forks jobs calling Stage::predict.

# --------------------------------------

sub fork_predict {

	# Arguments:
	my ( $hash_config, $array_fastq ) = @_;

	# Data structures:
	my @array_fastq = @$array_fastq;

	# Variables:
	my $var_predict = $hash_config->{forks}->{predict};
	my $obj_predict = Parallel::ForkManager->new($var_predict);

	# ------------------------------

	# Set $var_classify to number of jobs if larger:
	if ( scalar @array_fastq < $var_predict ) {

		$var_predict = scalar @array_fastq;

	}

	# Define threads to be used per job:
	my $var_threads	= $hash_config->{resources}->{threads};
	print "$var_threads $var_predict\n";
	$var_threads	= int($var_threads/$var_predict + 0.5);

	# ------------------------------

	Check::command ( $hash_config, "rgi" );

	# ------------------------------

	# Load RGI database:
	Wrapper::rgi_load ( $hash_config, "wildcard" );

	# Iterate through @array_fastq assigning predict jobs:
	foreach my $var_name_ID ( @array_fastq ) {

		# Start forking:
		$obj_predict->start and next;

		# Call predict subroutine:
		Stage::predict ( $hash_config, $var_name_ID, $var_threads );

		# End forking:
		$obj_predict->finish;

	}

	# Wait for all children from $obj_predict:
	$obj_predict->wait_all_children;

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	fork_process
# Description:	Forks jobs calling Stage::process

# --------------------------------------

sub fork_process {

	# Arguments:
	my ( $hash_config, $array_fastq ) = @_;

	# Data structures:
	my @array_fastq	= @$array_fastq;

	# Variables:
	my $var_process	= $hash_config->{forks}->{process};
	my $obj_process	= Parallel::ForkManager->new($var_process);

	# ------------------------------

	# Set $var_classify to number of jobs if larger:
	if ( scalar @array_fastq < $var_process ) {

		$var_process = scalar @array_fastq;

	}

	# Define threads to be used per job:
	my $var_threads	= $hash_config->{resources}->{threads};
	$var_threads	= int($var_threads/$var_process + 0.5);

	# ------------------------------

	# Iterate through @array_fastq assigning process jobs:
	foreach my $var_name_ID ( @array_fastq ) {

		# Start forking:
		$obj_process->start and next;

		# Call process subroutine:
		Stage::process ( $hash_config, $var_name_ID, $var_threads );

		# End forking:
		$obj_process->finish;

	}

	# Wait for all children from $obj_process:
	$obj_process->wait_all_children;

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

1;
