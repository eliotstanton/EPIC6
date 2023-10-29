package Check;

use strict;
use warnings;
use Cwd 'abs_path';

# ---------------------------------------------------------------------------- #

# Module name:	Check
# Created by:	Eliot Stanton (estanton@wisc.edu)	
# Created on:	09 September, 2023
# Modified:	11 September, 2023
# Description:	Subroutines used for checking that directories, software, and
#		databases/indexes are all present before workflow launches.

# ---------------------------------------------------------------------------- #

# Subroutines:
# - Check::unify:	Runs the other subroutines to perform a full check out
# - Check::dir:		Checks if all or a specified directory is present
# - Check::command:	Checks if all or a specified software command is working
# - Check::compare_files Checks if list of files and files in directory are same
# - Check::container	Verifies if container is present and downloads/builds
# - Check::download_db	Downloads new database if required
# - Check::databases:	Checks if all or a specific database is present
# - Checl::genome:	Checks if all or a specific genome index is present
# - Check::setup_rgi	Handles RGI database setup
# - Check::wrapper:	Verifies if wrappers are present

# ---------------------------------------------------------------------------- #

# Subroutine:	Check::unify
# Description:	Runs all other subroutines in this module to perform a full
#		check.

# --------------------------------------

sub unify {

	# Arguments:
	my ( $hash_config ) = @_;

	# ------------------------------

	# Check that directories are present, if missing create them:
	dir ( $hash_config, "" );

	# Check that all needed commands are present, if missing download
	# containers and creater wrappers:
	command ( $hash_config, "" );

	# Check that all databases are present, if missing download and unpack
	# them:
	database ( $hash_config, "" );

	# Check that all genome indexes are present, if missing download and
	# index them:
	genome ( $hash_config, "" );

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Check::dir
# Description:	Checks if a specified directory is present, and creates it if
#		missing, defaults to checking all directories and
#		subdirectories.

# --------------------------------------

sub dir {

	# Arguments:
	my ( $hash_config, $var_dir ) = @_;

	# ------------------------------

	# If $var_dir is defined look up its path in $hash_config and if it is
	# present.
	if ( $var_dir ) {

		# Look up filepath for directory in %hash_config:
		my $dir_out = $hash_config->{directories}{ $var_dir };

		# Look up filepath for subdirectories if $dir_out is undefined:
		$dir_out = $hash_config->{subdirectories}{ $var_dir } unless $dir_out;

		# Look up filepath for databases if $dir_out is undefined:
		$dir_out = $hash_config->{databases}{ $var_dir } unless $dir_out;

		# Create directory if $dir_out doesn't exist:
		system ( "mkdir -p $dir_out" ) if $dir_out and ! -d $dir_out;

	}

	# Loop through directories checking and creating them as needed:
	else { 

		for my $dir_out ( values %{$hash_config->{directories}} ) {

			system ( "mkdir -p $dir_out" ) unless -d $dir_out;

		}

	}

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Check::commands
# Description:	Checks if specified command is returns exit code, downloads and
#		installs software if error, defaults to checking all software. 

# --------------------------------------

sub command {

	# Arguments:
	my ( $hash_config, $var_command ) = @_;

	# Data structures:

	# File paths:
	my $dir_wrapper = abs_path($hash_config->{directories}->{wrapper_home});

	# ------------------------------

	# Add $var_wrapper_home to path if absent:
	$ENV{PATH} .= ":$dir_wrapper" unless ( $ENV{PATH} =~ /$dir_wrapper:/ );

	# Loop through all commands:
	for my $var_software ( keys %{$hash_config->{software}} ) {

		next if $var_command && $var_command ne $var_software;

		container ( $hash_config, $var_software );

		wrapper ( $hash_config, $var_software );

	}

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Check::container
# Description:	Verifies that docker containers are present and downloads/builds
#		them as required.

# --------------------------------------

sub container {

	# Arguments:
	my ( $hash_config, $var_command ) = @_;

	# Data structures:
	my $hash_command	= $hash_config->{software}->{$var_command};

	# File paths:
	my $dir_container	= abs_path($hash_config->{directories}->{container_home});
	my $dir_def		= abs_path($hash_config->{subdirectories}->{def_files});

	# Variables:
	my $var_tag		= $hash_command->{tag};
	my $var_repo		= $hash_command->{repository};
	my $var_cmd		= $hash_command->{command};
	my $file_def		= $hash_command->{def_file};

	# ------------------------------

	# Define string for building container:
	my $var_string1 = "singularity build \\\n";
	$var_string1	.= "	--force \\\n";

	# Define container name:
	my $file_container	= "$var_command\_$var_tag.sif";

	# Handle requiring building the container using a deffile:
	if ( $file_def ) {

		$var_string1	.= "	--fakeroot \\\n";
#		$var_string1	.= "	--sandbox \\\n";
		$var_string1	.= "	$dir_container/$file_container \\\n";
		$var_string1	.= "	$dir_def/$file_def";

	}

	else {

		$var_string1	.= "	$dir_container/$file_container \\\n";
		$var_string1	.= "	docker://";
		$var_string1	.= "$var_repo/" if $var_repo;
		$var_string1	.= "$var_command:$var_tag";

	}

	# Download container if missing:
	unless ( -s "$dir_container/$file_container" ) {

		# Print command:
		Comms::print_command ( $var_string1 );

		# Run command:
		system ( "$var_string1" );

	}

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Database::compare_files
# Description:	This subroutine takes a file containing a list of files and
#		it with the files present in a directory to see if any are
#		missing.

# --------------------------------------

sub compare_files {

	# Arguments:
	my ( $dir_in, $file_in, $dir_check ) = @_;

	# Data structures:
	my @array_in; my %hash_in;
	my @array_check; my %hash_check;

	# Variables:
	my $var_return	= 0;

	# ------------------------------

	# Import $file_list into @array_list:
	@array_in	= @{ Edit::file_to_array ( "$dir_in/$file_in" ) };

	# Import contents of $dir_check to @array_check:
	@array_check	= @{ Edit::dir_to_array ( $dir_check ) };

	# Remove ./ from files in @array_in:
	foreach my $var_in ( @array_in ) { $var_in = substr $var_in, 2 if $var_in =~ /\.\//; };
	foreach my $var_in ( @array_in ) { $var_in = substr $var_in, 0, -3 if $var_in =~ /.gz$/; };

	# Move data from arrays to hashes:
	foreach my $var_in ( @array_in ) { next if $var_in eq ""; $hash_in{$var_in} = 1; };
	foreach my $var_check ( @array_check ) { $hash_check{$var_check} = 1; };

	# Compare contents of %hash_in and %hash_check by removing entries from 
	# %hash_in that are also present in %hash_check:
	foreach my $file_in ( sort keys %hash_in ) {

		delete $hash_in{$file_in} if $hash_check{$file_in};

	}

	# If files in %hash_in aren't accounted for set var_return to 1:
	$var_return = 1 unless keys %hash_in == 0;

	# ------------------------------

	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Check::database
# Description:	Checks if specified database is present, downloads and installs
#		the database if missing, defaults to checking all databases.

# --------------------------------------

sub database {

	# Arguments:
	my ( $hash_config, $var_database ) = @_;

	# Define filepaths:
	my $dir_ref		= $hash_config->{directories}->{ref};
	my $dir_chocophlan_db	= $hash_config->{databases}->{chocophlan};
	my $dir_humann_db	= $hash_config->{databases}->{humann};
	my $dir_kraken2_db	= $hash_config->{databases}->{kraken2};
	my $dir_metaphlan_db	= $hash_config->{databases}->{metaphlan};
	my $dir_rgi_card	= $hash_config->{databases}->{rgi_card};
	my $dir_rgi_wild	= $hash_config->{databases}->{rgi_wildcard};
	my $dir_uniref_db	= $hash_config->{databases}->{uniref};

	# Variables:
	my $var_humann		= $hash_config->{versions}->{humann};
	my $var_chocophlan	= $hash_config->{versions}->{chocophlan};
	my $var_kraken2		= $hash_config->{versions}->{kraken2};
	my $var_metaphlan	= $hash_config->{versions}->{metaphlan};
	my $var_rgi_card	= $hash_config->{versions}->{rgi_card};
	my $var_rgi_wildcard	= $hash_config->{versions}->{rgi_wildcard};
	my $var_uniref		= $hash_config->{versions}->{uniref};

	# Define URLs:
	my $var_chocophlan_URL	= $hash_config->{URLs}->{chocophlan};
	my $var_humann_URL	= $hash_config->{URLs}->{humann};
	my $var_kraken2_URL	= $hash_config->{URLs}->{kraken2};
	my $var_metaphlan_URL	= $hash_config->{URLs}->{metaphlan};
	my $URL_rgi_card	= $hash_config->{URLs}->{rgi_card};
	my $URL_rgi_wild	= $hash_config->{URLs}->{rgi_wildcard};
	my $var_uniref_URL	= $hash_config->{URLs}->{uniref};

	# ------------------------------

	# Loop through databases - if $var_db is defined only check for that
	# database otherwise check and install them all:
	for my $var_key ( sort keys %{$hash_config->{databases}} ) {

		# Unless $var_db is undefined, move on if $var_db doesn't match:
		next unless ( $var_database eq "" || $var_key eq $var_database );
		# Check that directory for database exists:
		Check::dir ( $hash_config, "$var_key" );

		# Check for Chocophlan database:
		if ( $var_key eq "chocophlan" ) {

			# Check database integrity and download if needed:
			download_db ( $dir_chocophlan_db, "$var_chocophlan\.list",
			$dir_chocophlan_db, "$var_chocophlan\.tar.gz", $var_chocophlan_URL );

			# Handle compressed genomic files:
			if ( my @files = glob ("$dir_chocophlan_db/*ffn.gz") ) {

				# Decompress files:
				system ( "pigz -d $dir_chocophlan_db/*ffn.gz" );

				# Amend list to reflect decompression of ffn files:
				system ( "tar -tf $dir_chocophlan_db/$var_chocophlan\.tar.gz | rev | cut -f 2- -d '.' | rev > $dir_chocophlan_db/$var_chocophlan\.list" );

			}

		}

		# Check for database used by Humann:
		elsif ( $var_key eq "humann" ) {

			download_db ( $dir_humann_db, "$var_humann\.list", $dir_humann_db,
			"$var_humann\.tar", $var_humann_URL );

			download_db ( $dir_humann_db, "$var_humann\_bt2.list", $dir_humann_db,
			"$var_humann\_bt2.tar", "$var_humann_URL/bowtie2_indexes" );

		}

		# Check for Kraken2 database:
		elsif ( $var_key eq "kraken2" ) {

			download_db ( $dir_kraken2_db, "$var_kraken2\.list", $dir_kraken2_db,
			"$var_kraken2\.tar.gz", $var_kraken2_URL );

		}

		# Check for Metaphlan database:
		elsif ( $var_key eq "metaphlan" ) {

			download_db ( $dir_metaphlan_db, "$var_metaphlan\.list", $dir_metaphlan_db,
			"$var_metaphlan\.tar", $var_metaphlan_URL );

			download_db ( $dir_metaphlan_db, "$var_metaphlan\_bt2.list", $dir_metaphlan_db,
			"$var_metaphlan\_bt2.tar", "$var_metaphlan_URL/bowtie2_indexes" );

		}
		
		# Check for RGI databases:
		elsif ( $var_key eq "rgi_card" ) {

			# Check that RGI is available:
			Check::command ( $hash_config, "rgi" );

			# Download CARD database:
			download_db ( $dir_rgi_card, "$var_rgi_card\.list", $dir_rgi_card, 
			"broadstreet-v$var_rgi_card\.tar.bz2", "$URL_rgi_card/0" );

			# Download WildCARD database:
			download_db ( $dir_rgi_wild, "$var_rgi_wildcard\.list",
			$dir_rgi_wild, "prevalence-v$var_rgi_wildcard\.tar.bz2", "$URL_rgi_wild/6" );

			# Set up RGI databases if needed:
			setup_rgi ( $hash_config );

		}

		# Check for Uniref databases:
		elsif ( $var_key eq "uniref" ) {

			download_db ( $dir_uniref_db, "$var_uniref\.list", $dir_uniref_db,
			"$var_uniref\.tar.gz", $var_uniref_URL );

		}

	}

	# ------------------------------
	
	# End subroutine:
	return;

}



# ---------------------------------------------------------------------------- #

# Subroutine:	download_db
# Description:	This subroutine handles verifying database presence/integrity
#		and downloads new files if required.

# --------------------------------------

sub download_db {

	# Arguments:
	my ( $dir_db, $file_list, $dir_out, $file_db, $var_URL ) = @_;

	# Variables:
	my $var_return	= 1;

	# ------------------------------

	# Define string to download files:
	my $var_string1 = "wget \\\n";
	$var_string1	.= "	--quiet \\\n";
	$var_string1	.= "	$var_URL";
	$var_string1	.= "/$file_db \\\n";
	$var_string1	.= "	-P $dir_db";

	# Define string to unarchive compressed tar file:
	my $var_string2 = "pigz -dc $dir_db/$file_db \\\n";
	$var_string2	.= "	| tar -xvC $dir_db";

	# Define string to unarchive tar file:
	my $var_string3 = "tar -xvf \\\n";
	$var_string3	.= "	$dir_db/$file_db \\\n";
	$var_string3	.= "	-C $dir_db";

	# ------------------------------

	# If list is present compare with the files present in directory
	if ( -s "$dir_db/$file_list" ) {

		# Compare contents of list and directory:
		$var_return = compare_files ( $dir_db, $file_list, $dir_out );

	}

	# If list is absent or files are missing ($var_return != 0):
	if ( ! -s "$dir_db/$file_list" || $var_return != 0 ) {

		# If database is present check database integrity:
		if ( -s "$dir_db/$file_db" ) {

			# If file is compressed:
			if ( $file_db =~ /\.tar\.gz$/ ) {

				$var_return = system ( "pigz -qdc $dir_db/$file_db | tar -t >> /dev/null" );

			}

			# If file is a regular tar archive:
			else {

				$var_return = system ( "tar -tf $dir_db/$file_db >> /dev/null" );

			}

			# If database is corrupted remove it:
			unlink "$dir_db/$file_db" if $var_return != 0;

		}

		# If database is absent download it:
		unless ( -s "$dir_db/$file_db" ) {

			Comms::print_command ( $var_string1 );
			system ( "$var_string1" );

		}

		# Uncompress or unarchive database and save list of files:
		if ( $file_db =~ /\.tar\.gz$/ ) {

			Comms::print_command ("$var_string2 > $dir_db/$file_list");
			system ( "$var_string2 > $dir_db/$file_list" );

		}

		else {

			Comms::print_command ("$var_string3 > $dir_db/$file_list");
			system ( "$var_string3 > $dir_db/$file_list" );

		}

	}

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Check::genome
# Description:	Checks if a specified genome index is present, downloads and 
#		indexes the genome file if missing, defaults to checking all
#		indexes.

# --------------------------------------

sub genome {

	# Arguments:
	my ( $hash_config, $var_index ) = @_;

	# Data structures:
	my $hash_genomes	= $hash_config->{genomes};

	# File paths:
	my $dir_ref		= $hash_config->{directories}->{ref};
	my $dir_bwa_index	= $hash_config->{subdirectories}->{bwa_index};

	# Variables:
	my $var_string		= "";

	# ------------------------------

	# Check that directories are present:
	Check::dir ( $hash_config, "ref" );
	Check::dir ( $hash_config, "bwa_index" );

	# Chck that BWA is available:
	Check::command ( $hash_config, "bwa" );

	# Iterate through genomes listed:
	for my $file_fasta ( values %{$hash_config->{genomes}} ) {

		# Add file to $var_string:
		$var_string	.= "$dir_ref/$file_fasta ";

		# Remove extension from $file_fasta:
		my $file_stem	= ( split /\./, $file_fasta)[0];

		# Define URL for downloading file:
		my $var_URL	= $hash_config->{URLs}->{$file_stem};

		# Define downloaded FASTA file:
		my $file_gzip	= ( split /\//, $var_URL)[-1];

		# Define string to hold command for downloading FASTA file:
		my $var_string1 = "wget $var_URL -P $dir_ref/temp";

		# Define string to hold command to uncompress file:
		my $var_string2 = "gunzip -f $dir_ref/temp/$file_gzip";

		 # Define string to hold command for moving and renaming hg38:
		my $var_string3 = "mv $dir_ref/temp/* $dir_ref/$file_fasta";

		# If $file_fasta is missing download it:
		if ( ! -e "$dir_ref/$file_fasta" ) {

			# Print message if $file_fasta is missing:
			Comms::file_missing ( $file_fasta ) unless -s "$dir_ref/$file_fasta";

			# Create temporary directory to hold FASTA file:
			system ( "mkdir -p $dir_ref/temp" );

			# Run command to download, uncompress, and move/rename
			# FASTA file:
			system ( "$var_string1" );
			system ( "$var_string2" );
			system ( "$var_string3" );

			# Remove temporary directory:
			system ( "rm -rf $dir_ref/temp" );

		}

	}

	# Generate concatenated FASTA file if it doesn't already exist:
	unless ( -e "$dir_ref/concatenated.fasta" ) {

		system ( "cat $var_string > $dir_ref/concatenated.fasta" );

	}

	# ------------------------------

	# Generate BWA index from concatenated hg38 and phiX file if it doesn't
	# already exist:
	unless ( -s "$dir_bwa_index/concatenated.fasta.amb"
		&& -s "$dir_bwa_index/concatenated.fasta.ann"
		&& -s "$dir_bwa_index/concatenated.fasta.bwt"
		&& -s "$dir_bwa_index/concatenated.fasta.pac"
		&& -s "$dir_bwa_index/concatenated.fasta.sa" ) {

		# Define string to hold bwa index command:
		my $var_string = "bwa index \\\n";
		$var_string .= "	-p $dir_bwa_index/concatenated.fasta \\\n";
		$var_string .= "	-b 400000000 \\\n";
		$var_string .= "	$dir_ref/concatenated.fasta";

		# Print command:
		Comms::print_command ( $var_string );

		# Run command:
		system ( "$var_string" );

	}

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Database::setup_rgi
# Description:	Handles setting up CARD and WildCARD databases.

# --------------------------------------

sub setup_rgi {

	# Arguments:
	my ( $hash_config ) = @_;

	# Filepaths:
	my $dir_rgi_card	= $hash_config->{databases}->{rgi_card};
	my $dir_rgi_wild	= $hash_config->{databases}->{rgi_wildcard}; 

	# Variables:
	my $var_rgi_card	= $hash_config->{versions}->{rgi_card};
	my $var_rgi_wildcard	= $hash_config->{versions}->{rgi_wildcard};

	# ------------------------------

	# Define string for command to download rgi_db card database:
	my $var_string1 = "rgi load \\\n";
	$var_string1	.= "	--card_json $dir_rgi_card/card.json \\\n";
	$var_string1	.= "	--local";

	# Define string for command to run CARD annotation:
	my $var_string2 = "rgi card_annotation \\\n";
	$var_string2	.= "	-i $dir_rgi_card/card.json && \\\n";
	$var_string2	.= "	mv card* $dir_rgi_card";

	# Define string for making directory for wildCARD:
	my $var_string4 = "mkdir -p $dir_rgi_wild";

	# Define string for decompressing wildCARD:
	my $var_string5 = "pigz -dfq $dir_rgi_wild/*.gz";

	# Define string for running wildCARD annotation:
	my $var_string6 = "rgi wildcard_annotation \\\n";
	$var_string6	.= "	-i $dir_rgi_wild \\\n";
	$var_string6	.= "	--card_json $dir_rgi_wild/card.json \\\n";
	$var_string6	.= "	-v $var_rgi_wildcard && \\\n";
	$var_string6	.= "	mv wildcard_database*.fasta $dir_rgi_wild";

	# ------------------------------

	# If CARD database files are missing:
	unless ( -s "$dir_rgi_card/card_database_v$var_rgi_card\_all.fasta"
		&& -s "$dir_rgi_card/card_database_v$var_rgi_card\.fasta" ) {

		# LOAD CARD:
		Comms::print_command ( "$var_string1" );
		system ( "$var_string1" );

		# Annotate CARD:
		Comms::print_command ( "$var_string2" );
		system ( "$var_string2" );

	}

	# Make directory for WILDCARD:
	system ( "$var_string4" );

	# If WildCARD database files are missing:
	unless ( -s "$dir_rgi_wild/wildcard_database_v$var_rgi_wildcard\_all.fasta"
	&& -s "$dir_rgi_wild/wildcard_database_v$var_rgi_wildcard\.fasta" ) {

		# Uncompress WILDCARD files:
		Comms::print_command ( "$var_string5" );
		system ( "$var_string5" );

		# Annotate WILDCARD:
		Comms::print_command ( "$var_string6" );
		system ( "$var_string6" );

	}

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Check::wrapper
# Description:	Looks for wrapper for command and creates one if required.

# --------------------------------------

sub wrapper {

	# Arguments:
	my ( $hash_config, $var_command ) = @_;

	# Data structures:
	my $hash_command	= $hash_config->{software}->{$var_command};
	my $hash_directories	= $hash_config->{directories};
	my $hash_resources	= $hash_config->{resources}->{slurm}->{$var_command};

	# File paths:
	my $dir_wrapper		= abs_path($hash_directories->{wrapper_home});
	my $dir_container	= abs_path($hash_directories->{container_home});

	# Variables:
	my $var_time		= $hash_resources->{time};
	my $var_threads		= $hash_resources->{cpus};
	my $var_memory		= $hash_resources->{memory};
	my $var_tag		= $hash_command->{tag};
	my $var_cmd		= $hash_command->{command};
	my $file_container	= "$var_command\_$var_tag.sif";
	use Cwd;
	my $var_cwd		= getcwd;

	# ------------------------------

	# Define string for wrapper:
	my $var_string1 .= "#!/usr/bin/bash\n\n";
	$var_string1 .= "#SBATCH --partition=local\n";
	$var_string1 .= "#SBATCH --time=\"$var_time\"\n" if $var_time;
	$var_string1 .= "#SBATCH --mem=\"$var_memory\"\n" if $var_memory;
	$var_string1 .= "CONTAINER_DIR=\"$dir_container\"\n";
	$var_string1 .= "IMAGE_NAME=\"$file_container\"";
	$var_string1 .= "\n\n";
	$var_string1 .= "CMD=$var_cmd\nARGS=\"\$@\"\n\n";
	$var_string1 .= "singularity run -B $var_cwd \$CONTAINER_DIR/";
	$var_string1 .= "\$IMAGE_NAME \$CMD \$ARGS";

	# Fix if wrapper is absent or not executable:
	unless ( -X "$dir_wrapper/$var_command" ) {

		# Write $var_string2 to $file_wrapper:
		Edit::write_string ( "$dir_wrapper/$var_command", $var_string1 );

		# Change permissions for wrapper:
		system ("chmod +rwx $dir_wrapper/$var_command" );

	}

	# ------------------------------
	
	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

1;
