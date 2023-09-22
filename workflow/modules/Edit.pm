package Edit;

use strict;
use warnings;

# ---------------------------------------------------------------------------- #

# Module name:	Edit
# Created by:	Eliot Stanton (estanton@wisc.edu)
# Created on:	08 September, 2023
# Modified:	12 September, 2023
# Description:	Subroutines used for exporting/importing data to/from external
#		files.

# ---------------------------------------------------------------------------- #

# Subroutines:
# - Edit::append_string		Add a string to the end of a file
# - Edit::dir_to_array:		Returns the contents of a directory as an array
# - Edit::keys_to_array:	Moves keys from has to an array
# - Edit::file_to_array:	Reads contents of a file into an array
# - Edit::trim_contigs		Trims contigs below established length
# - Edit::write_string:		Writes a string to file

# ---------------------------------------------------------------------------- #

# Subroutine:	Edit::append_string
# Description:	This subroutine appends a string to file.

# --------------------------------------

sub append_string {

	# Arguments:
	my ( $file_out, $var_string ) = @_;

	# Variables:
	my $file_write;

	# ------------------------------

	# Open $file_out:
	open ( $file_write, '>>', "$file_out" ) or die "Unable to open $file_out!";

	print $file_write $var_string;

	# Close $file_out:
	close $file_write;

	# ------------------------------

	# End subroutine:
	return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Edit::dir_to_array
# Description:	Imports the contents of a directory into an array

# --------------------------------------

sub dir_to_array {

	# Arguments:
	my ( $dir_in ) = @_;

	# Data structures:
	my @array_out;

	# ------------------------------

	# Open $dir_in to get list of files:
	opendir my $dir_read, "$dir_in" or die "Unable to open $dir_in!\n";

	# Iterate through files storing details in %hash_in and @array_in:
	while ( readdir $dir_read ) {

		# Skip placeholder entries:
		next if ( $_ eq ".." || $_ eq "." );

		push @array_out, $_;

	}

	# Close $dir_in:
	closedir $dir_read;

	# ------------------------------

	# return refeence to @array_out and subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Edit::file_to_array
# Description:	Converts a file into a array of strings

# --------------------------------------

sub file_to_array {

	# Arguments:
	my ( $file_in ) = @_;

	# Data structures:
	my @array_out;

	# ------------------------------

	# Open $file_in filehandle:
	open ( my $file_read, '<', "$file_in" ) or die "Unable to open $file_in!";

	while (<$file_read>) {

		my $var_line = $_;

		chomp ( $var_line );

		push @array_out, $var_line;

	}

	# Close $file_in filehandle:
	close $file_read;

	# ------------------------------

	# Return reference to @array_out:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Edit::keys_to_array
# Description:	Moves keys in hash to elements in a sorted array.

# --------------------------------------

sub keys_to_array {

	# Arguments:
	my ( $hash_in ) = @_;

	# Data structures:
	my @array_out;

	# ------------------------------

	# Move keys in %hash_in to @array_out:
	foreach my $var_name_ID ( keys %{$hash_in} ) {

		push @array_out, $var_name_ID;

	}

	# Sort @array_out:
	@array_out	= sort {$a cmp $b} @array_out;

	# ------------------------------

	# Return @array-out and end subroutine:
	return \@array_out;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Edit::trim_contigs
# Description:	This subroutine trims a genome assembly file based upon the
#		specified pre-determined minimum length

# --------------------------------------

sub trim_contigs {

	# Arguments:
	my ( $hash_config, $file_contigs, $file_spades, $var_min_contig,
	$var_name_ID ) = @_;

	# Variables:
	my $var_return;

	# ------------------------------

	# Add subject ID to beginning of contig headers and remove contigs that
	# are less than $var_min_contig_length:
	open ( my $file_read, '<', "$file_contigs" ) or die $!;
	open ( my $file_write, '>', "$file_spades" ) or die $!;

	# Open $file_out and iterate through each line:
	while ( <$file_read>) {

		# Define variable holding line from file:
		my $var_line = $_;

		# If $var_line is a FASTA header, split to find contig length:
		if ( $var_line =~ /^>/ ) {

			# Define contig length:
			my $var_contig_length = ( split ( /_/, $var_line ) )[3];

			# If contig length is below minimum end loop:
			last if ( $var_contig_length < $var_min_contig );# {

			# Otherwise revise header to remove > sign:
			$var_line		= ( split (/>/, $var_line ) )[1];

			# Create new header with sample ID:
			$var_line		= ">$var_name_ID\_$var_line";

		}

		# Print $var_line to $file_out.temp:
		print $file_write $var_line;

	}

	# Close file handles:
	close $file_read;
	close $file_write;

	# ------------------------------

	# Return $var_return and end subroutine:
	return $var_return;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Edit::write string
# Description:	Writes a string to file.

# --------------------------------------

sub write_string {

	# Arguments:
	my ( $file_out, $var_string ) = @_;

	# ------------------------------

	# Open $file_out:
	open ( my $file_write, '>', "$file_out" ) or die "Unable to open $file_out!";

	# Print string to file:
	print $file_write $var_string;

	# Close $file_out:
	close $file_write;

	# ------------------------------

	# End subroutine;
	return;

}

# ---------------------------------------------------------------------------- #

1;
