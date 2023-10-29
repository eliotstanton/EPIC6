#!/usr/bin/perl

# ---------------------------------------------------------------------------- #

use strict;
use warnings;

use lib './workflow/modules';
use Modules;

# ---------------------------------------------------------------------------- #

# Script name:  Run
# Created by:   Eliot Stanton (estanton@wisc.edu)
# Created on:   08 September, 2023
# Modified:     08 October, 2023
# Description:  Wrapper script used for trimming spades contigs when jobs
#		submitted to slurm.

# ---------------------------------------------------------------------------- #

my $hash_config;
my $file_contigs	= $ARGV[0];
my $file_spades		= $ARGV[1];
my $var_min_contig	= $ARGV[2];
my $var_name_ID		= $ARGV[3];

# Trim contigs file:
Edit::trim_contigs ( $hash_config, $file_contigs, $file_spades, $var_min_contig,
$var_name_ID );
