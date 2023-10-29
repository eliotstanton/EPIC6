# File name:		abundance.R
# Date created:		15 September, 2023
# Date modified:	15 September, 2023
# Created by:		Eliot Stanton
# Description:		Script for calculating alpha and beta diversity from
#			a relative abundance table.

# Define initial script setup --------------------------------------------------

# Source functions from other files:
source( "workflow/stats/functions2.R" )

# Define required packages:
#if ( !require(vegan) ) { install.packages("vegan") }
#if ( !require(stringr) ) { install.packages("stringr") }
#if (!require("BiocManager") ) { install.packages("BiocManager") }
#if ( !require("mia") ) { BiocManager::install("mia") }
#if ( !require("ANCOMBC") ) { BiocManager::install("ANCOMBC") }
#if ( !require("ALDEx2") ) { BiocManager::install("ALDEx2") }
if ( !require("Maaslin2") ) { BiocManager::install("Maaslin2") }

# Define arguments from command line:
args = commandArgs(trailingOnly=TRUE)

# Define file for analysis:
file_in <- args[2]

# Define metadata file:
file_metadata <- args[1]

# Define output file:
file_out <- args[3]

# Define relative abundance file for Maaslin2:
file_relative <- args[4]

# Define output directory for Maaslin2:
dir_out <- args[5];

# Import file_in as dataframe:
df_in <- read.table ( file_in, header = TRUE )

# Import file_relative as dataframe:
df_relative <- read.table ( file_relative, header = TRUE )

# Import file_metadata as dataframe:
df_metadata <- read.table ( file_metadata, header = TRUE )

# Format df_metadata:
df_metadata <- format_metadata ( df_metadata )

# Format df_in:
df_in <- format_df_in ( df_in )

# Format df_relative:
df_relative <- format_df_in ( df_relative )

# Define list to hold output:
list_out <- list ()

# Perform differential abundance analysis --------------------------------------

# Define formula:
formula <- "Daycare_attendance + Income"

# Create phyloseq objects for original and stratified data:
physeq <- create_phyloseq ( df_in, df_metadata )

# Run differential analysis using ANCOMBC2:
ancombc2_out <- test_ancombc2 ( physeq, formula, "Species" )

# Run differential analysis using ALDeX2:
#aldex2_out <- test_aldex2 ( df_in, df_metadata, df_metadata$Daycare_attendance );

# Run differential analysis using Maaslin2:
maaslin2_out <- test_maaslin2 ( df_in, df_metadata, "Daycare_attendance,Income", dir_out )

# Add output to list_out
list_out <- list ( ancombc2_out ) #, aldex2_out )

# Save data to files -----------------------------------------------------------

# Print list_out:
#print (list_out)

#summary(warnings())

# Save list_out to file_out:
save_text ( file_out, list_out )
