# File name:		beta_div.R
# Date created:		15 September, 2023
# Date modified:	16 September, 2023
# Created by:		Eliot Stanton
# Description:		Script for calculating and plotting beta diversity.

# Define initial script setup --------------------------------------------------

# Source functions from other file:
source( "workflow/stats/functions2.R" )

# Define required packages:
if ( !require(ggplot2) ) { install.packages("ggplot2") }
#if ( !require(vegan) ) { install.packages("vegan") }

# Define arguments from command line:
args = commandArgs(trailingOnly=TRUE)

# Define file for analysis:
file_in <- args[2]

# Define metadata file:
file_metadata <- args[1]

# Define alpha file:
file_out <- args[3]

# Define beta plot file:
file_plot <- args [4]

# Import file_in as dataframe:
df_in <- read.table ( file_in, header = TRUE )

# Import file_metadata as dataframe:
df_metadata <- read.table ( file_metadata, header = TRUE )

# Format df_metadata:
df_metadata <- format_metadata ( df_metadata )

# Format df_in:
df_in <- format_df_in ( df_in )

# Define dataframe to hold beta diversity output:
df_beta <- df_metadata

# Define list to hold beta diversity output:
list_beta <- list ()

# Perform beta diversity calculations ------------------------------------------

# Calculate Bray-Curtis dissimilarity:
dist_bray <- beta_index ( df_in )

# Calculate dispersion and determine if variance differs by groups:
bray_dispersion <- beta_dispersion ( dist_bray, df_metadata )
bray_dispersion <- beta_permutest ( bray_dispersion )

# Perform ANOSIM:
df_anosim <- beta_anosim ( dist_bray, df_metadata )

# Define formula for PERMANOVA:
formula_new <- "dist_in ~ Daycare_attendance + Income"

# Perform PERMANOVA:
df_permanova <- permanova ( dist_bray, df_metadata, formula_new)

# Add output to list_beta:
list_beta <- list ( dist_bray, bray_dispersion, df_anosim, formula_new, df_permanova )

# Generate beta diversity scatter plot -----------------------------------------

# Create scatterplot:
plot_out <- beta_plot ( dist_bray, df_metadata, "Pastel1" )

# Save data to files -----------------------------------------------------------

# Print list_beta:
print (list_beta)

# Save list_beta to file_beta
save_text ( file_out, list_beta )

# Save plot_beta to file_plot_beta
save_pdf ( file_plot, plot_out, 7, 7 )
