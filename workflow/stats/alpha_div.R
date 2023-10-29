# File name:		alpha_div.R
# Date created:		15 September, 2023
# Date modified:	16 September, 2023
# Created by:		Eliot Stanton
# Description:		Script for calculating and plotting alpha diversity.

# Define initial script setup --------------------------------------------------

# Source functions from other file:
source( "workflow/stats/functions2.R" )

# Define required packages:
if ( !require(cowplot) ) { install.packages("cowplot") }
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

# Define dataframe to hold alpha diversity output:
df_alpha <- df_metadata

# Define list to hold alpha diversity output:
list_alpha <- list ()

# Define list to hold boxplots:
list_box_plots <- list ()

# Perform alpha diversity calculations -----------------------------------------

# Define list of indexes to be used:
list_indexes <- list( "chao1", "diversity_inverse_simpson", "diversity_shannon")

# Loop through list_indexes and calculate for each metric:
for ( i in list_indexes ) {

	# Calculate alpha index using df_in:
	df_index <- alpha_index ( df_in, i )

	# Merge df_index with df_alpha:
	df_alpha <- cbind( df_index, df_alpha)

}

# Loop through and perform testing:
for ( i in list_indexes ) {

	# Define formula for testing:
	formula <- paste ( i, "Daycare_attendance", sep = " ~ " )

	# Perform Wilcoxon testing:
	htest_out <- test_wilcox ( df_alpha, formula)

	# Create dataframe with output of interest:
	df_temp <-	data.frame(
				"formula" = c(formula),
				"statistic" = c(htest_out$statistic),
				"p-value" = c(htest_out$p.value)
			)

	# Store data in list_alpha:
	list_alpha <- list ( list_alpha, htest_out );

}

# Generate alpha diversity box plot --------------------------------------------

# Define variable for storing plots in list:
index <- 1

# Loop through list_index and create box plots:
for ( i in list_indexes ) {

	list_box_plots[[index]] <- box_plot ( df_alpha, "Daycare_attendance", i )

	index <- index + 1

}

# Plot boxplots for indexes:
plot_out <-	cowplot::plot_grid (
			list_box_plots[[1]] +
				coord_cartesian( ylim = c(1, 400) ),
			list_box_plots[[2]] +
				coord_cartesian( ylim = c(1, 40) ),
			list_box_plots[[3]] +
				coord_cartesian( ylim = c(1, 4.5) ),
			ncol = 1
		)

# Save data to files -----------------------------------------------------------

# Print list_alpha:
#print (list_alpha)

# Save list_alpha to file_alpha
save_text ( file_out, list_alpha )

# Save plot_alpha to file_plot_alpha
save_pdf ( file_plot, plot_out, 7, 7 )

