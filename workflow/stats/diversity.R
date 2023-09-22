# File name:		diversity.R
# Date created:		14 July, 2023
# Date modified:	14 September, 2023
# Created by:		Eliot Stanton
# Description:		Script for calculating alpha and beta diversity from
#			a relative abundance table.

# Define initial script setup --------------------------------------------------

# Source functions from other files:
source( "workflow/stats/functions.R" )

# Set seed for reproducibility:
#set.seed(20210725)

# Define arguments from command line:
args = commandArgs(trailingOnly=TRUE)

# Define file for analysis:
file_in <- args[1]

# Define metadata file:
file_metadata <- args[2]

# Define alpha file:
file_alpha <- args[3]

# Define beta file:
file_beta <- args[4]

# Define alpha plot file:
file_plot_alpha <- args [5]

# Define beta plot file:
file_plot_beta <- args [6]

# Import file_in as dataframe:
df_in <- read.table ( file_in, header = TRUE )

# Import file_metadata as dataframe:
df_metadata <- read.table ( file_metadata, header = TRUE )

# Define rejected IDs:
rejected_IDs <- c("E44_S5","E48_S8")

# Remove columns with rejected IDs from df_in:
df_in <- df_in [,!(names(df_in) %in% rejected_IDs)]

# Remove row names with rejected IDs from metadata_df:
df_metadata <- df_metadata[!(row.names(df_metadata) %in% rejected_IDs),]

# Convert data in metadata to factors:
regression_df <- data.frame ( lapply ( df_metadata, as.factor ), stringsAsFactors = TRUE )

# Shift some data back to numeric:
row.names ( regression_df ) <- row.names ( df_metadata )
regression_df$Age <- as.numeric( df_metadata$Age )
regression_df$Household_size <- as.numeric ( df_metadata$Household_size )
regression_df$Income <- as.numeric ( df_metadata$Income )

# Rename back to df_metadata:
df_metadata <- regression_df
 
# Define dataframe to hold alpha diversity output:
df_alpha <- df_metadata

# Define list to hold alpha diversity output:
list_alpha <- list ()

# Define list to hold beta diversity output:
list_beta <- list ()

# Define list to hold box plots:
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

index <- 1

# Loop through list_index and create box plots:
for ( i in list_indexes ) {

	list_box_plots[[index]] <- box_plot ( df_in, "Daycare_attendance", i )

	index <- index + 1

}

# Plot boxplots for indexes:
plot_alpha <-	cowplot::plot_grid (
			list_box_plots[1] + 
				coord_cartesian( ylim = c(0, 350) ),
			list_box_plots[2] + 
				coord_cartesian( ylim = c(0, 35) ),
			list_box_plots[3] + 
				coord_cartesian( ylim = c(0, 5) ),
			ncol = 1
		)

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
plot_beta <- beta_plot ( dist_bray, df_metadata, "Pastel1" )

# Save data to files -----------------------------------------------------------

# Print list_alpha:
print (list_alpha)

# Save list_alpha to file_alpha
save_text ( file_alpha, list_alpha )

# Print list_beta:
print (list_beta)

# Save list_beta to file_beta
save_text ( file_beta, list_beta )

# Save plot_alpha to file_plot_alpha
save_pdf ( file_plot_alpha, plot_alpha, 7, 7 )

# Save plot_beta to file_plot_beta
save_pdf ( file_plot_beta, plot_beta, 7, 7 )
