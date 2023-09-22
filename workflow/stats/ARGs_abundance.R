# File name:            ARGs_abundance.R
# Date created:         19 September, 2023
# Date modified:        19 September, 2023
# Created by:           Eliot Stanton
# Description:          Script for generating ARG abundance bar chart from a
#                       normalized abundance table.

# Define initial script setup --------------------------------------------------

# Source functions from other files:
source( "workflow/stats/functions2.R" )

# Define required packages:
#if ( !require(stringr) ) { install.packages("stringr") }
if ( !require(forcats) ) { install.packages("forcats") }
if ( !require(ggplot2) ) { install.packages("ggplot2") }
#if ( !require(vegan) ) { install.packages("vegan") }
#if ( !require(microshades) ) { remotes::install_github("KarstensLab/microshades") }
#if ( !require(speedyseq) ) { remotes::install_github("mikemc/speedyseq") }

# Define arguments from command line:
args = commandArgs(trailingOnly=TRUE)

# Define file for analysis:
file_in <- args[1]

# Define metadata file:
file_metadata <- args[2]

# Define output file:
file_plot <- args[3]

# Import file_in as dataframe:
df_in <- read.table ( file_in, header = TRUE )

# Import file_metadata as dataframe:
df_metadata <- read.table ( file_metadata, header = TRUE )

# Format df_metadata:
df_metadata <- format_metadata ( df_metadata )

# Format df_in:
df_in <- format_df_in ( df_in )

# Generate figure -------------------------------------------------------------

# Define formula used:
formula <- "Daycare_attendance"

# Create phyloseq object:
physeq <- create_phyloseq ( df_in, df_metadata )

#print ( (physeq@tax_table))

# Generate plot:
plot_out <- phyloseq::plot_bar( physeq, fill = "Species" )

# Create legend:
legend <- cowplot::get_legend ( plot_out )

# Place plots on grid:
plot_out <-	cowplot::plot_grid (
			plot_out +
				geom_bar(
					aes(
						fill = fct_reorder( 
							Species,
							Abundance,
							.desc = TRUE ), 
						color = Species
					), 
					stat = "identity",
					position="stack"
				) +
				theme(
					legend.position="none",
					axis.ticks.x=element_blank()
				) +
				facet_grid(
					formula,
					scales="free",
					labeller=label_both
				),
			legend,
			ncol = 1,
			rel_heights = c(4, 1)
			)

# Save data to file ------------------------------------------------------------

# Save plot_out to file_plot:
save_pdf ( file_plot, plot_out, 8.5, 11 )
