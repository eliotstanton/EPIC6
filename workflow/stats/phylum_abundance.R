# File name:		phylum_abundance.R
# Date created:		14 September, 2023
# Date modified:	16 September, 2023
# Created by:		Eliot Stanton
# Description:		Script for generating abundance bar chart from a
#			relative abundance table.

# Define initial script setup --------------------------------------------------

# Source functions from other files:
source( "workflow/stats/functions2.R" )

# Define required packages:
#if ( !require(stringr) ) { install.packages("stringr") }
if ( !require(ggplot2) ) { install.packages("ggplot2") }
#if ( !require(vegan) ) { install.packages("vegan") }
if ( !require(microshades) ) { remotes::install_github("KarstensLab/microshades") }
if ( !require(speedyseq) ) { remotes::install_github("mikemc/speedyseq") }

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

# Generate figure --------------------------------------------------------------

# Create phyloseq object:
physeq <- create_phyloseq ( df_in, df_metadata )

# Prepare phyloseq object for mapping colors:
mdf_prep <- prep_mdf( physeq ) #, subgroup_level = "Class" )

# Define phyla to be used:
groups <-	c(
			"p__Actinobacteria",
			"p__Bacteroidetes",
			"p__Firmicutes",
			"p__Proteobacteria",
			"p__Verrucomicrobia"
		)

# Define colors to be used on phyla:
colors <-	c(
			"micro_cvd_green",
			"micro_cvd_blue",
			"micro_cvd_orange",
			"micro_cvd_purple",
			"micro_cvd_turquoise"
		)

# Generate abundance sorted dataframe:
color_obj <-	create_color_dfs (
			mdf_prep,
			selected_groups = groups,
			top_n_subgroups = 4,
			group_level = "Phylum",
			subgroup_level = "Genus",
			cvd = TRUE
		)

# Extract data:
mdf_GP <- color_obj$mdf
cdf_GP <- color_obj$cdf

# Reassign colors to phyla:
cdf_GP <-	color_reassign ( 
			cdf_GP,
			group_assignment = groups,
			color_assignment = colors
		)

# Reorder samples based upon Firmicutes:
reordered <-	reorder_samples_by ( 
			mdf_GP,
			cdf_GP,
			order = "p__Firmicutes",
			group_level = "Phylum",
			subgroup_level = "Genus",
			sink_abundant_groups = TRUE
		)

	# Extract data again:
	mdf_GP_reordered <- reordered$mdf
	cdf_GP_reordered <- reordered$cdf

	# Generate initial plot:
	plot_1 <-	plot_microshades(
				mdf_GP_reordered,
				cdf_GP_reordered
			)

	# Generate legend:
	legend_1 <-	custom_legend(
				mdf_GP_reordered,
				cdf_GP_reordered,
				legend_key_size = 0.3,
				legend_text_size = 7
			)

	# Customize initial plot:
	plot_out_1 <-	plot_1 + 
				scale_y_continuous( 
					labels = scales::percent, 
					expand = expansion(0)
				) +
				theme ( 
					legend.key.size = unit(0.4, "cm"), 
					text = element_text(size=7),
					axis.text.x = element_text(size = 6),
					legend.position = "none"
				)

	# Place plots on grid:
	plot_out <-	cowplot::plot_grid(
				plot_out_1,
				legend_1,
				rel_widths = c(4, 1)
			)

# Save data to file ------------------------------------------------------------

# Save plot_out to file_plot:
save_pdf ( file_plot, plot_out, 7, 7 )
