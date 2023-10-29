# File name:            lfc.R
# Date created:         09 October, 2023
# Date modified:        09 October, 2023
# Created by:           Eliot Stanton
# Description:          Script for visualizing log fold changes

# Define initial script setup --------------------------------------------------

# Source functions from other files:
source( "workflow/stats/functions2.R" )

# Define required packages:
if ( !require(ggplot2) ) { install.packages("ggplot2") }
if ( !require(dplyr) ) { install.packages("dplyr") }
#if ( !require(vegan) ) { install.packages("vegan") }
#if ( !require(stringr) ) { install.packages("stringr") }
#if (!require("BiocManager") ) { install.packages("BiocManager") }
#if ( !require("mia") ) { BiocManager::install("mia") }
#if ( !require("ANCOMBC") ) { BiocManager::install("ANCOMBC") }
#if ( !require("ALDEx2") ) { BiocManager::install("ALDEx2") }
#if ( !require("Maaslin2") ) { BiocManager::install("Maaslin2") }

# Define arguments from command line:
args = commandArgs(trailingOnly=TRUE)

# Define file for analysis:
file_in <- args[1]

# Define output file:
file_plot <- args[2]

# Import file_in as dataframe:
df_in <- read.table ( file_in, header = TRUE )

# Import file_metadata as dataframe:
#df_metadata <- read.table ( file_metadata, header = TRUE )

# Format df_metadata:
#df_metadata <- format_metadata ( df_metadata )

# Format df_in:
df_in <- format_df_in ( df_in )

# Define list to hold output:
list_out <- list ()

# Generate waterfall plot ------------------------------------------------------

# Sort df_in according to lfc:
df_in <- arrange( df_in, lfc_Daycare_attendance1 )

#print (df_in)

# Create vector based upon lfc sign:
df_sign <- factor ( sign(df_in$lfc_Daycare_attendance1), labels = c( "lfc Negative", "lfc Positive" ) )

# Create vector based upon sensitivity pass:
#df_passed <- replace ( df_in$passed_ss_Daycare_attendance1, df_in$passed_ss_Daycare_attendance1==FALSE, "red")
#df_passed <- replace ( df_passed, df_passed==TRUE, "black")

df_passed <- replace	( 
				df_in$passed_ss_Daycare_attendance1,
				df_in$passed_ss_Daycare_attendance1==TRUE & df_in$diff_Daycare_attendance1==TRUE,
				"black"
			)
df_passed <- replace ( df_passed, df_passed==FALSE | df_passed==TRUE, "red" )

#print (df_passed)

# Create new dataframe with important information:
df_out <- data.frame	(	df_in$taxon,
				df_in$lfc_Daycare_attendance1,
				df_in$se_Daycare_attendance1,
				df_passed,
				df_sign
			)

# Add column names to df_out:
colnames(df_out) <- c( "taxon", "lfc_Daycare", "se_Daycare", "ss_passed", "sign" )

# Sort df_out according to lfc:
df_out <- arrange( df_out, lfc_Daycare )

#print (df_out)

# Define plot:
plot_out <- ggplot ( df_out, aes (
					x = reorder ( taxon, -desc(lfc_Daycare)), 
					y = lfc_Daycare,
					fill = sign
				)
			) +
			geom_bar (
				stat = "identity",
				width = 0.7,
				color = "black",
				position = position_dodge(width = 0.4)
			) +
			geom_errorbar(
				aes(
					ymin = lfc_Daycare - se_Daycare,
					ymax = lfc_Daycare + se_Daycare
				), 
				width = 0.2,
				position = position_dodge(0.05),
				color = "black"
			) +
			labs(
				x = NULL,
				y = "Log fold change", 
				title = "Log fold changes as one unit increase of daycare attendance"
			) + 
			theme_classic() +
			theme (
				plot.title = element_text(hjust = 0.5),
				panel.grid.minor.y = element_blank(),
#				axis.text.x = element_blank()
				axis.text.x = element_text(
					angle = 90,
					hjust = 1,
					color = as.character ( df_passed )
				)
			) +
			scale_fill_brewer(
				palette = "Pastel1"
			) +
			coord_cartesian(ylim = c(-15, 15))


# Save data to file ------------------------------------------------------------

# Determine number of rows in df_out:
rows <- nrow( df_out )/4

# Save plot_out to file_plot:
save_pdf ( file_plot, plot_out, rows, 7 )
