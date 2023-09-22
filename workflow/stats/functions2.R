# File name:		functions.R
# Date created:		15 September, 2023
# Date modified:	19 September, 2023
# Created by:		Eliot Stanton
# Description:		General functions used for analysis.

# ------------------------------------------------------------------------------

# Set seed for reproducibility:
set.seed(20210725)

# Functions for formatting data ------------------------------------------------

# Define function for formating metadata:
format_metadata <- function ( df_metadata ) {

	# Define rejected IDs:
	rejected_IDs <- c("E44_S5","E48_S8")

	# Remove row names with rejected IDs from metadata_df:
	df_metadata <- df_metadata[!(row.names(df_metadata) %in% rejected_IDs),]

	# Convert data in metadata to factors:
	regression_df <- data.frame ( lapply ( df_metadata, as.factor ), stringsAsFactors = TRUE )

	# Shift some data back to numeric:
	row.names ( regression_df ) <- row.names ( df_metadata )
	regression_df$Age <- as.numeric( df_metadata$Age )
	regression_df$Household_size <- as.numeric ( df_metadata$Household_size )
	regression_df$Income <- as.numeric ( df_metadata$Income )

	# Return formatted data frame:
	return ( regression_df )

}

# Define function for formatting df_in:
format_df_in <- function ( df_in ) {

	# Define rejected IDs:
	rejected_IDs <- c("E44_S5","E48_S8")

	# Remove columns with rejected IDs from df_in:
	df_in <- df_in [,!(names(df_in) %in% rejected_IDs)]

}

# Functions for saving data to file --------------------------------------------

# Define function for saving plots as PDFs:
save_pdf <- function ( file_out, plot_in, var_width, var_height) {

	# Define width and height:
	var_width <- as.integer(var_width)
	var_height <- as.integer(var_height)

	# Save pdf to file_out:
	pdf(
		file = file_out,
		width = var_width,
		height = 10,
		colormodel = "rgb",
		family = "ArialMT"
	)

	plot(plot_in)
	dev.off()

}

# Define function for saving data in text files:
save_text <- function (file_out, table_out ) {

	# Set maximum lines to print:
	options(max.print=999999)

	# Save text to file_out:
	sink(file = file_out)
	print( table_out, width = "1000")
	sink( file = NULL )

}

# Define alpha diversity functions ---------------------------------------------

# Define function wrapper for calculating alpha diversity:
alpha_index <- function ( df_in, index_name ) {

	# Calculate diversity:
	df_out <-	microbiome::alpha (
				df_in,
				index = index_name
			)

	# Return table_out:
	return ( df_out )

}

# Define function for plotting diversity as a boxplot:
box_plot <- function ( table_in, variable1, variable2 ) {

	# Redefine variables:
	variable1 <- variable1
	variable2 <- variable2
	variable1 <- ggplot2::ensym(variable1)
	variable2 <- ggplot2::ensym(variable2)

#	print (variable1)
#	print (variable2)

	plot_out <-	ggplot2::ggplot(
				table_in,
				aes(
					x=as.factor(!!variable1),
					y=!!variable2)
				) +
				geom_boxplot(
					outlier.colour = "BLUE",
					show.legend = FALSE
				)

	# Return plot_out:
	return (plot_out)

}

# Define function wrapper for performing Wilcoxon testing:
test_wilcox <- function ( df_in, formula ) {

	# Extract formula from variable: 
	formula <- as.formula(formula)

	# Perform Ranked Sum Test:
	object_out <-	wilcox.test(
				formula,
				data = df_in,
				exact = FALSE
			)

	# Return object_out:
	return ( object_out )

}

# Define beta diversity functions ----------------------------------------------

# Define function wrapper calculating Bray-Curtis dissimilarity:
beta_index <- function (df_in) {

	# Transpose the data for vegdist:
	df_in <- t(df_in)

	# Make distance matrix use Bray-Curtis method:
	dist_out <-	vegan::vegdist (
				df_in,
				method="bray",
				binary= FALSE,
				na.rm = TRUE
			)

	# Return dist_out:
	return ( dist_out )

}

# Define function wrapper for calculating dispersion:
beta_dispersion <- function ( dist_bray, df_metadata ) {

	df_out <-	vegan::betadisper (
				dist_bray,
				df_metadata$Daycare_attendance
			)

	# Return df_out:
	return ( df_out )

}

# Define function wrapper for calculating if dispersion varies by group:
beta_permutest <- function ( bray_dispersion ) {

	df_out <-	vegan::permutest (
				bray_dispersion,
				pairwise=TRUE,
				permutations=999
			)

	# Return df_out:
	return ( df_out )

}


# Define function wrapper for calculating if variation between groups exists:
beta_anosim <- function ( dist_bray, df_metadata ) {

	df_out <-	vegan::anosim (
				dist_bray,
				df_metadata$Daycare_attendance,
				permutations = 999
			)

	# Return df_out:
	return ( df_out )

}

# Define function wrapper for performing permanova:
permanova <- function ( dist_in, metadata, formula ) {

	# Format as formula:
	formula <- as.formula ( formula )

	# Perform testing:
	df_out <-	vegan::adonis2 (
				formula,
				data = metadata,
				na.action = na.exclude,
				permutations = 9999,
				method = bray,
				by = "margin"
			)

	# Return df_out:
	return ( df_out )

}

# Define function for generating beta diversity scatterplots:
beta_plot <- function ( df_in, metadata, palette ) {

	# Perform PCoA ordination:
	pcoa_out <- ape::pcoa( df_in )

	# Calculate percent explained variance:
	explain_1 <- ((pcoa_out$values[1])[1,])
	explain_2 <- ((pcoa_out$values[1])[2,])
	explain_sum <- (sum(pcoa_out$values[1]))
	percent_1 <- round( (explain_1 / explain_sum), digits = 4 )
	percent_2 <- round( (explain_2 / explain_sum), digits = 4 )

	percent_1 <- paste ( "pcoa1 % exp.", percent_1, sep = " " )
	percent_2 <- paste ( "pcoa2 % exp.", percent_2, sep = " " )

	# Convert to data-frame:
	df_pcoa <- data.frame( pcoa1 = pcoa_out$vectors[,1], pcoa2 = pcoa_out$vectors[,2] )

	# Add row names from metadata
	rownames(df_pcoa) <- row.names( metadata )

	# Merge df_out and metadata:
	df_pcoa <- data.frame ( df_pcoa, metadata )

	# Perform nMDS ordination:
#	df_nmds <- vegan::metaMDS ( df_in, distance="bray", k=2, trymax=1000 )

	# Convert to dataframe:
#	df_nmds <- as.data.frame ( df_nmds$points )

	# Merge with columns from metadata: 
#	df_nmds <- cbind ( df_nmds, metadata )

	# Convert Daycare to a character in dataframes:
	df_pcoa$Daycare_attendance <- as.factor(df_pcoa$Daycare_attendance)
#	df_nmds$Daycare_attendance <- as.character(df_nmds$Daycare_attendance)

	# Generate PCoA scatterplot:
	plot_out <-	ggplot( 
				df_pcoa,
				aes(
					x = pcoa1,
					y = pcoa2,
					color = Daycare_attendance
				)
				) +
			geom_point (
				size = 2
			) +
			scale_color_brewer(
				palette = palette
			) +
			theme_classic () +
			labs (
				x = percent_1,
				y= percent_2
			)

#	Generate nMDS scatterplot:
#	plot_nmds <- scatterplot ( df_nmds, "MDS1", "MDS2", "Daycare_attendance", palette)

	# Create legend:
#	legend <- get_legend(plot_pcoa)

	# Merge plots and legend:
#	plot_out <- cowplot::plot_grid(
#	plot_pcoa + theme(legend.position="none") + labs( x = percent_1, y= percent_2 ) + coord_cartesian(xlim = c($
#	plot_nmds + theme(legend.position="none"),
#	legend,
#	rel_widths = c(1,0.3),
#	ncol = 2 )

#	plot(plot_out)

	return ( plot_out )

}

# Helpers ----------------------------------------------------------------------


# Define function for creating phyloseq objects:
create_phyloseq <- function ( df_in, metadata ) {

	# Create tables for dummy OTUs and taxonomy:
	otumat <- df_in

	# Store rownames:
	df <- rownames ( df_in )

	# Split rownames by |:
	taxmat <- data.frame( stringr::str_split_fixed(df, "\\|", 7 ))

	# Move first column over to last if empty:
	if (length(unique(taxmat[,7]))==1) {

		taxmat[7] <- taxmat[1]
		taxmat[1] <- taxmat[2]

	}

	# Define column names for taxmat:
	colnames (taxmat) = c ( "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species" )

	# Remove row names from otumat:
	rownames(otumat) <- NULL

	# Convert otumat and taxmat to matrices:
	otumat <- data.matrix( otumat )
	taxmat <- as.matrix( taxmat )

	# Create phyloseq classes:
	OTU = phyloseq::otu_table( otumat, taxa_are_rows = TRUE )
	TAX = phyloseq::tax_table( taxmat )
	META = phyloseq::sample_data( metadata )

	# Add classes to phyloseq object:
	physeq_out = phyloseq::phyloseq( OTU, TAX, META )

	# Return phyloseq object:
	return ( physeq_out )

}

# Differential abundance analysis ----------------------------------------------

# Define function for performing DAA with ANCOMBC2:
test_ancombc2 <- function ( object_in, formula, tax_var ) {

	# Run ANCOMBC2 function:
	object_out <- ANCOMBC::ancombc2 (
				object_in,
				tax_level = tax_var,
				fix_formula = formula,
				p_adj_method = "fdr"
			)

#	print (object_out)

	# Return object_out:
	return (object_out)

}

# Define function for performing DAA with ALDex2:
test_aldex2 <- function ( df_in, df_metadata, variable ) {
 
df_out <-  ALDEx2::aldex(
                        df_in,
                        variable,
                        mc.samples=128,
                        test="t",
                        effect=TRUE,
                        include.sample.summary=TRUE,
                        verbose=FALSE,
                        paired.test=FALSE,
                        denom="all",
			useMC=true
               )
   
	# Define list of arguments for ALDEx2:
#	args =	list ( 
#			df_in, 
#			variable, 
#			mc.samples=1000, 
#			test="t", 
#			effect=TRUE, 
#			include.sample.summary=TRUE, 
#			verbose=FALSE,
#			paired.test=FALSE, 
#			denom="all"
#		)
    
	# Run ALDEx2:
#	df_out <- do.call ( ALDEx2::aldex, args)
    
#	print (df_out)
    
	# Return df_out:
	 return (df_out)
    
}

# Define function for performing DAA with Maaslin2:
test_maaslin2 <- function ( df_in, df_metadata, variable, output_dir ) {

	df_out <-	Maaslin2 ( 
				input_data = df_in,
				input_metadata = df_metadata,
				output = output_dir,
				fixed_effects = c("Daycare_attendance","Income"),
				normalization='NONE'
			)



	# Return df_out:
	return ( df_out )

}
