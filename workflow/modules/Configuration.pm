package Configuration;

use strict;
use warnings;

# ---------------------------------------------------------------------------- #

# Module name:	Configuration
# Created by:	Eliot Stanton (estanton@wisc.edu)
# Created on:	08 September, 2023
# Modified:	19 September, 2023
# Description:	Subroutines used for initial configuration of workflow

# ---------------------------------------------------------------------------- #

# Subroutines:
# - Configuration::unify:	Coordinates other subroutines and returns
#				completed hash with all information needed.
# - Configuration::default:	Load default config keys and values into a hash
# - Configuration::fastq	Identify fastq files, check paired-end, return
#				hash with sequence IDs and files.
# - Configuration::filenames	Generate filenames for each sample and final
#				analysis.
# - Configuration::verify	Verify data in hash by checking for missing keys
#				and ensuring matching values.

# ---------------------------------------------------------------------------- #

# Subroutine:	Config::unify
# Description:	Runs subroutines for loading config file and adding fastq and
#		output files together to return completed hash_config.

# --------------------------------------

sub unify {

	my ( $hash_opts, $dir_in, $file_config ) = @_;

	# Data structures:
	my $hash_config;
	my $hash_default;
	my $hash_fastq;
	my $hash_filenames;

	# ------------------------------

	# Load defaults into $hash_default:
	$hash_default = default ();

	# If $file_config is present, import contents to $hash_config and
	# verify contents:
	if ( -s "$dir_in/$file_config" ) {

		# Load data from file:
		$hash_config = YAML::XS::LoadFile ( "$dir_in/$file_config" );

		# Verify data:
		$hash_config = verify ( $hash_config, $hash_default );

	}

	# If $file_config doesn't exist write one to $dir_in/$file_config using
	# and $hash_default and store as $hash_config:
	else {

		# Print error message to user:
		Comms::file_missing ( "$dir_in/$file_config" );

		# Convert $hash_default to a string:
		my $var_string = YAML::XS::Dump( $hash_default );

		# Write string to file:
		Edit::write_string ( "$dir_in/$file_config", $var_string );

		# Point $hash_default to $hash_config:
		$hash_config = $hash_default;

	}

	# Add $hash_opts to $hash_config:
	$hash_config->{opts} = $hash_opts;

	# Identify fastq files and add to $hash_config:
	$hash_config->{fastq} = fastq ( $hash_config );

	# Generate filenames for each sample and final analysis and add to
	# $hash_config:
	$hash_config->{filenames} = filenames ( $hash_config );

	# ------------------------------

	# Return $hash_config and end subroutine:
	return $hash_config;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Config::default
# Description:	Generates and returns a hash with default configuration
#		keys/values.

# -------------------------------------- 

sub default {

	# Data structures:
	my $hash_default;

	# ------------------------------

	$hash_default = {

	databases=>{
		chocophlan=>'ref/humann_db/chocophlan',
		humann=>'ref/humann_db/mpa_vJan21',
		kraken2=>'ref/kraken2_db',
		metaphlan=>'ref/metaphlan_db',
		rgi_card=>'ref/rgi_card',
		rgi_wildcard=>'ref/rgi_card/wildcard',
		uniref=>'ref/humann_db/uniref'
	},

	directories=>{
		analysis=>'analysis',
		container_home=>'singularity',
		fastq=>'fastq',
		logs=>'logs',
		ref=>'ref',
		work=>'work',
		wrapper_home=>'bin'
	},

	forks=>{
		assemble=>'4',
		classify=>'8',
		count_reads=>'8',
		predict=>'8',
		process=>'16',
		postprocess=>'48',
		subsample=>'48'
	},

	genomes=>{
		hg38=>'hg38.fa',
		phiX=>'phiX.fa'
	},

	resources=>{
		slurm=>{
			bwa=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			},
			fastqc=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			},
			humann=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			},
			kraken2=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			},
			metaphlan=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			},
			prokka=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			},
			quast=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			},
			'r-base'=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			},
			rgi=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			},
			samtools=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			},
			seqtk=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			},		
			spades=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			},
			trimmomatic=>{
				time=>'1:00:00',
				cpus=>'4',
				memory=>'4GB'
			}
		},
		threads=>'120',
		memory=>'896'
	},

	software=>{
		bwa=>{
			repository=>'staphb',
			container=>'bwa',
			tag=>'0.7.17',
			command=>'bwa'
		},
		fastqc=>{
			repository=>'staphb',
			container=>'fastqc',
			tag=>'0.11.9',
			command=>'fastqc'
		},
		humann=>{
			repository=>'biobakery',
			container=>'humann',
			tag=>'3.8',
			command=>'humann',
			def_file=>'humann_3.8.def'
		},
		kraken2=>{
			repository=>'staphb',
			container=>'kraken2',
			tag=>'2.1.2-no-db',
			command=>'kraken2'
		},
		metaphlan=>{
			repository=>'biobakery',
			container=>'metaphlan',
			tag=>'4.0.6',
			command=>'metaphlan',
			def_file=>'metaphlan_4.0.6.def'
		},
		prokka=>{
			repository=>'staphb',
			containter=>'prokka',
			tag=>'1.14.5',
			command=>'prokka'
		},
		'r-base'=>{
			repository=>'r-base',
			container=>'r-base',
			tag=>'4.3.1',
			command=>'Rscript',
			def_file=>'r-base_4.3.1.def'
		},
		quast=>{
			repository=>'staphb',
			container=>'quast',
			tag=>'5.0.2',
			command=>'quast.py'
		},
		rgi=>{
			repository=>'finlaymaguire',
			container=>'rgi',
			tag=>'latest',
			command=>'rgi'
		},
		samtools=>{
			repository=>'staphb',
			container=>'samtools',
			tag=>'1.16.1',
			command=>'samtools'
		},
		seqtk=>{
			repository=>'staphb',
			container=>'seqtk',
			tag=>'1.3',
			command=>'seqtk'
		},
		spades=>{
			repository=>'staphb',
			container=>'spades',
			tag=>'3.15.5',
			command=>'spades.py'
		},
		trimmomatic=>{
			repository=>'staphb',
			container=>'trimmomatic',
			tag=>'0.39',
			command=>'trimmomatic'
		}
	},

	subdirectories=>{
		amrfinder=>'work/amrfinder',
		bbtools_index=>'ref/bbtools',
		bwa=>'work/bwa',
		bwa_index=>'ref/bwa',
		dedupe=>'work/dedupe',
		def_files=>'workflow/def_files',
		fastqc=>'work/fastqc',
		humann=>'work/humann',
		kraken2=>'work/kraken2',
		metaphlan=>'work/metaphlan',
		prokka=>'work/prokka',
		quast=>'work/quast',
		resfinder=>'work/resfinder',
		rgi=>'work/rgi',
		seqtk=>'work/seqtk',
		shortbred=>'work/shortbred',
		spades=>'work/spades',
		trimmomatic=>'work/trimmomatic'
	},

	variables=>{
		adapter=>'NexteraPE-PE',
		avgqual=>'30',
		leading=>'20',
		min_contig=>'900',
		min_reads=>'100',
		minlen=>'150',
		stats=>'/usr/bin/time -f \Mem:%M\KB\ \CPU:%P\ \Time:%E',
		sub=>'0.02',
		threads=>'120',
		trailing=>'20'
	},

	URLs=>{
		chocophlan=>'https://huttenhower.sph.harvard.edu/humann_data/chocophlan',
		hg38=>'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz',
		humann=>'http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases',
		kraken2=>'https://genome-idx.s3.amazonaws.com/kraken',
		metaphlan=>'http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases',
		phiX=>'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Sinsheimervirus_phiX174/latest_assembly_versions/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz',
		rgi_card=>'https://card.mcmaster.ca/download',
		rgi_wildcard=>'https://card.mcmaster.ca/download',
		uniref=>'https://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated'

	},

	versions=>{
		chocophlan=>'full_chocophlan.v201901_v31',
		humann=>'mpa_vJan21_CHOCOPhlAnSGB_202103',
		kraken2=>'k2_standard_16gb_20230314',
		metaphlan=>'mpa_vOct22_CHOCOPhlAnSGB_202212',
		rgi_card=>'3.2.7',
		rgi_wildcard=>'4.0.1',
		uniref=>'uniref90_annotated_v201901b_full',
	}

	};

	# ------------------------------

	# Return $hash_default and end subroutine:
	return $hash_default;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Config::fastq
# Description:	Identifies paired-end fastq files in directory. Returns hash
#		with sequence IDs pointing to arrays with paired-end filenames.

# -------------------------------------- 

sub fastq {

	# Arguments:
	my ( $hash_config ) = @_;

	# Data structures:
	my $hash_fastq;
	my @array_fastq;

	# Filepaths:
	my $dir_fastq	= $hash_config->{directories}->{fastq};

	# Variables:
	my $var_ext1	= ".fastq.gz";
	my $var_ext2	= ".fastq";

	# ------------------------------

	# Check that $dir_fastq directory is present:
	Comms::dir_missing ( $dir_fastq ) and die unless -e $dir_fastq;

	# Import files in $dir_fastq into @array_fastq:
	@array_fastq	= @{ Edit::dir_to_array ( $dir_fastq )};

	# Iterate through @array_fastq and store FASTQ files in %hash_fastq:
	for ( my $i = 0; $i < scalar @array_fastq; $i++ ) {

		# Define FASTQ file:
		my $var_file	= $array_fastq[$i];

		# Remove from @array_fastq any files that don't end with
		# $var_ext1 or $var_ext2:
		unless ( $var_file =~ /$var_ext1$/ || $var_file =~ /$var_ext2$/) {

			# Remove element from array:
			splice @array_fastq, $i, 1;

			# Backup $i and advance loop:
			$i-- and next;

		}

		# Convert $var_file to an array splitting on underscores:
		my @array_file = split /_/, $var_file;

		# Iterate through @array_file determine sample name by locating
		# read number:
		for ( my $j = 0; $j < scalar @array_file; $j++ ) {

			# Define chunk:
			my $var_flag	= $array_file[$j];

			# If read number is recognized remove remainder of
			# elements in the array and end loop:
			if ( $var_flag =~ /R1/ || $var_flag =~ /R2/ ) {

				splice @array_file, $j;

				last;

			}

		}

		# Join sample name together:
		my $var_sample = join ( "_", @array_file );

		# Push filename to array refernce in %hash_fastq:
		push ( @{ $hash_fastq->{"$var_sample"} }, "$var_file" );

	}

	# ------------------------------

	# Check that all files in %hash_fastq are pairs:
	foreach my $var_sample ( sort keys %{$hash_fastq} ) {

		my @array_fastq = @{$hash_fastq->{ $var_sample }};

		unless ( scalar @array_fastq %2 == 0 ) {

			Comms::pair_missing ( $hash_fastq->{$var_sample}->[0] );

			delete $hash_fastq->{$var_sample};

		}

	}

	# Print warning if array is empty:
	Comms::dir_empty ( $dir_fastq ) if scalar $hash_fastq == 0;

	# ------------------------------

	# Return $hash_fastq and end subroutine:
	return $hash_fastq;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Config::filenames
# Description:	Generates output filenames for each sample and top-level 
#		analysis. Returns hash with filenames.

# -------------------------------------- 

sub filenames {

	# Arguments:
	my ( $hash_config ) = @_;

	# Data structures:
	my $hash_files;
	my $hash_fastq		= $hash_config->{fastq};
	my $hash_analysis;

	# File paths:
	my $dir_fastq		= $hash_config->{directories}->{fastq};
	my $dir_log		= $hash_config->{directories}->{logs};
	my $dir_analysis	= $hash_config->{directories}->{analysis};

	my $dir_bwa		= $hash_config->{subdirectories}->{bwa};
	my $dir_bwa_index	= $hash_config->{subdirectories}->{bwa_index};
	my $dir_dedupe		= $hash_config->{subdirectories}->{dedupe};
	my $dir_humann		= $hash_config->{subdirectories}->{humann};
	my $dir_kraken2		= $hash_config->{subdirectories}->{kraken2};
	my $dir_metaphlan	= $hash_config->{subdirectories}->{metaphlan};
	my $dir_prokka		= $hash_config->{subdirectories}->{prokka};
	my $dir_quast		= $hash_config->{subdirectories}->{quast};
	my $dir_rgi		= $hash_config->{subdirectories}->{rgi};
	my $dir_resfinder	= $hash_config->{subdirectories}->{resfinder};
	my $dir_seqtk		= $hash_config->{subdirectories}->{seqtk};
	my $dir_spades		= $hash_config->{subdirectories}->{spades};
	my $dir_trim		= $hash_config->{subdirectories}->{trimmomatic};

	# Variables:
	my ( $var_date, $var_time)		= Comms::date_time ();
	my $flag_wildcard			= $hash_config->{opts}->{flag_wildcard};

	# ------------------------------

	# Define concatenated FASTA index file:
	$hash_files->{file_concatenated}		= "$dir_bwa_index/concatenated.fasta";

	# Define read quantity output file:
	$hash_files->{file_reads}			= "$dir_analysis/surviving_reads.txt";

	# Define phylum abundance file:
	$hash_files->{file_phylum_abundance}		= "$dir_analysis/phylum_abundance.pdf";

	# Define phylogetic tree file:
	$hash_files->{file_phylo_pdf}			= "$dir_analysis/phylo_tree.pdf";

	# Define metadata file:
	$hash_files->{file_metadata}			= "metadata.txt";

	# ------------------------------

	# Iterate through each sample in $hash_fastq:
	foreach my $var_name_ID ( keys %{$hash_fastq} ) {

		my $hash_name_ID;

		# Define original FASTQ files:
		my $file_fastq_R1			= $hash_fastq->{$var_name_ID}->[0];
		my $file_fastq_R2			= $hash_fastq->{$var_name_ID}->[1];

		# Add original FASTQ files to $hash_name_ID:
		$hash_name_ID->{file_fastq_R1}		= "$dir_fastq/$file_fastq_R1";
		$hash_name_ID->{file_fastq_R2}		= "$dir_fastq/$file_fastq_R2";

		# Add subsampled FASTQ files to $hash_name_ID:
		$hash_name_ID->{file_sub_R1}		= "$dir_seqtk/$var_name_ID\_R1_sub.fastq.gz";
		$hash_name_ID->{file_sub_R2}		= "$dir_seqtk/$var_name_ID\_R2_sub.fastq.gz";

		# Add trimmed FASTQ files to $hash_name_ID:
		$hash_name_ID->{file_trimmed_R1}	= "$dir_trim/$var_name_ID\_R1_trimmed.fastq.gz";
		$hash_name_ID->{file_trimmed_R2}	= "$dir_trim/$var_name_ID\_R2_trimmed.fastq.gz";
		$hash_name_ID->{file_unpaired_R1}	= "$dir_trim/$var_name_ID\_R1_unpaired.fastq.gz";
		$hash_name_ID->{file_unpaired_R2}	= "$dir_trim/$var_name_ID\_R2_unpaired.fastq.gz";

		# Add filtered FASTQ files to $hash_name_ID:
		$hash_name_ID->{file_filtered_R1}	= "$dir_bwa/$var_name_ID\_R1_filtered.fastq.gz";
		$hash_name_ID->{file_filtered_R2}	= "$dir_bwa/$var_name_ID\_R2_filtered.fastq.gz";

		# Add filtered BAM file to $hash_name_ID:
		$hash_name_ID->{file_bam}		= "$dir_bwa/$var_name_ID\_filtered.bam";

		# Add deduplicated files to $hash_name_ID:
		$hash_name_ID->{file_dedupe}		= "$dir_dedupe/$var_name_ID\_dedupe.fastq.gz";
		$hash_name_ID->{file_dedupe_R1}		= "$dir_dedupe/$var_name_ID\_R1_dedupe.fastq.gz";
		$hash_name_ID->{file_dedupe_R2}		= "$dir_dedupe/$var_name_ID\_R2_dedupe.fastq.gz";

		# Add RGI and resfinder files to $hash_name_ID:
		$hash_name_ID->{file_rgi}		= "$dir_rgi/$var_name_ID/$var_name_ID\_rgi";
		$hash_name_ID->{file_rgi}		= "$dir_rgi/$var_name_ID\_wildcard/$var_name_ID\_rgi" if $flag_wildcard;

		$hash_name_ID->{file_allele}		= "$dir_rgi/$var_name_ID/$var_name_ID\_rgi.allele_mapping_data.txt";
		$hash_name_ID->{file_allele}		= "$dir_rgi/$var_name_ID\_wildcard/$var_name_ID\_rgi_wildcard.allele_mapping_data.txt" if $flag_wildcard;

		$hash_name_ID->{file_stats}		= "$dir_rgi/$var_name_ID/$var_name_ID\_rgi.overall_mapping_stats.txt";
		$hash_name_ID->{file_stats}		= "$dir_rgi/$var_name_ID\_wildcard/$var_name_ID\_rgi_wildcard.overall_mapping_stats.txt" if $flag_wildcard;

		$hash_name_ID->{file_resfinder}		= "$dir_resfinder/$var_name_ID/ResFinder_results.txt";

		# Add Metaphlan files to $hash_name_ID:
		$hash_name_ID->{file_bowtie}		= "$dir_metaphlan/$var_name_ID.bowtie2out.txt";
		$hash_name_ID->{file_metaphlan}		= "$dir_metaphlan/$var_name_ID\_metaphlan.tsv";

		# Add Humann files to $hash_name_ID:
		$hash_name_ID->{file_merged_fastq}	= "$dir_humann/$var_name_ID\_merged.fastq.gz";
		$hash_name_ID->{file_humann_gene}	= "$dir_humann/$var_name_ID\_merged_genefamilies.tsv";
		$hash_name_ID->{file_humann_path}	= "$dir_humann/$var_name_ID\_merged_pathabundance.tsv";
		$hash_name_ID->{file_humann_gene_cond}	= "$dir_humann/$var_name_ID\_merged_genefamilies_cond.tsv";
		$hash_name_ID->{file_humann_path_cond}	= "$dir_humann/$var_name_ID\_merged_pathabundance_cond.tsv";
		$hash_name_ID->{file_humann_regroup}	= "$dir_humann/$var_name_ID\_merged_genefamilies_regroup.tsv";
		$hash_name_ID->{file_humann_gene_relab}	= "$dir_humann/$var_name_ID\_merged_genefamilies_relab.tsv";
		$hash_name_ID->{file_humann_gene_cpm}	= "$dir_humann/$var_name_ID\_merged_genefamilies_cpm.tsv";
		$hash_name_ID->{file_humann_path_relab} = "$dir_humann/$var_name_ID\_merged_pathabundance_relab.tsv";
		$hash_name_ID->{file_humann_path_cpm}	= "$dir_humann/$var_name_ID\_merged_pathabundance_cpm.tsv";

		# Add Kraken2 files to $hash_name_ID:
		$hash_name_ID->{file_kraken2_out}	= "$dir_kraken2/$var_name_ID\_kraken2_out.txt";
		$hash_name_ID->{file_kraken2_report}	= "$dir_kraken2/$var_name_ID\_kraken2_report.txt";

		# Add SPAdes files to $hash_name_ID:
		$hash_name_ID->{file_contigs}		= "$dir_spades/$var_name_ID/contigs.fasta";
		$hash_name_ID->{file_spades}		= "$dir_spades/$var_name_ID/$var_name_ID\_contigs.fasta";

		# Add quast file to $hash_name_ID:
		$hash_name_ID->{file_quast}		= "$dir_quast/$var_name_ID/report.txt";

		# Add prokka files to $hash_name_ID:
		$hash_name_ID->{file_prokka}		= "$dir_prokka/$var_name_ID/$var_name_ID.gbk";

		# Add log file for seqtk to $hash_name_ID:
		$hash_name_ID->{file_seqtk_log}		= "$dir_log/$var_name_ID/$var_date\_seqtk.log";

		# Add log files for Trimmomatic:
		$hash_name_ID->{file_trim_log}		= "$dir_log/$var_name_ID/$var_date\_trimmomatic.log";

		# Add log file for BWA:
		$hash_name_ID->{file_bwa_log}		= "$dir_log/$var_name_ID/$var_date\_bwa.log";

		# Add log files for dedupe:
		$hash_name_ID->{file_dedupe_log}	= "$dir_log/$var_name_ID/$var_date\_dedupe.log";

		# Add log files for FastQC:
		$hash_name_ID->{file_fastqc_log}	= "$dir_log/$var_name_ID/$var_date\_fastqc.log";

		# Add log files for RGI:
		$hash_name_ID->{file_rgi_log}		= "$dir_log/$var_name_ID/$var_date\_rgi.log";

		# Add log files for SPAdes:
		$hash_name_ID->{file_spades_log}	= "$dir_log/$var_name_ID/$var_date\_spades.log";

		# Add log file for quast:
		$hash_name_ID->{file_quast_log}		= "$dir_log/$var_name_ID/$var_date\_quast.log";

		# Add log file for prokka:
		$hash_name_ID->{file_prokka_log}	= "$dir_log/$var_name_ID/$var_date\_prokka.log";

		# Add log files for metaphlann:
		$hash_name_ID->{file_metaphlan_log}	= "$dir_log/$var_name_ID/$var_date\_metaphlan.log";

		# Add log file for humann:
		$hash_name_ID->{file_humann_log}	= "$dir_log/$var_name_ID/$var_date\_humann.log";

		# Add log file for kraken2:
		$hash_name_ID->{file_kraken2_log}	= "$dir_log/$var_name_ID/$var_date\_kraken2.log";

		# Store $hash_name_ID in $hash_files:
		$hash_files->{$var_name_ID}		= $hash_name_ID;

	}

	# ------------------------------

	my @array_analysis = ( "metaphlan_tax", "rgi_ref_seq", "rgi_ARO_term", "rgi_AMR_fam",
	"rgi_drug_class", "rgi_res_mech", "humann_gene_families", "humann_path_abundance" );

	# Iterate through @array_analysis and generate output files:
	foreach my $var_term ( @array_analysis ) {

		#print "$var_term\n";

		# Define merged relative abundance output:
		$hash_analysis->{ $var_term }->{ "file_merged_rel" }	= "$dir_analysis/merged_$var_term\_rel.txt";

		# Define merged absolute abundance output:
		$hash_analysis->{ $var_term }->{ "file_merged_abs" }	= "$dir_analysis/merged_$var_term\_abs.txt" unless $var_term =~ /humann/;

		# Define merged normalized output:
		$hash_analysis->{ $var_term }->{ "file_merged_norm" }	= "$dir_analysis/merged_$var_term\_norm.txt" unless $var_term =~ /metaphlan_tax/;

		# Define merged normalized cpm output:
		$hash_analysis->{ $var_term }->{ "file_merged_cpm" }	= "$dir_analysis/merged_$var_term\_cpm.txt" if $var_term =~ /humann/;

		# Define alpha diversity output files:
		$hash_analysis->{ $var_term }->{ "file_alpha_div" }	= "$dir_analysis/alpha_div_$var_term\.txt";
		$hash_analysis->{ $var_term }->{ "file_alpha_div_pdf" }	= "$dir_analysis/alpha_div_$var_term\.pdf";

		# Define beta diversity output files:
		$hash_analysis->{ $var_term }->{ "file_beta_div" }	= "$dir_analysis/beta_div_$var_term\.txt";
		$hash_analysis->{ $var_term }->{ "file_beta_div_pdf" }	= "$dir_analysis/beta_div_$var_term\.pdf";

		# Define abundance files:
		$hash_analysis->{ $var_term }->{ "file_abundance_pdf" }	= "$dir_analysis/abundance_$var_term\.pdf";

		# Define DAA output files:
		$hash_analysis->{$var_term}->{ "file_DAA" }		= "$dir_analysis/DAA_$var_term\.tsv";

	}

	# Append wildcard to output files if $flag_wildcard is defined:
	if ( $flag_wildcard ) {

		# Iterate through terms in @array_analysis:
		foreach my $var_term ( @array_analysis ) {

			# Move on unless term involves RGI:
			next unless $var_term =~ /rgi/;

			# Append wildcard to filename:
			foreach my $var_tag ( keys %{$hash_analysis->{$var_term}} ) {

				my $var_file	= $hash_analysis->{$var_term}->{$var_tag};

				my @array_temp = ( split /\./, $var_file );

				$var_file = "$array_temp[0]\_wildcard.$array_temp[1]";

				$hash_analysis->{$var_term}->{$var_tag}	= $var_file;

			}			

		}

	}

	# Store $hash_analysis in $hash_files:
	$hash_files->{analysis} = $hash_analysis;

	# ------------------------------

	# Return $hash_files and end subroutine:
	return $hash_files;

}

# ---------------------------------------------------------------------------- #

# Subroutine:	Configuration::verify
# Description:	Verifies hash_config contains required values, adds defaults if
#		missing and checks for matching values as needed.

# --------------------------------------

sub verify {

	# Arguments:
	my ( $hash_config, $hash_default ) = @_;

	# ------------------------------

	# Iterate through keys in $hash_default and check if present in
	# $hash_config:
	foreach my $var_key1 ( keys %{$hash_default} ) {

		# Warn user if $var_key1 is absent and add data to $hash_config:
		unless ( $hash_config->{$var_key1} ) {

			# Print warning message:
			Error::config_missing ( $var_key1 );

			print "hi\n";

			# Add missing data:
			$hash_config->{$var_key1} = $hash_default->{$var_key1};

		}

		# Iterate through keys in each sub hash
		foreach my $var_key2 ( keys %{$hash_default->{$var_key1}} ) {

			# Warn user if $var_key2 is absent and add data to
			# $hash_config:
			unless ( $hash_config->{$var_key1}->{$var_key2} ) {

				# Print warning message:
				Comms::config_missing ( "$var_key1:$var_key2" );

				print "hi\n";

				# Add missing data:
				$hash_config->{$var_key1}->{$var_key2} = $hash_default->{$var_key1}->{$var_key2};

			}

		}

	}

	# ------------------------------

	# Entries in databases should have a match in URLs:
	foreach my $var_databases ( keys %{$hash_config->{databases} }) {

		# Warn user if $var_database is absent in URLs:
		unless ( $hash_config->{URLs}->{$var_databases} ) {

			# Print warning message:
			Comms::config_missing ( $var_databases );

		}

	}

	# Entries in genomes should have a match in URLs:
	foreach my $var_genomes ( keys %{$hash_config->{genomes} } ) {

		# Warn user if $var_genomes is absent in URLs:
		unless ( $hash_config->{URLs}->{$var_genomes} ) {

			# Print warning message:
			Comms::config_missing ( $var_genomes );

		}

	}

	# Entries in URLs should have a match in either databases or genomes: 
	foreach my $var_URLs ( keys %{$hash_config->{URLs} } ) {

		# Warn user if $var_URLs is absent in databases or genomes:
		unless ( $hash_config->{databases}->{$var_URLs} || $hash_config->{genomes}->{$var_URLs}) {

			# Print warning message:
			Comms::config_missing ( $var_URLs );

		}

	}

	# Entries in versions should have a match with databases:
	foreach my $var_versions ( keys %{$hash_config->{versions} } ) {

		# Warn user if $var_URLs is absent in databases or genomes:
		unless ( $hash_config->{versions}->{$var_versions} ) {

			# Print warning message:
			Comms::config_missing ( $var_versions );

		}

	}

	# Entries in software should have a match in resources->slurm:
	foreach my $var_software ( keys %{$hash_config->{software} } ) {

		# Warn user if $var_software is absent in slurm:
		unless ( $hash_config->{resources}->{slurm}->{$var_software} ) {

			# Print warning message:
			Comms::config_missing ( $var_software );

		}

	}

	# ------------------------------
	
	# End subroutine:
	return $hash_config;

}

# ---------------------------------------------------------------------------- #

1;
