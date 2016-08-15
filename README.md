shRNA_Screening_Analysis requries the following packages:
MASS, reshape2, and preprocessCore (a bioconductor package)
To use shRNA_Screening_Analysis simply save the file to a location which R has permission to access.
source the file in R 
# source("/path/to/shRNA_Screening_Analysis.r")
And run the function below with appropriate parameters.
# shRNA_Screening_Analysis()
Parameter options include:
# dat1 = "control.csv"	# Control data
# dat2 = "treat.csv"	# Treated data
# genes = "targets.txt"	# List of genes with targeting shRNAs
# BootRep = 10000		# Number of boot strap repetitions 
# thresh = 10			# Minimum number of shRNA reads to include in analysis. shRNAs with less than the threshold number of counts are converted to NA 
# passages = 5			# Number of passages (could be rounds of any sort of selection)
# bio.reps = 2			# Number of biological replicates
# adjust = "BH"			# Method to use for adjusting p-value for multiple testing. Argument is passed to p.adjust{stats}. see ?p.adjust for available options.
