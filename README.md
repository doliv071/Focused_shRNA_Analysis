# R Code for Focused shRNA Screening Analysis
# Authors
David Oliver, Piaomu Liu, Edsel Pe\~{n}a
### shRNA_Screening_Analysis.r requires:

+ MASS
+ reshape2
+ preprocessCore (bioconductor)

### shRNA_Screening_Analysis.r "install":

`git clone https://github.com/doliv071/Focused_shRNA_Analysis.git`

In R, source the file

`source(shRNA_Screening_Analysis.r)`

### shRNA_Screening_Analysis.r usage:

Run the function 

`shRNA_Screening_Analysis()`

required parameters:

+ dat1 = "control.csv"	# Control data
+ dat2 = "treat.csv"		# Treated data
+ genes = "targets.txt"	# List of genes with targeting shRNAs, 1 per line

**_\*\*NOTE\*\*:_** data files must be csv 

**_\*\*NOTE\*\*:_** rownames for the data must be gene-specific shRNA names that match the gene names specified by targets.txt. (e.g. For 3 shRNAs targeting RPL35A, the rownames should be RPL35A1, RPL35A2, RPL35A3) 

Additional parameters that can be specified:

+ BootRep = 10000		# Number of boot strap repetitions 
+ thresh = 10			# Minimum number of shRNA reads to include in analysis. shRNAs with less than the threshold number of counts are converted to NA 
+ passages = 5		# Number of passages (could be rounds of any sort of selection)
+ bio.reps = 2		# Number of biological replicates
+ adjust = "BH"		# Method to use for adjusting p-value for multiple testing. Argument is passed to p.adjust{stats}. see ?p.adjust for available options.

