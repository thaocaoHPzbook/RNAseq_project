Table of Content
# 1.Preparation
## 1.1.Software setup for RNA-seq preprocessing and analysis
## 1.2.Get the public RNA-seq data from SRA

# 2.Preprocessing of RNA-seq data 
## 2.1.Quality control of RNA-seq data
## 2.2.Read mapping/pseudomapping and quantification
### 2.2.1.Read mapping with STAR and data quantification
### 2.2.2.Read pseudomapping with kallisto and data quantification
## 2.3.Cross-species comparison

# 3.Analyze and compare RNA-seq data
## 3.1.Import data to R
## 3.2.Comparison of transcriptomic profiles across samples
## 3.3.Differential expression analysis - DESeq2
## 3.4.Grouping of the identified DEGs
## 3.5.Making sense of the genes
## 3.6Other analysis



# 1.Preparation
## 1.1.Software setup for RNA-seq preprocessing and analysis
### 1.1.2. Install the required tools with the help from conda
Now you have access to the server/cluster and hopefully also know the basics of using it via the command line. The next step is to set up the tools required for the following data preprocessing and analysis.

Below is a summary of the main software that will be introduced and/or used throughout the workflow.

Software	Link	Function	Compatible OS
SRA Toolkit	GitHub	Retrieve sequencing data from the NCBI Sequence Read Archive (SRA)	UNIX/Unix-like, Windows
SRA Run Selector	NCBI SRA Run Selector	Interactive filtering and selection of SRA entries to obtain metadata and accession numbers	Online
FastQC	FastQC	Quality control of FASTQ files	UNIX/Unix-like, Windows
Cutadapt	Cutadapt Documentation	Find and remove unwanted sequences, such as adapters, from sequencing reads	UNIX/Unix-like, Windows
STAR	GitHub	RNA-seq read alignment to a reference genome	UNIX/Unix-like
kallisto	kallisto	RNA-seq transcript-level pseudoalignment and quantification	UNIX/Unix-like
Samtools	Samtools	View, process, and manipulate SAM/BAM files	UNIX/Unix-like
RSEM	RSEM	Gene and transcript expression quantification	UNIX/Unix-like
R	R Project	Programming language and statistical computing environment widely used for downstream analysis	UNIX/Unix-like, Windows
DESeq2	Bioconductor	Differential gene expression analysis	R package
