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

## Software

| Software | Link | Function | Compatible OS |
|---|---|---|---|
| SRA-Toolkit | https://github.com/ncbi/sra-tools/wiki | Retrieve data from SRA | UNIX/Unix-like, Win |
| SRA Run Selector | https://www.ncbi.nlm.nih.gov/Traces/study/ | Interactive filter and selection of SRA entries to obtain their metadata and accessions | Online |
| FastQC | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ | Quality control for the FASTQ files | UNIX/Unix-like, Win |
| Cutadapt | https://cutadapt.readthedocs.io/en/stable/index.html | Find and remove unwanted sequence from sequencing reads | UNIX/Unix-like, Win |
| STAR | https://github.com/alexdobin/STAR | RNA-seq read mapping | UNIX/Unix-like |
| kallisto | https://pachterlab.github.io/kallisto/ | RNA-seq read pseudomapping | UNIX/Unix-like |
| Samtools | http://www.htslib.org/ | View and manipulate SAM/BAM files | UNIX/Unix-like |
| RSEM | https://deweylab.github.io/RSEM/ | Expression quantification | UNIX/Unix-like |
| R | https://www.r-project.org/ | Commonly used programming language and analytical framework for statistics | UNIX/Unix-like, Win |
| DESeq2 | https://bioconductor.org/packages/release/bioc/html/DESeq2.html | Differential expression analysis | R package |


## RNA-seq Environment Setup

Create a dedicated Conda environment for RNA-seq preprocessing and analysis.

```bash
# Create environment
conda create -n rnaseq -c conda-forge -c bioconda \
    sra-tools \
    fastqc \
    cutadapt \
    star \
    kallisto \
    samtools \
    rsem \
    r-base \
    multiqc \
    -y
```

## Get the public RNA-seq data from SRA
The dataset contains **25 RNA-seq samples** obtained from the NCBI Sequence Read Archive (SRA).  
The corresponding SRA run accession numbers are stored in **SRR_Acc_List.txt**.

### Download the raw sequencing data in FASTQ format via SRA Toolkit

Create directories for the FASTQ files and temporary files:

```bash
mkdir -p rawdata rawdata/tmp
```

Download the SRA runs using `fasterq-dump`.  
A maximum of **2 samples are downloaded simultaneously**, with **8 threads per sample**.

```bash
while read SRR; do
    echo "Starting $SRR ..."

    mkdir -p "rawdata/tmp/$SRR"

    fasterq-dump "$SRR" \
        --split-files \
        --outdir rawdata \
        --temp "rawdata/tmp/$SRR" \
        --threads 8 &

    while [ "$(jobs -rp | wc -l)" -ge 2 ]; do
        sleep 5
    done

done < SRR_Acc_List.txt

wait
```

Compress all downloaded FASTQ files:

```bash
gzip rawdata/*.fastq
```

The final compressed FASTQ files will be stored in **rawdata** folder
```


