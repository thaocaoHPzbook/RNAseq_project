<img width="1891" height="1112" alt="image" src="https://github.com/user-attachments/assets/4225334d-7b65-4c3c-b43f-0ea0ef2510c8" />Table of Content
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

## 1.2. Get the public RNA-seq data from SRA
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

gzip rawdata/*.fastq
```

The final compressed FASTQ files will be stored in **rawdata/** folder.


# 2. Preprocessing of RNA-seq data
## 2.1. Quality control of RNA-seq data
Before processing the RNA-seq data, the quality of the raw FASTQ files should be assessed.

### Run FastQC

Create directories for FastQC and MultiQC results:

```bash
mkdir -p qc/fastqc qc/multiqc
```

Run **FastQC** on all raw FASTQ files:

```bash
fastqc rawdata/*.fastq \
    -o qc/fastqc \
    -t 16
```

### Summarize QC results with MultiQC

Combine all FastQC reports into a single summary report:

```bash
multiqc qc/fastqc \
    -o qc/multiqc
```
At this stage, the main QC metrics to inspect are:

- **Per-base sequence quality / quality score**
<img width="1881" height="1125" alt="image" src="https://github.com/user-attachments/assets/e7ee5b98-0f5b-4b18-8132-1378bef207f6" />

- **Adapter content**
<img width="1912" height="883" alt="image" src="https://github.com/user-attachments/assets/157b74ad-9fc8-4ae6-a767-9b3662300ff2" />
- **GC content**
<img width="1891" height="1112" alt="image" src="https://github.com/user-attachments/assets/772a06e6-634c-4bfa-ab15-dd00e24fc88b" />

Other QC metrics can also be examined in more detail in the MultiQC report **qc/multiqc/multiqc_report.html**

It is important to mention that some of the grades are made assuming the data to be whole-genome DNA sequencing. For instance, the "Per sequence GC content" compares the distribution of G/C bases proportion per read to a theoretical distribution derived from the whole genome, which is expected to be different from the GC content of transcriptome. From my personal experience, the "Per base sequence content" and "Per sequence GC content" are the two sections that easily get the warning for failed grade for RNA-seq data, but can be ignored if other sections are fine. In addition, the "Sequence Duplication Levels" is another section that could give out warning of RNA-seq data, while it may or may not be a problem that needs to be solved later.

Meanwhile, the sections that I would suggest to pay attention to for RNA-seq data include "Per base sequence quality", "Sequence Duplication Levels", "Overrepresented sequences" and "Adapter Content". They represent potential We will need to try to fix the problem if they get a failed grade:

- If any read locus shows low quality (e.g. median <20 or even <10) in the "Per base sequence quality" section, especially at the two ends, we should try to trim them if the low-quality part is large (>10 bases), either by all reads removing the same number of bases or different number per read based on the quality scores.
- Since different transcripts have very different abundance, to make sure that very lowly expressed transcripts are also detected, it is possible that the highly expressed transcripts are over-amplified and/or over-sequenced, resulting in warning of "Sequence Duplication Levels". In this case, a de-duplication step may be wanted to collapse the identical reads into one.
- For standard mRNA-seq data with oligoT enrichment, problems of "Overrepresented sequences" and "Adapter Content" often come together and represent the adapter ligation issue mentioned above. We can try to cut the adapter sequences from reads later.

For the example shown in the screenshot above, we don't need to do anything as it looks all good.

IMPORTANT NOTE: It is not always necessary to do anything here even if problems were found, especially those related to base quality. For instance, many up-to-date software being used later for read mapping (e.g. STAR) has implemented a soft trimming mechanism to deal with low-quality bases at the end of a read.


## Adapter trimming

If adapter contamination is detected during the initial FastQC/MultiQC assessment, remove the adapter sequences before read alignment and quantification.

Create a directory for the trimmed reads:

```bash
cd /media/admin1/DATA/RNAseq_project

mkdir -p trimmed
```

Trim the Illumina adapter sequence using **Cutadapt**:

```bash
for file in rawdata/*.fastq; do
    sample=$(basename "$file" .fastq)

    echo "Trimming: $sample"

    cutadapt \
        --adapter=AGATCGGAAGAG \
        --minimum-length=25 \
        -j 4 \
        -o "trimmed/${sample}_trimmed.fastq.gz" \
        "$file"

done
```

Here:

- `--adapter=AGATCGGAAGAG` removes the Illumina adapter sequence from the 3' end of reads.
- `--minimum-length=25` removes reads shorter than 25 bp after trimming.
- `-j 4` uses 4 CPU threads.
- The trimmed reads are compressed automatically and stored in folder **trimmed/SRRxxxxxxx_trimmed.fastq.gz**

### Quality control after trimming

After adapter trimming, run FastQC again to confirm that adapter contamination has been removed and that the remaining reads retain acceptable sequence quality.

```bash
mkdir -p qc/trimmed_fastqc qc/trimmed_multiqc

fastqc trimmed/*.fastq.gz \
    -o qc/trimmed_fastqc \
    -t 16

multiqc qc/trimmed_fastqc \
    -o qc/trimmed_multiqc
```

The post-trimming MultiQC report can be found at:

```text
qc/trimmed_multiqc/multiqc_report.html
```

Compare this report with the initial raw-read QC report to verify the improvement in **Adapter Content** and **Per base sequence quality** before proceeding to read mapping.
<img width="1823" height="1094" alt="image" src="https://github.com/user-attachments/assets/2a6f1fdf-e34a-4d55-9659-e190d59cec25" />

