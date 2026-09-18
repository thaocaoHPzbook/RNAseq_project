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

conda activate rnaseq
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


## Read duplication and deduplication

High sequence duplication can be observed in RNA-seq data because transcripts are expressed at very different abundance levels. Highly expressed genes may naturally generate many identical or near-identical reads.

Therefore, a high **Sequence Duplication Level** in FastQC does not necessarily indicate a technical problem.

For standard bulk RNA-seq data, **deduplication is generally not recommended unless there is a clear technical reason to do so**, such as:

- strong evidence of PCR over-amplification,
- a library preparation protocol specifically designed for duplicate removal,
- or the presence of UMI (Unique Molecular Identifier) information.

Without UMI or other supporting evidence, identical reads cannot be reliably distinguished between true biological signal and PCR duplicates. Removing them may bias gene expression estimates, especially for highly expressed transcripts.

For this workflow, no deduplication step is performed unless a specific technical issue is identified.

## 2.2. Read mapping/pseudomapping and quantification
2-2-1 Read mapping with STAR and data quantification
(Back to top)
Once the quality of the data is confirmed, we need to convert those millions of reads per sample into the gene- or transcript-level quantification. This would need the assignment of reads to genes or transcripts. To do this, the mostly common first step is for each read, to look for the genomic region that match with the read, given the complete genomic sequences. The identified region is then most likely the region being transcribed and generate the sequenced read in the end. This step of looking for the matched genomic regions for reads is called read genome mapping or alignment.

There are different tools, or aligners, that have been developed for this purpose. The most famous examples include Tophat/Tophat2/HISAT2 and STAR. As the commonly used modern aligners, HISAT2 and STAR shares quite some features, such as their high-efficiency, and their support of soft-trimming for low-quality bases at the ends of reads. They also have their own adventage and disadventage. HISAT2 uses fewer computational resource than STAR (particularly memory) and has better support for SNPs (single-nucleotide polymorphism) that in the same locus on the genome different individuals can have different nucleotides. On the other hand, STAR is suggested to provide more accurate alignment results. It also supports varied ways for the next step to quantify transcript abundance. In this tutorial, we will use STAR to map the FASTQ files we retreived from SRA to the human genome.

Brief introduction to the STAR aligner
Before STAR was developed and got widely acknowledged, there had been other aligners being developed, with the most commonly used example being Tophat by Cole Trapnell when he was in his PhD in University of Maryland (he is now an Associate Professor in University of Washington). Those tools were great and used by many studies using RNA-seq whichhad shown its great potential but was not yet fully mature as it is today. The main problem of those tools before STAR was their speed. They might be good enough when there were several or dozens of samples to process, but not for the research project with huge consotia effort such as ENCODE (Encyclopedia of DNA Elements), which generated RNA-seq data for hundreds or even thousands of samples.

To solve the speed issue was one of the major motivations that STAR (Spliced Transcripts Alignment to Reference) was developed in 2009 by Alexander Dobin in Cold Spring Harbor Laboratory (he is now an Assistant Professor in CSHL). The STAR paper was published in 2013 in Bioinformatics. Until now, the tool is still under active improvement and maintenance.

In brief, STAR uses a two-step procedure to achieve the high-speed alignment. The first step is seed search. For each read, STAR firstly searches for its longest sub-sequence that perfectly matches at least one locus on the reference genome. This fragment is called maximal mappable prefixes (MMP). After getting the MMP for the read, STAR searches MMP again but only for the unmapped portion (i.e. the parts outside of the MMP) of the read. These two MMPs obtained by the sequential search are also altogether called seeds, and that's why this step is named "seed search". This sequential search not only make it straightforward to deal with splice junctions where a read contains sequences from two separated exons in the genome, but also greatly speed up the alignment as searching for the entire read sequence is no longer needed.

<img width="795" height="601" alt="image" src="https://github.com/user-attachments/assets/fed751e7-46c3-4f69-adc0-2bc351ed5d05" />

Figure 1 in the STAR paper

To further deal with the possible mismatches (due to sequencing errors on reads, SNPs, point mutations, errors in the reference genome, etc.), when the MMP search doesn't reach the end of the read, the MMPs will serve as anchors in the genome that can be extended to allow for alignments with mismatches. If the extension procedure does not yield a good genomic alignment (due to poor quality at the ends of reads, poly-A tails, adapter sequence ligation, etc.), a soft trimming is applied (the remaining read sequence is ignored although not physically cut off).

After the seed search is done, STAR applies the second step, which is clustering, stitching and scoring. Seeds are firstly clustered based on proximity to a set of ‘anchor’ seeds. Then, seeds that map close enough around the anchors (so that it can be still considered to be an intron) are stitched together. In this way, different seeds of a read which are from different exons can be stitched. The stitching is guided by a local alignment scoring scheme to penaltze mismatches, insertions, deletions and splice junction gaps. The stitched combination with the highest score is chosen as the best alignment of a read.

More details information are available in the STAR paper (technical details in its Supplementary Materials).


### Download the reference genome

Before read alignment, download the reference genome sequence and prepare it for STAR indexing.

Create a directory for the genome files:

```bash
mkdir -p genome
cd genome
```

Download the human reference genome (hg38) from the UCSC Genome Browser:

Download the corresponding GENCODE v50 primary assembly annotation:

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gtf.gz
```

Decompress both files:

```bash
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v50.primary_assembly.annotation.gtf.gz
```

The resulting reference files are:

```text
GRCh38.primary_assembly.genome.fa
gencode.v50.primary_assembly.annotation.gtf
```
### Reference genome sources

The UCSC Genome Browser is not the only source for reference genome sequences. Other commonly used databases include:

- [Ensembl](https://www.ensembl.org/) – provides genome sequences, gene annotations, and comparative genomics data for many species.
- [GENCODE](https://www.gencodegenes.org/) – provides highly curated genome annotations and corresponding reference sequences for human and mouse.
- [FlyBase](https://flybase.org/) – provides genomic and genetic resources for *Drosophila* species.
- [WormBase](https://wormbase.org/) – provides genome sequences, annotations, and genetic information for *Caenorhabditis elegans* and related nematodes.

> **Note:** The reference genome FASTA and gene annotation GTF should be compatible with each other. They should use the same genome assembly, such as GRCh38 for human, and preferably come from the same database and release.  
>
> For example, in this workflow, both the reference genome and annotation are obtained from **GENCODE Release 50**:
>
> ```text
> GRCh38.primary_assembly.genome.fa
> gencode.v50.primary_assembly.annotation.gtf
> ```
>
> Using the genome and annotation from the same release helps avoid inconsistencies in chromosome names, genomic coordinates, scaffolds, and transcript annotations during read mapping and expression quantification.
>
> For standard human RNA-seq analysis, the **primary assembly** is commonly used because it avoids additional alternate loci and haplotype sequences that may increase ambiguous or multi-mapping reads.

### Build the STAR genome index

Create a directory for the STAR genome index, and Generate the STAR index:

```bash
mkdir -p star-index
STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir star-index \
    --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile gencode.v50.primary_assembly.annotation.gtf \
    --sjdbOverhang 99
```
**Note:** 
`--sjdbOverhang` is generally set to:

```text
read length - 1
```

For 100-bp reads:

```text
sjdbOverhang = 99
```

