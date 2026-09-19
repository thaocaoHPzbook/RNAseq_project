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
## Mapping with STAR

The trimmed RNA-seq reads were aligned to the **GENCODE v50 GRCh38 primary assembly** using the STAR genome index generated in the previous step.

Create a directory for the STAR mapping results:

```bash
mkdir -p mapping_count
```

Map all trimmed FASTQ files:

```bash
for file in trimmed/*.fastq.gz; do

    sample=$(basename "$file" .fastq.gz)
    sample=${sample%_trimmed}

    echo "======================================"
    echo "Mapping: $sample"
    echo "Input:   $file"
    echo "======================================"

    mkdir -p "mapping_count/$sample"

    STAR \
        --genomeDir genome/star-index \
        --runThreadN 16 \
        --readFilesIn "$file" \
        --readFilesCommand zcat \
        --quantMode GeneCounts \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "mapping_count/$sample/"

done
```

Here:

- `--genomeDir genome/star-index` uses the STAR index generated from the **GENCODE v50 GRCh38 primary assembly**.
- `--readFilesCommand zcat` allows STAR to read compressed `.fastq.gz` files directly.
- `--quantMode GeneCounts` generates gene-level read counts based on the annotation incorporated during genome indexing.
- `--outSAMtype BAM SortedByCoordinate` generates coordinate-sorted BAM files.
- Each sample is stored in a separate directory under:

```text
mapping_count/
```

For each sample, STAR will generate files including:

```text
mapping_count/SRRxxxxxxx/
├── Aligned.sortedByCoord.out.bam
├── Log.final.out
├── Log.out
├── Log.progress.out
├── ReadsPerGene.out.tab
└── SJ.out.tab
```


## Brief introduction to SAM/BAM format

**SAM (Sequence Alignment/Map)** is a text-based format used to store sequencing reads aligned to a reference genome.  
**BAM (Binary Alignment/Map)** contains the same information as SAM but in a compressed binary format, requiring less storage space.

A SAM/BAM file contains two main sections:

1. **Header section** – contains metadata such as the reference genome, sorting order, and file format version.
2. **Alignment section** – contains the alignment information for individual sequencing reads.

Each alignment record contains **11 mandatory fields**:

| Col | Field | Type | Description |
|---|---|---|---|
| 1 | QNAME | String | Query template name |
| 2 | FLAG | Integer | Bitwise alignment flag |
| 3 | RNAME | String | Reference sequence name |
| 4 | POS | Integer | 1-based leftmost mapping position |
| 5 | MAPQ | Integer | Mapping quality |
| 6 | CIGAR | String | CIGAR alignment string |
| 7 | RNEXT | String | Reference name of the mate/next read |
| 8 | PNEXT | Integer | Position of the mate/next read |
| 9 | TLEN | Integer | Observed template length |
| 10 | SEQ | String | Read sequence |
| 11 | QUAL | String | Phred-scaled base quality |

Additional optional fields may appear after these 11 mandatory fields.

### View BAM files with Samtools

BAM files are binary files and therefore cannot be read directly with a normal text editor.  
Use **Samtools** to inspect and manipulate BAM files.

For example:

```bash
samtools view mapping_count/SRRxxxxxxx/Aligned.sortedByCoord.out.bam | head -5
```

> **Note:** Although SAM files are plain text and can technically be opened or edited using a text editor, manual editing is not recommended because it may break the required SAM format. Samtools should preferably be used for viewing and manipulating both SAM and BAM files.

### Common optional alignment fields

STAR may include additional fields in each alignment record, such as:

- `NH` – number of loci to which the read maps.
- `HI` – alignment index for multi-mapped reads.
- `AS` – alignment score.
- `nM` – number of mismatches between the read and the reference sequence.

Two particularly important fields for interpreting alignments are **FLAG** and **CIGAR**.

### FLAG

The `FLAG` field stores multiple alignment properties as a bitwise integer.

Common FLAG values include:

| Integer | Meaning |
|---:|---|
| 1 | Read is paired |
| 2 | Read is mapped in a proper pair |
| 4 | Read is unmapped |
| 8 | Mate is unmapped |
| 16 | Read is mapped to the reverse strand |
| 32 | Mate is mapped to the reverse strand |
| 64 | First read in a pair |
| 128 | Second read in a pair |
| 256 | Secondary alignment |
| 512 | Alignment fails quality checks |
| 1024 | PCR or optical duplicate |
| 2048 | Supplementary alignment |

Because FLAG is bitwise encoded, a single integer can represent several properties simultaneously.

### CIGAR

The **CIGAR string** describes how a read is aligned to the reference genome.

It consists of combinations of:

```text
<number><operation>
```

For example:

```text
100M
1S99M
33M1685N66M1S
```

Common CIGAR operations include:

| Operation | Description | Consumes read | Consumes reference |
|---|---|---|---|
| `M` | Alignment match or mismatch | Yes | Yes |
| `I` | Insertion relative to reference | Yes | No |
| `D` | Deletion relative to reference | No | Yes |
| `N` | Skipped region in reference | No | Yes |
| `S` | Soft clipping | Yes | No |
| `H` | Hard clipping | No | No |
| `P` | Padding | No | No |
| `=` | Sequence match | Yes | Yes |
| `X` | Sequence mismatch | Yes | Yes |

For RNA-seq, the most commonly encountered operations are:

```text
M, I, D, N, S
```

For example:

```text
1S99M
```

means that:

- the first base is **soft-clipped**;
- the remaining 99 bases are aligned to the reference.

A splice-junction alignment may look like:

```text
33M1685N66M1S
```

which means:

- 33 bases align to the reference;
- 1685 bases are skipped on the reference, typically representing an intron;
- another 66 bases align;
- the final base is soft-clipped.

The `N` operation is particularly important in RNA-seq because it allows splice-aware aligners such as STAR to represent reads spanning exon-exon junctions.

## RNA-seq quality assessment with RSeQC

After STAR alignment, additional RNA-seq–specific quality control can be performed using **RSeQC**.

RSeQC requires:

- the coordinate-sorted **BAM files** generated by STAR;
- a transcript annotation in **BED12** format;
- indexed BAM files for some analyses.

In this workflow, the BED12 annotation is generated from the same **GENCODE v50 primary assembly annotation** used for STAR indexing.

### Convert GENCODE GTF annotation to BED12

First convert the GTF annotation to `genePred` format and then to BED12 using UCSC utilities:

```bash
gtfToGenePred \
    genome/gencode.v50.primary_assembly.annotation.gtf \
    genome/gencode.v50.primary_assembly.annotation.genePred

genePredToBed \
    genome/gencode.v50.primary_assembly.annotation.genePred \
    genome/gencode.v50.primary_assembly.annotation.bed
```

The resulting annotation file used by RSeQC is:

```text
genome/gencode.v50.primary_assembly.annotation.bed
```

> **Note:** The BED annotation should be derived from the same genome build and annotation release used for STAR mapping. Here, both STAR and RSeQC use the **GENCODE v50 GRCh38 primary assembly** annotation.

> **Requirement:** `gtfToGenePred` and `genePredToBed` are UCSC utilities and must be installed or available in the environment before running these commands.

---

### Run RSeQC for all samples

Create an output directory, and Run RSeQC on all STAR-aligned BAM files:

```bash
mkdir -p rseqc
for bam in mapping_count/*/Aligned.sortedByCoord.out.bam; do

    sample=$(basename "$(dirname "$bam")")
    outdir="rseqc/$sample"

    echo "=========================================="
    echo "RSeQC: $sample"
    echo "=========================================="

    mkdir -p "$outdir"

    # 1. BAM index
    if [ ! -f "${bam}.bai" ]; then
        echo "[1/4] Creating BAM index..."

        samtools index -@ 8 "$bam"

    else
        echo "[1/4] BAM index already exists -> skip"
    fi


    # 2. Infer library strandedness
    if [ ! -s "$outdir/infer_experiment.txt" ]; then
        echo "[2/4] infer_experiment..."

        infer_experiment.py \
            -r genome/gencode.v50.primary_assembly.annotation.bed \
            -i "$bam" \
            > "$outdir/infer_experiment.txt" 2>&1

    else
        echo "[2/4] infer_experiment already exists -> skip"
    fi


    # 3. Read distribution
    if [ ! -s "$outdir/read_distribution.txt" ]; then
        echo "[3/4] read_distribution..."

        read_distribution.py \
            -i "$bam" \
            -r genome/gencode.v50.primary_assembly.annotation.bed \
            > "$outdir/read_distribution.txt" 2>&1

    else
        echo "[3/4] read_distribution already exists -> skip"
    fi


    # 4. Gene body coverage
    if [ ! -s "$outdir/geneBody.geneBodyCoverage.txt" ]; then
        echo "[4/4] geneBody_coverage..."

        geneBody_coverage.py \
            -r genome/gencode.v50.primary_assembly.annotation.bed \
            -i "$bam" \
            -o "$outdir/geneBody"

    else
        echo "[4/4] Gene body coverage already exists -> skip"
    fi


    echo "Finished: $sample"
    echo

done
```

### Main RSeQC analyses

The workflow performs three RNA-seq–specific QC analyses:

- **`infer_experiment.py`** – infers the library strandedness by examining how mapped reads overlap annotated genes.
- **`read_distribution.py`** – evaluates how reads are distributed across genomic features such as CDS, UTRs, introns, and intergenic regions.
- **`geneBody_coverage.py`** – evaluates read coverage along the gene body from the 5' to 3' end and can reveal potential coverage bias.

Before running these analyses, each BAM file is indexed using:

```bash
samtools index
```

The results for each sample will be stored separately:

```text
rseqc/
├── SRRxxxxxxx/
│   ├── infer_experiment.txt
│   ├── read_distribution.txt
│   ├── geneBody.geneBodyCoverage.txt
│   ├── geneBody.geneBodyCoverage.r
│   └── ...
├── SRRxxxxxxx/
│   └── ...
└── ...
```

### Interpretation notes

`infer_experiment.py` is especially important because the inferred strandedness should be checked before downstream gene-expression quantification. Using an incorrect strandedness setting can substantially affect read counting.

`read_distribution.py` helps determine whether the mapped reads are distributed as expected for an RNA-seq library. For example, standard mRNA-seq data are generally expected to show substantial enrichment in annotated exonic regions.

`geneBody_coverage.py` is useful for detecting strong **5' or 3' coverage bias**, which may indicate RNA degradation, library preparation bias, or sequencing-related effects.

RSeQC results should be interpreted together with the previous **FastQC/MultiQC** and **STAR mapping statistics**, rather than using any single QC metric alone.

### Summarize RSeQC results with MultiQC

After running RSeQC for all samples, use **MultiQC** to combine the QC results into a single report.

Create an output directory, and Run MultiQC on the RSeQC output directory:

```bash
mkdir -p qc/rseqc_multiqc
multiqc rseqc \
    -o qc/rseqc_multiqc
```

MultiQC will automatically detect supported RSeQC outputs from all samples, including results from:

- `infer_experiment.py`
- `read_distribution.py`
- `geneBody_coverage.py`

The combined report will be generated at:

```text
qc/rseqc_multiqc/multiqc_report.html
```

This report provides a convenient overview of RNA-seq–specific QC metrics across all samples and helps identify potential outliers or systematic biases.
<img width="1787" height="1135" alt="image" src="https://github.com/user-attachments/assets/8fb66a52-7168-43fa-93e8-795eee36ea42" />
<img width="1804" height="1027" alt="image" src="https://github.com/user-attachments/assets/b0503888-f43c-42be-9c2b-727b74179a20" />
<img width="1808" height="1109" alt="image" src="https://github.com/user-attachments/assets/32f3e541-672f-4e8f-9965-6e3f5b6a6b39" />
### Interpretation of RSeQC results

The RSeQC results indicated that the libraries were predominantly **unstranded**, with approximately equal proportions of sense and antisense reads across samples. Therefore, downstream expression quantification should be performed using an unstranded configuration.

Read distribution showed that most mapped reads were located within annotated exonic regions (CDS and UTRs), while only small proportions mapped to intronic or intergenic regions, consistent with a typical mRNA-seq dataset.

Gene body coverage profiles were generally similar across samples, although a moderate **3' coverage bias** was observed. One sample showed a stronger deviation from the overall coverage pattern and should be examined together with its FastQC and STAR mapping statistics.

Overall, no major RNA-seq–specific QC problem was identified that would prevent downstream expression analysis.


### Optional: Gene biotype composition

As an additional RNA-seq quality assessment, gene counts can also be summarized according to **GENCODE gene biotypes**, such as:

- `protein_coding`
- `lncRNA`
- `pseudogene`
- `rRNA`
- `snRNA`
- `snoRNA`

This analysis is not required for the main RNA-seq workflow, but it can provide additional information about the RNA composition of the sequencing library and help identify unexpected enrichment of particular RNA classes.

Its interpretation depends on the sample type and library preparation method. For example, whole-blood RNA-seq may additionally require examination of highly abundant globin transcripts, whereas total-RNA libraries may contain larger proportions of non-coding RNAs.

Therefore, biotype composition is used here mainly as an **optional library-characterization and QC step** rather than as a mandatory preprocessing step.

## Gene expression quantification with RSEM

Read alignment is an intermediate step in RNA-seq analysis. For downstream analyses, the main objective is usually to estimate the expression level of genes and transcripts.

Raw read counts are affected by several factors, including:

- sequencing depth;
- transcript abundance;
- transcript length;
- alternative transcript isoforms.

**RSEM (RNA-Seq by Expectation-Maximization)** estimates gene- and transcript-level expression using an expectation-maximization (EM) algorithm. This allows reads that may originate from multiple transcript isoforms to be probabilistically assigned and enables estimation of:

- expected counts;
- effective transcript/gene length;
- TPM;
- FPKM.

In this workflow, RSEM uses **STAR** internally for read alignment.

> **Note:** The RSEM reference should be generated using the same reference genome and annotation used throughout the workflow. Here, the **GENCODE v50 GRCh38 primary assembly genome and annotation** are used.

---

### Build the RSEM reference

Create a directory for the RSEM reference:

```bash
mkdir -p genome/rsem_ref
```

Prepare the RSEM reference using the GENCODE v50 annotation and GRCh38 primary assembly:

```bash
rsem-prepare-reference \
    --gtf genome/gencode.v50.primary_assembly.annotation.gtf \
    --star \
    -p 16 \
    genome/GRCh38.primary_assembly.genome.fa \
    genome/rsem_ref/GRCh38
```

Here:

- `--gtf` specifies the GENCODE gene annotation.
- `--star` generates the STAR index required for RSEM to use STAR as the aligner.
- `-p 16` uses 16 CPU threads.
- `genome/GRCh38.primary_assembly.genome.fa` is the reference genome.
- `genome/rsem_ref/GRCh38` is the prefix of the generated RSEM reference.

---

### Run RSEM for all samples

Create an output directory:

```bash
mkdir -p rsem
```

Run RSEM on all trimmed RNA-seq samples:

```bash
for file in trimmed/*.fastq.gz; do

    sample=$(basename "$file" .fastq.gz)
    sample=${sample%_trimmed}

    echo "======================================"
    echo "RSEM:  $sample"
    echo "Input: $file"
    echo "======================================"

    mkdir -p "rsem/$sample"

    rsem-calculate-expression \
        --star \
        --star-gzipped-read-file \
        --append-names \
        --no-bam-output \
        -p 16 \
        "$file" \
        genome/rsem_ref/GRCh38 \
        "rsem/$sample/sample"

done
```

The main options are:

- `--star` – uses STAR for read alignment.
- `--star-gzipped-read-file` – allows STAR to directly read compressed `.fastq.gz` files.
- `--append-names` – appends gene or transcript names to the output when annotation information is available.
- `--no-bam-output` – avoids retaining BAM alignment files generated internally by RSEM, reducing storage usage.
- `-p 16` – uses 16 CPU threads.

> **Note:** This workflow assumes **single-end RNA-seq data**. For paired-end data, `--paired-end` must be specified and both read files must be provided.

---

### RSEM output

For each sample, RSEM generates gene- and transcript-level expression estimates:

```text
rsem/
└── SRRxxxxxxx/
    ├── sample.genes.results
    ├── sample.isoforms.results
    └── ...
```

The gene-level result file can be inspected using:

```bash
head rsem/SRRxxxxxxx/sample.genes.results
```

A typical `sample.genes.results` file contains:

```text
gene_id                 transcript_id(s)        length    effective_length    expected_count    TPM     FPKM
ENSG00000000003.15      ENST00000373020.9,ENST00000494424.1,ENST00000496771.5,ENST00000612152.4,ENST00000614008.4    1803.81    1704.81    10.00    5.50    3.79
ENSG00000000005.6       ENST00000373031.5,ENST00000485971.1    873.50    774.50    0.00    0.00    0.00
ENSG00000000419.14      ENST00000371582.8,ENST00000371584.9,ENST00000371588.10,ENST00000413082.1,ENST00000466152.5,ENST00000494752.1,ENST00000681979.1,ENST00000682366.1,ENST00000682713.1,ENST00000682754.1,ENST00000683010.1,ENST00000683048.1,ENST00000683466.1,ENST00000684193.1,ENST00000684628.1,ENST00000684708.1    1056.68    957.68    26.00    25.46    17.53
ENSG00000000457.14      ENST00000367770.5,ENST00000367771.11,ENST00000367772.8,ENST00000423670.1,ENST00000470238.1    2916.00    2817.00    8.00    2.66    1.83
```

The main quantities are:

- `expected_count` – estimated number of reads/fragments assigned to the gene.
- `length` – estimated gene/transcript length.
- `effective_length` – effective length used by RSEM during abundance estimation.
- `TPM` – Transcripts Per Million.
- `FPKM` – Fragments Per Kilobase Million.

RSEM estimates transcript abundance first and then summarizes transcript isoforms to obtain gene-level expression estimates.
> **Important:** TPM and FPKM should not be used directly as input for differential expression analysis with DESeq2. For differential expression, count-based expression estimates should be used instead. RSEM results can be imported into DESeq2 using tools such as `tximport`, while TPM can be used for expression visualization and descriptive comparisons.

## Transcript-level quantification with kallisto

In addition to genome alignment-based quantification using STAR/RSEM, RNA-seq expression can also be quantified using **kallisto**.

Unlike STAR, which aligns reads to the reference genome, kallisto performs **pseudoalignment directly against the reference transcriptome**. Instead of determining the exact genomic alignment position of each read, kallisto identifies the set of transcripts that are compatible with the read sequence.

This approach is generally faster and requires less computational memory than conventional genome alignment, making it useful as an alternative RNA-seq quantification strategy.

In this workflow, the reference transcriptome is obtained from the same **GENCODE v50** release used for the genome annotation.

---

### Download the reference transcriptome and build the kallisto index

Create a directory for the transcriptome reference:

```bash
mkdir -p transcriptome
cd transcriptome
```

Download the GENCODE v50 human transcriptome:

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.transcripts.fa.gz
```

Create the kallisto index directory:

```bash
mkdir -p kallisto_index
```

Build the kallisto transcriptome index:

```bash
kallisto index \
    -i kallisto_index/grch38_gencode50 \
    gencode.v50.transcripts.fa.gz
```

Return to the project directory:

```bash
cd ..
```

---

### Run kallisto for all samples

Create the output directory:

```bash
mkdir -p kallisto
```

Run kallisto pseudoquantification for all trimmed FASTQ files:

```bash
for file in trimmed/*.fastq.gz; do

    sample=$(basename "$file" .fastq.gz)
    sample=${sample%_trimmed}

    echo "======================================"
    echo "Kallisto: $sample"
    echo "Input:    $file"
    echo "======================================"

    mkdir -p "kallisto/$sample"

    kallisto quant \
        -i transcriptome/kallisto_index/grch38_gencode50 \
        -o "kallisto/$sample" \
        --single \
        -l 400 \
        -s 40 \
        --threads=10 \
        "$file"

done
```

The main options are:

- `-i` – specifies the kallisto transcriptome index.
- `-o` – specifies the output directory.
- `--single` – indicates that the dataset contains single-end reads.
- `-l 400` – specifies the estimated mean fragment length.
- `-s 40` – specifies the estimated standard deviation of fragment length.
- `--threads=10` – uses 10 CPU threads.

> **Note:** For paired-end RNA-seq, kallisto can estimate the fragment length distribution directly from the read pairs. For single-end data, however, the mean fragment length (`-l`) and its standard deviation (`-s`) must be provided manually.
>
> In this workflow, `-l 400` and `-s 40` are used based on the expected fragment-size distribution of the library preparation protocol. If reliable library-specific fragment length information is available, those values should be used instead.

---

### kallisto output

For each sample, kallisto generates files such as:

```text
kallisto/
└── SRRxxxxxxx/
    ├── abundance.tsv
    ├── abundance.h5
    └── run_info.json
```

The main result file is:

```text
abundance.tsv
```

which contains transcript-level abundance estimates.

### Example kallisto output
| target_id | length | eff_length | est_counts | tpm |
|---|---:|---:|---:|---:|
| ENST00000456328.2\|ENSG00000223972.5\|OTTHUMG00000000961.2\|OTTHUMT00000362751.1\|DDX11L1-202\|DDX11L1\|1657\|processed_transcript\| | 1657 | 1258 | 0 | 0 |
| ENST00000450305.2\|ENSG00000223972.5\|OTTHUMG00000000961.2\|OTTHUMT00000002844.2\|DDX11L1-201\|DDX11L1\|632\|transcribed_unprocessed_pseudogene\| | 632 | 233 | 0 | 0 |
| ENST00000488147.1\|ENSG00000227232.5\|OTTHUMG00000000958.1\|OTTHUMT00000002839.1\|WASH7P-201\|WASH7P\|1351\|unprocessed_pseudogene\| | 1351 | 952 | 11.4439 | 9.70548 |
| ENST00000619216.1\|ENSG00000278267.1\|-\|-\|MIR6859-1-201\|MIR6859-1\|68\|miRNA\| | 68 | 5.21251 | 0 | 0 |
| ENST00000473358.1\|ENSG00000243485.5\|OTTHUMG00000000959.2\|OTTHUMT00000002840.1\|MIR1302-2HG-202\|MIR1302-2HG\|712\|lncRNA\| | 712 | 313 | 0 | 0 |

Unlike RSEM, kallisto directly reports expression estimates at the **transcript level** rather than providing a separate gene-level result table.

Transcript-level estimates can later be summarized to the gene level using tools such as **tximport**, together with the transcript-to-gene relationship from the GENCODE annotation.

> **Important:** TPM values are useful for describing relative transcript abundance but should not be used directly as input for DESeq2 differential expression analysis. For DESeq2, kallisto quantifications can be imported and summarized using `tximport`.

# 3. Analyze and compare RNA-seq data
## Analyze and compare RNA-seq data

After preprocessing and expression quantification, the RNA-seq data can be imported into **R** for downstream analyses, including:

- data normalization and transformation;
- sample-level quality assessment;
- PCA and clustering;
- batch-effect assessment;
- differential expression analysis;
- gene annotation;
- functional enrichment analysis;
- pathway and gene-set analysis.

### Install required R packages

Install the required CRAN packages:

```r
install.packages(c(
  "tidyverse",
  "ggrepel",
  "pbapply",
  "gplots",
  "pheatmap",
  "matrixStats",
  "patchwork",
  "msigdbr",
  "WGCNA",
  "BiocManager"
))
```

Install the required Bioconductor packages:

```r
BiocManager::install(c(
  "DESeq2",
  "edgeR",
  "limma",
  "sva",
  "biomaRt",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "tximport",
  "apeglm",
  "clusterProfiler",
  "fgsea",
  "EnhancedVolcano",
  "ComplexHeatmap"
))
```

Load the main packages:

```r
library(tidyverse)
library(DESeq2)
library(edgeR)
library(limma)
library(sva)
library(biomaRt)
library(tximport)
library(ggrepel)
library(pheatmap)
library(matrixStats)
library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)
```

### Import expression data into R

To begin the downstream analysis, import the gene expression estimates generated by **RSEM** for all 25 samples.

Here, the **TPM values** from `sample.genes.results` are combined into a single expression matrix, with genes as rows and samples as columns.

```r
samples <- list.files("rsem")

expr <- sapply(samples, function(sample){

  file <- paste0(
    "rsem/",
    sample,
    "/sample.genes.results"
  )

  quant <- read.csv(
    file,
    sep = "\t",
    header = TRUE
  )

  tpm <- setNames(
    quant$TPM,
    quant$gene_id
  )

  return(tpm)

})
```

Check the dimensions of the expression matrix:

```r
dim(expr)
```

```text
61852    25
```

The resulting matrix contains **61,852 genes across 25 RNA-seq samples**.

Inspect the first few rows:

```r
head(expr)
```

```text
                              SRR2815952 SRR2815954 SRR2815957 SRR2815958 SRR2815961 ...
ENSG00000000003.15_TSPAN6          6.76       8.26       3.44       2.90       5.07 ...
ENSG00000000005.6_TNMD             0.00       0.00       0.33       0.15       0.43 ...
ENSG00000000419.14_DPM1           25.92      39.52      36.24      35.12      37.71 ...
ENSG00000000457.14_SCYL3           1.70       2.85       4.74       6.28       5.52 ...
ENSG00000000460.17_C1orf112        0.00       0.90       0.64       0.49       0.33 ...
ENSG00000000938.13_FGR             6.80       2.60       1.91       2.05       2.47 ...
```

> **Note:** TPM values are useful for exploratory analyses, visualization, clustering, and comparison of relative expression patterns. For differential expression analysis with **DESeq2**, count-based expression estimates should be used instead of TPM.

### Import sample metadata

In addition to the expression matrix, the corresponding sample metadata are required for downstream analyses.

When retrieving the SRA accessions, the metadata table was also downloaded from the SRA Run Selector as:

```text
SraRunTable.csv
```

Additional sample information, including the cortical layer of each sample, was stored separately in:

```text
meta_additional.tsv
```

The two metadata tables can be merged using the sample identifier:

```r
meta <- read.csv(
  "meta_additional.tsv",
  sep = "\t",
  header = TRUE
) %>%

  dplyr::inner_join(
    read.csv("SraRunTable.csv", header = TRUE),
    by = c("Sample" = "Sample.Name"),
    suffix = c("", ".y")
  ) %>%

  dplyr::select(
    Run,
    Individual,
    AGE,
    Sample,
    Layer
  ) %>%

  dplyr::rename(
    Age = AGE
  )
```

Check the metadata:

```r
head(meta)
```

```text
         Run Individual  Age    Sample Layer
1 SRR2815952     DS1_H1 16.7 DS1_H1_01    L1
2 SRR2815954     DS1_H1 16.7 DS1_H1_03    L2
3 SRR2815957     DS1_H1 16.7 DS1_H1_06    L3
4 SRR2815958     DS1_H1 16.7 DS1_H1_07    L4
5 SRR2815961     DS1_H1 16.7 DS1_H1_10    L5
6 SRR2815964     DS1_H1 16.7 DS1_H1_13    L6
```

Check the dimensions:

```r
dim(meta)
```

```text
25 5
```

The final metadata table contains **25 samples** and the following variables:

### Gene annotation with biomaRt

Remove the gene symbol appended by RSEM:

```r
rownames(expr) <- sub("_.*$", "", rownames(expr))
```

Retrieve gene annotation from Ensembl:

```r
library(biomaRt)
library(dplyr)

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

meta_genes <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "ensembl_gene_id_version",
    "hgnc_symbol",
    "description",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand"
  ),
  filters = "ensembl_gene_id_version",
  values = rownames(expr),
  mart = ensembl
) %>%
  right_join(
    data.frame(
      ensembl_gene_id_version = rownames(expr)
    ),
    by = "ensembl_gene_id_version"
  ) %>%
  distinct(
    ensembl_gene_id_version,
    .keep_all = TRUE
  )
```

Match the expression matrix to the annotation table:

```r
expr <- expr[
  meta_genes$ensembl_gene_id_version,
]
```

### Check expression matrix and gene annotation

Inspect the first rows of the expression matrix:

```r
head(expr)
```

| Gene | SRR2815952 | SRR2815954 | SRR2815957 | SRR2815958 | SRR2815961 | SRR2815964 | SRR2815969 | SRR2815970 | SRR2815971 | SRR2815975 | ⋯ | SRR2815992 | SRR2815993 | SRR2815997 | SRR2816000 | SRR2816005 | SRR2816008 | SRR2816011 | SRR2816014 | SRR2816017 | SRR2816023 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ENSG00000000003.15_TSPAN6 | 6.76 | 8.26 | 3.44 | 2.90 | 5.07 | 4.88 | 3.26 | 24.11 | 12.07 | 4.83 | ⋯ | 5.18 | 4.40 | 3.91 | 3.61 | 4.36 | 4.58 | 3.84 | 2.97 | 4.51 | 4.01 |
| ENSG00000000005.6_TNMD | 0.00 | 0.00 | 0.33 | 0.15 | 0.43 | 0.85 | 0.98 | 0.32 | 0.35 | 0.41 | ⋯ | 0.40 | 0.31 | 0.31 | 0.00 | 0.00 | 0.66 | 0.18 | 0.00 | 0.00 | 0.00 |
| ENSG00000000419.14_DPM1 | 25.92 | 39.52 | 36.24 | 35.12 | 37.71 | 50.59 | 32.04 | 18.59 | 36.26 | 43.54 | ⋯ | 25.89 | 27.19 | 45.77 | 35.15 | 29.80 | 34.11 | 33.35 | 33.59 | 34.54 | 25.31 |
| ENSG00000000457.14_SCYL3 | 1.70 | 2.85 | 4.74 | 6.28 | 5.52 | 1.95 | 5.06 | 6.16 | 3.41 | 5.41 | ⋯ | 4.71 | 4.61 | 3.69 | 6.69 | 4.84 | 5.97 | 5.06 | 4.40 | 4.78 | 5.23 |
| ENSG00000000460.17_C1orf112 | 0.00 | 0.90 | 0.64 | 0.49 | 0.33 | 1.58 | 2.41 | 3.94 | 0.51 | 1.49 | ⋯ | 0.44 | 0.66 | 0.52 | 2.35 | 3.45 | 0.49 | 0.62 | 0.56 | 1.24 | 1.72 |
| ENSG00000000938.13_FGR | 6.80 | 2.60 | 1.91 | 2.05 | 2.47 | 1.26 | 2.38 | 25.37 | 8.34 | 5.90 | ⋯ | 3.72 | 2.73 | 1.38 | 1.92 | 2.27 | 2.78 | 2.13 | 2.09 | 2.96 | 2.30 |

Inspect the gene annotation table:

```r
head(meta_genes)
```


| rsem_id | ensembl_gene_id_version | hgnc_symbol | gene_type | chromosome_name | start_position | end_position | strand | ensembl_gene_id |
|---|---|---|---|---|---:|---:|---|---|
| ENSG00000000003.15_TSPAN6 | ENSG00000000003.15 | TSPAN6 | protein_coding | chrX | 100627108 | 100639991 | - | ENSG00000000003 |
| ENSG00000000005.6_TNMD | ENSG00000000005.6 | TNMD | protein_coding | chrX | 100584936 | 100599885 | + | ENSG00000000005 |
| ENSG00000000419.14_DPM1 | ENSG00000000419.14 | DPM1 | protein_coding | chr20 | 50934867 | 50959140 | - | ENSG00000000419 |
| ENSG00000000457.14_SCYL3 | ENSG00000000457.14 | SCYL3 | protein_coding | chr1 | 169849631 | 169894267 | - | ENSG00000000457 |
| ENSG00000000460.17_C1orf112 | ENSG00000000460.17 | C1orf112 | protein_coding | chr1 | 169662007 | 169854080 | + | ENSG00000000460 |

### Comparison of transcriptomic profiles across samples

After importing the expression matrix into R, we can begin exploring the overall transcriptomic profiles across samples.

The expression matrix contains **61,852 genes across 25 samples**. However, many annotated genes may have very low or no detectable expression in this dataset. Therefore, we first examine the distribution of the **mean TPM across all samples**.

```r
library(ggplot2)
library(patchwork)

avg_expr <- rowMeans(expr)

df_expr <- data.frame(
  avg_expr = avg_expr,
  log_avg_expr = log10(avg_expr + 1)
)

p1 <- ggplot(df_expr, aes(x = avg_expr)) +
  geom_histogram(
    bins = 50,
    fill = "#B8E0D2",
    color = "#5C5C5C",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = median(avg_expr),
    linetype = "dashed",
    color = "#7A7A7A"
  ) +
  labs(
    title = "Distribution of mean\ngene expression",
    x = "Mean TPM across samples",
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14,
      margin = margin(b = 10)
    ),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 15, 10, 15)
  )

p2 <- ggplot(df_expr, aes(x = log_avg_expr)) +
  geom_histogram(
    bins = 50,
    fill = "#F6C6C6",
    color = "#5C5C5C",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = median(log10(avg_expr + 1)),
    linetype = "dashed",
    color = "#7A7A7A"
  ) +
  labs(
    title = "Log-transformed\ndistribution",
    x = expression(log[10]("Mean TPM + 1")),
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14,
      margin = margin(b = 10)
    ),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 15, 10, 15)
  )

options(repr.plot.width = 12, repr.plot.height = 5)

p1 + p2
```

<img width="1440" height="600" alt="image" src="https://github.com/user-attachments/assets/98c31f6d-ef1b-408a-ab15-fccd7b872bad" />


Because gene expression values are highly skewed, the distribution can also be visualized using log-transformed axes:

```r
ggplot(data.frame(avg_expr), aes(x = avg_expr)) +

  geom_histogram(
    bins = 50,
    fill = "#B8D8E8",
    color = "white",
    linewidth = 0.3
  ) +

  scale_x_continuous(
    breaks = c(0, 1, 10, 100, 1000, 10000, 20000),
    trans = "log1p",
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    trans = "log1p",
    expand = c(0, 0)
  ) +

  labs(
    title = "Distribution of Mean Gene Expression",
    x = "Mean TPM across samples",
    y = "Number of genes"
  ) +

  theme_minimal(base_size = 13) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    axis.title = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 15, 10, 15)
  )
```
<img width="1440" height="600" alt="image" src="https://github.com/user-attachments/assets/cddcc622-cf50-4a73-9cef-94957a271764" />


These plots help identify the large proportion of genes with very low expression and provide a basis for defining expressed genes before downstream PCA, clustering, and other transcriptomic analyses.


### Number of samples in which each gene is detected

In addition to the average expression level, we can examine how many samples each gene is detected in.

```r
library(ggplot2)

num_det <- rowSums(expr > 0)

ggplot(data.frame(num_det), aes(x = num_det)) +

  geom_histogram(
    bins = 25,
    fill = "#CDB4DB",
    color = "white",
    linewidth = 0.3
  ) +

  labs(
    title = "Number of Samples in Which Each Gene Is Detected",
    x = "Number of samples with TPM > 0",
    y = "Number of genes"
  ) +

  theme_minimal(base_size = 13) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 15, 10, 15)
  )
```
<img width="1440" height="600" alt="image" src="https://github.com/user-attachments/assets/585943f4-0d77-4dae-847a-5f4283f20c1e" />

This plot shows how consistently each gene is detected across the 25 samples and can help identify genes that are only expressed in a small number of samples.

### Filter unexpressed and lowly expressed genes

A gene is retained if it satisfies **at least one** of the following conditions:

- `TPM > 0` in at least **50% of the samples**, or
- the **mean TPM across all samples is ≥ 1**.

```r
expressed <- rowMeans(expr > 0) >= 0.5 |
             rowMeans(expr) >= 1
```

Here:
- `rowMeans(expr > 0) >= 0.5` – retains genes detected (`TPM > 0`) in at least **50% of samples**.
- `rowMeans(expr) >= 1` – retains genes with a **mean TPM ≥ 1** across all samples.
- `|` – means **OR**, so meeting either condition is sufficient.

### Gene biotype distribution

After filtering lowly expressed genes, the distribution of gene biotypes can be examined using the `gene_type` annotation.

```r
meta_genes <- meta_genes[
  expressed,
  ,
  drop = FALSE
]

meta_genes %>%
  count(gene_type, sort = TRUE) %>%
  slice_head(n = 15) %>%
  ggplot(aes(
    x = reorder(gene_type, n),
    y = n
  )) +
  geom_col(
    width = 0.75,
    fill = "#6C5CE7"
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Number of genes",
    title = "Gene biotype distribution"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )
```

- `count(gene_type, sort = TRUE)` – counts genes in each biotype.
- `slice_head(n = 15)` – keeps the 15 most abundant gene biotypes.
- `coord_flip()` – displays the bar plot horizontally for easier reading.
<img width="1440" height="600" alt="image" src="https://github.com/user-attachments/assets/2d1214bf-8122-4a39-b4ff-817bc3face94" />

### Sample correlation and hierarchical clustering

Pairwise sample correlations were calculated using both **Pearson correlation** and **Spearman rank correlation**.

```r
library(ggplot2)
library(ggdendro)
library(patchwork)

corr_pearson <- cor(
  log1p(expr),
  method = "pearson"
)

corr_spearman <- cor(
  expr,
  method = "spearman"
)
```

- `Pearson` – measures linear correlation between samples.
- `Spearman` – measures rank-based correlation and is more robust to the expression-value distribution.
- `log1p(expr)` – log-transforms TPM values before calculating Pearson correlation.

Hierarchical clustering was performed using `1 - correlation` as the distance:

```r
hcl_pearson <- hclust(
  as.dist(1 - corr_pearson)
)

hcl_spearman <- hclust(
  as.dist(1 - corr_spearman)
)
```

Visualize the clustering results:

```r
plot_dendrogram <- function(hcl, title, line_color) {

  dend <- ggdendro::dendro_data(hcl)

  ggplot() +

    geom_segment(
      data = dend$segments,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      color = line_color,
      linewidth = 0.7
    ) +

    geom_text(
      data = dend$labels,
      aes(
        x = x,
        y = y,
        label = label
      ),
      angle = 60,
      hjust = 1,
      size = 3.2
    ) +

    labs(
      title = title,
      x = NULL,
      y = "1 − Correlation"
    ) +

    scale_y_continuous(
      expand = expansion(mult = c(0.15, 0.05))
    ) +

    theme_minimal(base_size = 12) +

    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      axis.title.y = element_text(face = "bold"),
      plot.margin = margin(10, 15, 30, 15)
    )
}
```

```r
p1 <- plot_dendrogram(
  hcl_pearson,
  "Pearson Correlation",
  "#8FB9A8"
)

p2 <- plot_dendrogram(
  hcl_spearman,
  "Spearman Correlation",
  "#C6A6C9"
)

options(
  repr.plot.width = 14,
  repr.plot.height = 6
)

p1 + p2
```

Samples with more similar transcriptomic profiles cluster closer together in the dendrogram. The clustering can later be compared with sample metadata such as `Individual`, `Age`, and `Layer`.
<img width="1680" height="720" alt="image" src="https://github.com/user-attachments/assets/65367bb3-1ef0-48f3-a9ac-ca3921d623f0" />


### Hierarchical clustering with sample metadata

The Pearson- and Spearman-based dendrograms show similar clustering patterns. However, using only SRR accession numbers makes biological interpretation difficult.

Therefore, the sample labels can be replaced with metadata variables such as **Individual** and **Layer**. Here, we focus on the **Spearman correlation-based clustering**.

```r
library(ggplot2)
library(ggdendro)
library(patchwork)

# Ensure metadata follows the same sample order as the expression matrix
stopifnot(
  identical(
    colnames(expr),
    meta$Run
  )
)

plot_dendro_label <- function(
  hcl,
  meta,
  label_var,
  title,
  color
) {

  dend <- ggdendro::dendro_data(hcl)

  labels_df <- dend$labels

  labels_df$display_label <- meta[[label_var]][
    match(
      labels_df$label,
      meta$Run
    )
  ]

  ggplot() +

    geom_segment(
      data = dend$segments,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      color = color,
      linewidth = 0.7
    ) +

    geom_text(
      data = labels_df,
      aes(
        x = x,
        y = y,
        label = display_label
      ),
      angle = 60,
      hjust = 1,
      size = 3.5
    ) +

    labs(
      title = title,
      x = NULL,
      y = "1 − Spearman correlation"
    ) +

    scale_y_continuous(
      expand = expansion(
        mult = c(0.18, 0.05)
      )
    ) +

    theme_minimal(base_size = 12) +

    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      axis.title.y = element_text(
        face = "bold"
      ),
      plot.margin = margin(
        10,
        15,
        35,
        15
      )
    )
}
```

Visualize clustering according to **Individual** and **Layer**:

```r
p_individual <- plot_dendro_label(
  hcl_spearman,
  meta,
  "Individual",
  "Clustering by Individual",
  "#92BFB1"
)

p_layer <- plot_dendro_label(
  hcl_spearman,
  meta,
  "Layer",
  "Clustering by Layer",
  "#C6A6C9"
)

options(
  repr.plot.width = 14,
  repr.plot.height = 6
)

p_individual + p_layer
```

These plots help determine whether the major transcriptomic similarities among samples are associated with **individual differences** or **cortical layer**.
<img width="1680" height="720" alt="image" src="https://github.com/user-attachments/assets/9f846a87-3509-4802-86fd-697145dd6747" />

### Principal component analysis (PCA)

Another way to compare transcriptomic similarities between samples is through **dimension reduction**. PCA summarizes the expression profiles of thousands of genes into a smaller number of principal components (PCs), allowing the major sources of variation among samples to be visualized.

Since lowly expressed genes were already removed in the previous step, PCA is performed using the filtered expression matrix.

```r
# PCA using expressed genes
pca <- prcomp(
  log1p(t(expr)),
  center = TRUE,
  scale. = TRUE
)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)
```

### Variance explained by principal components

```r
library(ggplot2)

eigs <- pca$sdev^2
prop_var <- eigs / sum(eigs)

pca_var_df <- data.frame(
  PC = seq_along(prop_var),
  Proportion = prop_var
)

ggplot(pca_var_df, aes(x = PC, y = Proportion)) +
  geom_line(
    color = "#A8BFD1",
    linewidth = 0.8
  ) +
  geom_point(
    color = "#C8B6D9",
    size = 2.8
  ) +
  scale_x_continuous(
    breaks = c(5, 10, 15, 20, 25)
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.30, by = 0.05)
  ) +
  labs(
    title = "PCA Variance Explained",
    x = "PC",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
```

The first principal components capture the largest proportions of transcriptomic variation, allowing the relationships among samples to be visualized in a low-dimensional space.
<img width="1680" height="720" alt="image" src="https://github.com/user-attachments/assets/58dca788-1043-4ca8-accf-981304e4940c" />


### PCA visualization with sample metadata

Calculate the percentage of variance explained by each PC:

```r
var_explained <- (
  pca$sdev^2 / sum(pca$sdev^2)
) * 100
```

Ensure that the metadata and PCA samples are in the same order:

```r
stopifnot(
  identical(
    rownames(pca$x),
    meta$Run
  )
)
```

Combine PCA coordinates with sample metadata:

```r
pca_df <- data.frame(
  pca$x,
  meta
)
```

Visualize the first two principal components:

```r
library(ggplot2)

ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Layer,
    shape = Individual
  )
) +

  geom_point(
    size = 5.5,
    alpha = 1,
    stroke = 1.2
  ) +

  scale_color_manual(
    values = c(
      "L1" = "#6FB7BE",
      "L2" = "#7FAFE6",
      "L3" = "#9C7BC0",
      "L4" = "#E79BB7",
      "L5" = "#E86F9C",
      "L6" = "#D89B2D",
      "WM" = "#76BE8B"
    )
  ) +

  scale_shape_manual(
    values = c(
      "DS1_H1" = 16,
      "DS1_H2" = 17,
      "DS1_H3" = 15,
      "DS1_H4" = 3
    )
  ) +

  labs(
    title = "PCA of Transcriptomic Profiles",
    x = paste0(
      "PC1 (",
      round(var_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(var_explained[2], 1),
      "%)"
    ),
    color = "Layer",
    shape = "Individual"
  ) +

  theme_minimal(base_size = 13) +

  theme(
    aspect.ratio = 1,
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    axis.text = element_text(
      size = 12
    ),
    legend.title = element_text(
      face = "bold",
      size = 14
    ),
    legend.text = element_text(
      size = 12
    ),
    legend.key.size = unit(
      1.2,
      "lines"
    ),
    panel.grid.minor = element_blank()
  ) +

  guides(
    color = guide_legend(
      override.aes = list(size = 6)
    ),
    shape = guide_legend(
      override.aes = list(
        size = 6,
        color = "black"
      )
    )
  )
```

In the PCA plot, samples positioned closer together have more similar overall transcriptomic profiles. Colors represent cortical `Layer`, while point shapes represent `Individual`.
<img width="1680" height="720" alt="image" src="https://github.com/user-attachments/assets/84790409-268d-4faf-abd2-568f92060a69" />

#### Optional: Highly variable gene identification

Even after removing unexpressed and lowly expressed genes, the dataset still contains many genes. For an optional global analysis, we can further focus on **highly variable genes (HVGs)**, which show greater expression variability across samples than expected.

### Estimate gene variability

```r
estimate_variability <- function(expr){

  means <- apply(expr, 1, mean)
  vars  <- apply(expr, 1, var)

  cv2 <- vars / means^2

  minMeanForFit <- unname(
    median(means[which(cv2 > 0.3)])
  )

  useForFit <- means >= minMeanForFit

  fit <- glm.fit(
    x = cbind(
      a0 = 1,
      a1tilde = 1 / means[useForFit]
    ),
    y = cv2[useForFit],
    family = Gamma(link = "identity")
  )

  a0 <- unname(fit$coefficients["a0"])
  a1 <- unname(fit$coefficients["a1tilde"])

  df <- ncol(expr) - 1

  afit <- a1 / means + a0

  varFitRatio <- vars / (afit * means^2)

  pval <- pchisq(
    varFitRatio * df,
    df = df,
    lower.tail = FALSE
  )

  res <- data.frame(
    mean = means,
    var = vars,
    cv2 = cv2,
    useForFit = useForFit,
    pval = pval,
    padj = p.adjust(
      pval,
      method = "BH"
    ),
    row.names = rownames(expr)
  )

  return(res)
}
```

Since lowly expressed genes were already filtered, variability can be estimated directly from the filtered expression matrix:

```r
var_genes <- estimate_variability(expr)

highvar_ids <- rownames(var_genes)[
  var_genes$padj < 0.01
]

meta_genes$highvar <- meta_genes$rsem_id %in% highvar_ids
```

- `padj < 0.01` – identifies genes with significantly higher variability than expected.
- `highvar` – indicates whether each gene is classified as highly variable.

### Hierarchical clustering using highly variable genes

```r
library(ggplot2)
library(ggdendro)
library(patchwork)

expr_highvar <- expr[
  meta_genes$highvar,
  ,
  drop = FALSE
]

corr_spearman_highvar <- cor(
  expr_highvar,
  method = "spearman"
)

hcl_spearman_highvar <- hclust(
  as.dist(
    1 - corr_spearman_highvar
  )
)
```

Visualize the clustering with metadata labels:

```r
plot_dendro_highvar <- function(
  hcl,
  meta,
  label_var,
  title,
  line_color
){

  dend <- ggdendro::dendro_data(hcl)

  labels_df <- dend$labels

  labels_df$display_label <- meta[[label_var]][
    match(
      labels_df$label,
      meta$Run
    )
  ]

  ggplot() +

    geom_segment(
      data = dend$segments,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      color = line_color,
      linewidth = 0.9
    ) +

    geom_text(
      data = labels_df,
      aes(
        x = x,
        y = y,
        label = display_label
      ),
      angle = 60,
      hjust = 1,
      size = 4,
      fontface = "bold"
    ) +

    labs(
      title = title,
      x = NULL,
      y = "1 − Spearman correlation"
    ) +

    scale_y_continuous(
      expand = expansion(
        mult = c(0.18, 0.05)
      )
    ) +

    theme_minimal(base_size = 13) +

    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 15
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(
        face = "bold"
      ),
      panel.grid = element_blank(),
      plot.margin = margin(
        10,
        15,
        35,
        15
      )
    )
}
```

```r
p_individual <- plot_dendro_highvar(
  hcl_spearman_highvar,
  meta,
  "Individual",
  "Highly Variable Genes: Individual",
  "#5FA8B0"
)

p_layer <- plot_dendro_highvar(
  hcl_spearman_highvar,
  meta,
  "Layer",
  "Highly Variable Genes: Layer",
  "#9C7BC0"
)

options(
  repr.plot.width = 14,
  repr.plot.height = 6
)

p_individual + p_layer
```

This optional analysis focuses the clustering on genes that contribute the strongest variation across samples, which can make major transcriptomic patterns easier to identify.
<img width="1680" height="720" alt="image" src="https://github.com/user-attachments/assets/b6e0ce6b-df05-494d-922d-cd6a4264d25d" />


### PCA using highly variable genes

Perform PCA using only the highly variable genes:

```r
pca_highvar <- prcomp(
  log1p(t(expr_highvar)),
  center = TRUE,
  scale. = TRUE
)
```

Calculate the proportion of variance explained:

```r
var_highvar <- (
  pca_highvar$sdev^2 /
  sum(pca_highvar$sdev^2)
) * 100
```

Combine the PCA coordinates with sample metadata:

```r
pca_highvar_df <- data.frame(
  pca_highvar$x,
  meta
)
```

Visualize the first two principal components:

```r
ggplot(
  pca_highvar_df,
  aes(
    x = PC1,
    y = PC2,
    color = Layer,
    shape = Individual
  )
) +

  geom_point(
    size = 5.8,
    alpha = 1,
    stroke = 1.2
  ) +

  scale_color_manual(
    values = c(
      "L1" = "#5FA8B0",
      "L2" = "#6D9EEB",
      "L3" = "#8E6BBE",
      "L4" = "#E595B5",
      "L5" = "#D95D8A",
      "L6" = "#C98C1E",
      "WM" = "#5DAA72"
    )
  ) +

  scale_shape_manual(
    values = c(
      "DS1_H1" = 16,
      "DS1_H2" = 17,
      "DS1_H3" = 15,
      "DS1_H4" = 3
    )
  ) +

  labs(
    title = "PCA of Highly Variable Genes",
    x = paste0(
      "PC1 (",
      round(var_highvar[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(var_highvar[2], 1),
      "%)"
    ),
    color = "Layer",
    shape = "Individual"
  ) +

  theme_minimal(base_size = 14) +

  theme(
    aspect.ratio = 1,
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    axis.text = element_text(size = 12),
    legend.title = element_text(
      face = "bold",
      size = 14
    ),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  ) +

  guides(
    color = guide_legend(
      override.aes = list(size = 6)
    ),
    shape = guide_legend(
      override.aes = list(
        size = 6,
        color = "black"
      )
    )
  )
```
<img width="1680" height="720" alt="image" src="https://github.com/user-attachments/assets/ab6cde3c-0e4b-4ef7-832d-35e3c445faa5" />

