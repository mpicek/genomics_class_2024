## Exercise 1

### "Rezipping" the Reference genome and indexing it (rule: rezip_and_index_chrom)

The data format we were provided with was not unfortunately compatible with the tools we are supposed to use. So I had to first rezip it (unzip with gzip and zip it with bgzip). Then we indexed the genome so that we can map the reads.

**Input**: Compressed reference genome file (data/chr7.fa.gz).

**Output**: Rezipped genome and the files required by the BWA aligner (chr7.fa.gz.amb, .ann, .bwt, .pac, .sa).
**Command**: Uses gunzip and bgzip for unziping and zipping, bwa index to generate index files from the reference genome.

### Mapping Reads and Generating BAM Files (rule: bwa_map_to_BAM)

Here we align sequencing reads to the reference genome and directly output the alignments in BAM format for efficient storage.

**Input**: Reference genome and FASTQ files for each read (data/R1.fastq.gz, data/R2.fastq.gz).

**Output**: Unsorted BAM files for each read (mapped/R1.bam, mapped/R2.bam).

**Command**: Executes bwa mem to perform the alignment, then pipes (|) the SAM output to samtools view -Sb -, which converts SAM to BAM format directly, optimizing the use of disk space and processing time.

### Sorting BAM Files (rule: samtools_sort)

We want to sort BAM files by genomic coordinates. Sorting is necessary for further analyses.

**Input**: BAM files from the mapping step.

**Output**: Sorted BAM files (sorted/R1.bam, sorted/R2.bam).

**Command**: Uses samtools sort with temporary file naming based on the sample identifier. This improves the efficiency of the sorting process by managing temporary files more effectively.

### Indexing Sorted BAM Files (rule: index_bam)

Now we generate an index for each sorted BAM file, which allows for fast access to the data.

**Input**: Sorted BAM files.

**Output**: BAM index files (indexed_bam/R1.bam.bai, indexed_bam/R2.bam.bai).

**Command**: samtools index, which creates an index file for each sorted BAM for quick data retrieval.

### Annotating SNPs (rule: vep_annotate)

The `vep_annotate` rule is designed to annotate variant files (VCF format) using the Ensembl Variant Effect Predictor (VEP). I couldn't use the downloaded cache as it was too huge (I have barely any space left on my PC). VEP not only annotates the variants, it also includes prediction of the effects of variants on genes (that we will use later).

- **Input**: VCF files.

- **Output**: R1.vcf and R2.vcf - files with annotated variants.


The rule the following options (only the non-obvious one selected):

- **`--vcf`**: Output in .vcf file.
- **`--sift b`, `--polyphen b`**: These options activate SIFT and PolyPhen predictions, respectively, which assess the impact of amino acid changes on protein function, information on SNP effects.
- **`--database`**:  VEP will get the annotation data from Ensembl databases (as stated above, I couldn't use the cache).


## The Likely Causative SNP

To find the likely causative SNP, I searched for stop_gained in the outputs of vep (out/R1.vcf and out/R2.vcf).
That's because stop gained is highly significant as it leads to **a stop codon that could shorten the protein**, potentially leading to a loss of function.
They had the same one SNP on position 2915243, reference was G, the reads both had A. (quality 225).