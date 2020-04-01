# Bamkin: simple sample correlation and relationship analysis

## About

Bamkin estimates concordance and similarity within a group of sequenced samples based on SNP frequencies.
It is designed to be used for checking potential cross-sample contamination in sequencing studies.
This type of analysis is routinely done in genome-wide association studies (GWAS) and there are many advanced tools available when working in that context, but many experiments are not designed for genotyping.
Bamkin is a basic genotyping pipeline for a wide range of sequencing-based experiments with varying depth and distribution of coverage, such as RNA-seq or ChIP-seq.

Bamkin advantages:

- compatible with sequencing experiments beyond WGS or WES, such as RNA-seq or ChIP-seq
- minimal software prerequisites
- simple input requirements (no sample genotypes or relationships needed)
- compares all samples in a batch (helpful to diagnose sample swaps)
- reasonable runtime

## Requirements

Software dependencies:

- samtools 1.9+
- R 3.3+

Required inputs:

- reference genome FASTA file
- known population SNPs in BED format (see below for instructions on how to generate)
- BAM files for each sample (mapped to the same reference genome as the FASTA)

## Usage

Download the code from GitHub, which will create the `bamkin` directory in the current working directory:

```
git clone --depth 1 https://github.com/igordot/bamkin
```

Find all the BAM files in a specified directory and add them to the sample sheet:

```
./bamkin/gather-bams /path/to/BAMs/ > samplesheet.csv
```

The sample sheet can be modified to adjust sample names or remove samples.
The first column is the sample name and the second column is the BAM file path.

Run the analysis:

```
./bamkin/run -g genome.fa -b snps.bed -s samplesheet.csv
```

You can additionally specify the output directory (`-o`) and the number of processing threads (`-t`).

Runtime will depend on the sequencing depth and the distribution of reads, but can be a few minutes per sample for mammalian genomes with 10-50 million reads per sample and 1-5 million known SNPs.

## Output

* `*.png`: Hierarchical clustering, correlation, and PCA plots
* `snp.stats.csv`: SNP summary for every sample
* `snp.freq.csv`: SNP frequencies for all samples
* `snp.freq.filtered.csv`: variable SNPs

## Retrieving known SNPs

There is a variety of resources for known SNPs.
UCSC Genome Browser is a convenient resources for many commonly-used genomes.
It provides tables of reasonably common population SNPs (dbSNP SNPs that have a minor allele frequency of at least 1% and are mapped to a single location in the reference genome assembly).

To retrieve the human `hg38` SNPs:

```
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp150Common.txt.gz
```

To retrieve the human `hg19` SNPs:

```
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp150Common.txt.gz
```

To retrieve the mouse `mm10` SNPs:

```
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/snp142Common.txt.gz
```

The table can then be filtered to remove less useful variants. This will keep only point mutations on primary chromosomes with 5% population frequency and convert the output to BED format:

```
gunzip -c snpXXXCommon.txt.gz \
| grep -F "single" \
| grep -F "maf-5" \
| grep -Fiv "mismatch" \
| grep -Fv "random" \
| grep -Fv "chrUn" \
| cut -f 2,3,4,5 \
| LC_ALL=C sort -k1,1 -k2,2n \
> snpXXXCommon.snv.maf5.bed
```

## Alternative methods

Bamkin is a basic tool designed for versatility and simplicity.
There are a number of methods available for more rigorous analysis:

* [Ancestry and Kinship Tools (AKT)](http://illumina.github.io/akt/)
* [Conpair](https://github.com/nygenome/Conpair)
* [GATK CalculateContamination](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_contamination_CalculateContamination.php)
* [Genome-wide Complex Trait Analysis (GCTA)](http://cnsgenomics.com/software/gcta/#Overview)
* [peddy](https://github.com/brentp/peddy)
* [Pedigree Reconstruction and Identification of a Maximum Unrelated Set (PRIMUS)](https://primus.gs.washington.edu/primusweb/res/documentation.html)
* [verifyBamID](https://genome.sph.umich.edu/wiki/VerifyBamID)

