# Bamkin: simple SNP-based sample relationship and correlation analysis

## About

Bamkin estimates concordance and similarity within a group of samples.
It is designed to be used for checking potential cross-sample contamination in sequencing studies.
This type of analysis is routinely done in genome-wide association studies (GWAS) and there are many advanced tools available when working in that context, but many experiments are not designed for genotyping.
Bamkin is a basic genotyping pipeline for a wide range of sequencing-based experiments with varying depth and distribution of coverage, such as RNA-seq or ChIP-seq.
It has minimal prerequisites and takes BAM files as input, which would probably need to be generated for any other analysis.

## Usage

Dependencies: samtools 1.2+ and R 3.3+.

Download the code from GitHub, which will create the `bamkin` directory in the current working directory:

```
git clone --depth 1 https://github.com/igordot/bamkin
```

Find all the BAM files in a specified directory and add them to a sample sheet:

```
./bamkin/gather-bams /path/to/BAMs/ > samplesheet.csv
```

The sample sheet can be modified to adjust sample names or remove samples.
The first column is the sample name and the second column is the BAM file path.

Run the analysis:

```
./bamkin/run -g genome.fa -b snps.bed -s samplesheet.csv
```

The `genome.fa` is the reference genome FASTA and should match what the reads in the BAMs were aligned to.
The `snps.bed` are the known population SNPs.
You can additionally specify the output directory (`-o`) and the number of processing threads (`-t`).

Runtime varies greatly depending on the sequencing depth and the distribution of reads, but can be a few minutes per sample for mammalian genomes with 10-50 million reads per sample and 1-5 million known SNPs.

## Output

* `*.png`: Hierarchical clustering, correlation, and PCA plots
* `snp.stats.csv`: SNP summary for every sample
* `snp.freq.csv`: SNP frequencies for all samples
* `snp.freq.filtered.csv`: variable SNPs

## Retrieving known SNPs

There is a variety resources for known SNPs.
UCSC hosts tables of reasonably common population SNPs from dbSNP (SNPs that have a minor allele frequency of at least 1% and are mapped to a single location in the reference genome assembly) for many commonly-used genomes.

To retrieve the human `hg19` SNPs:

```
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp150Common.txt.gz
```

To retrieve the mouse `mm10` SNPs:

```
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/snp142Common.txt.gz
```

The table can then be filtered (keeping only point mutations on primary chromosomes with 5% population frequency) and converted to BED format:

```
zcat snpXXXCommon.txt.gz \
| grep -F "single" \
| grep -F "maf-5" \
| grep -v "random" \
| grep -v "chrUn" \
| cut -f 2,3,4,5 \
| LC_ALL=C sort -k1,1 -k2,2n \
> snpXXXCommon.snv.maf5.bed
```

## Alternative methods

Bamkin is a basic tool designed for versatility and simplicity.
There are a number of methods available for more rigorous analysis:

* [Ancestry and Kinship Tools (AKT)](http://illumina.github.io/akt/)
* [GATK CalculateContamination](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_contamination_CalculateContamination.php)
* [Genome-wide Complex Trait Analysis (GCTA)](http://cnsgenomics.com/software/gcta/#Overview)
* [peddy](https://github.com/brentp/peddy)
* [Pedigree Reconstruction and Identification of a Maximum Unrelated Set (PRIMUS)](https://primus.gs.washington.edu/primusweb/res/documentation.html)
* [verifyBamID](https://genome.sph.umich.edu/wiki/VerifyBamID)

