# ADratio
Python script for calculating normalized average depth ratio ("AD ratio") per-scaffold, given depth information for a pair of samples. The utility of this method is identifying candidate scaffolds or contigs belonging to sex chromosomes or autosomes, when chromosome-level assemblies for closely related organisms are not available to map your non-model genome to. 

The idea for this approach was taken from [Bidon et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4524476/) (citation below), which is the first place I've seen it. If you are aware of an earlier reference for this method, please let me know, but for now if you use it, please cite the original authors:

Bidon T, Schreck N, Hailer F, Nilsson MA, & Janke A. 2015. Genome-wide search identifies 1.9 Mb from the Polar Bear Y chromosome for evolutionary analysis. Genome Biology and Evolution. 7(7): 2010-2022.

Below I provide a description of options for running ADratio, as well as an example pipeline for preparing the necessary inputs using BWA, Samtools, Picard, and Bedtools. If you use those programs, please cite the original authors. 

## Installation 

ADratio is a Python3 program with the following dependencies:
- pandas 
- numpy
- seaborn

## Running ADratio 

COming soon

## ADratio inputs 

ADratio requires a very simple input file which can be parsed from a samtools mpileup file, or produced using the bedtools genomecov tool. All that is required is a tab-delimited field (for each of 2 samples) formatted like so:
```
AB1011.1  1 0
AB1011.1  2 0
AB1011.1  3 10
AB1011.1  4 10
AB1011.1  5 11
AB1011.1  6 11
...
...
...
```

The first field is the scaffold ID (which should match a corresponding header in the multi-fasta for the genome). The second field is the physical base position (starting from the left) on the scaffold, using *1-based indexing* (as is output by bedtools genomecov). The final field is the number of reads which piled up on that exact coordinate (i.e. the per-base read depth). 

Below I provide an example pipeline using standard bioinformatics tools for generating these inputs, starting from raw reads formatted as fastq. 

## An example pipeline 

There are various ways to go about producing the above input, but the general process is to: 1) Map reads for each sample to the given genome (after any number of optional quality filtering steps); 2) Count up how many reads pile up on each read position. In our example here, we will accomplish the first using [bwa mem](http://bio-bwa.sourceforge.net/), a popular short-read aligner, but you could use any number of others (such as [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [bbmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)). 

In my case, I downloaded two sets of raw sequence files from the NCBI SRA: 
```
SRR830685_1.fastq.gz SRR830685_2.fastq.gz SRR7813601_1.fastq.gz SRR7813601_2.fastq.gz
```

These represent paired-end Illumina HiSeq sequences for a known female (SRR830685) and male (SRR7813601) black bear individuals, sequenced originally by [Cahill et al 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3597504/#pgen.1003345.s013) and [Srivastava et al. 2019](https://academic.oup.com/dnaresearch/article/26/1/37/5161192), respectively.

I first used [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to remove unpaired reads, do some basic quality trimming, and also trim the female sequences to a maximum length of 125, to match the shorter length of the male reads (which was sequenced at 2x125 instead of 2x150). 

```
#trimming male sequence (note NO cropping)
java -jar /share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 32 SRR7813601_1.fastq.gz SRR7813601_2.fastq.gz male_R1.fq.gz male_unpaired_R1.fq.gz male_R2.fq.gz male_unpaired_R2.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15

#trimming female sequences (with crop to 125bp)
java -jar /share/apps/bioinformatics/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 32 SRR830685_1.fastq.gz SRR830685_2.fastq.gz female_R1.fq.gz female_unpaired_R1.fq.gz female_R2.fq.gz female_unpaired_R2.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15 CROP:125
```


