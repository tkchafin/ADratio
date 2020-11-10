# ADratio
Python program for calculating normalized average depth ratio ("AD ratio") per-scaffold, given depth information for a pair of samples. The utility of this method is identifying candidate scaffolds or contigs belonging to sex chromosomes or autosomes, when chromosome-level assemblies for closely related organisms are not available to map your non-model genome to. 

The idea for this approach was taken from [Bidon et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4524476/) (citation below), which is the first place I've seen it. If you are aware of an earlier reference for this method, please let me know, but for now if you use it, please cite the original authors:

Bidon T, Schreck N, Hailer F, Nilsson MA, & Janke A. 2015. Genome-wide search identifies 1.9 Mb from the Polar Bear Y chromosome for evolutionary analysis. Genome Biology and Evolution. 7(7): 2010-2022.

Below I provide a description of options for running ADratio, as well as an example pipeline for preparing the necessary inputs using BWA, Samtools, Picard, and Bedtools. If you use those programs, please cite the original authors. 

## Installation 

ADratio is a Python3 program with the following dependencies:
- numpy
- seaborn

Installation instructions coming soon.

## Allele depth ratios

Discussion coming soon

## Running ADratio 

All running options for ADratio can be viewed in the command-line menu by calling the program using the -h flag:
```
tyler:ADratio $ ./ADratio.py -h

Exiting because help menu was called.

ADratio.py

Author: Tyler K Chafin, University of Arkansas
Contact: tkchafin@uark.edu
Description: Computes allele depth ratios from pileup data

	Mandatory arguments:
		-r	: Reference FASTA file
		-1	: Sample 1 coverage file
		-2	: Sample 2 coverage file
	Optional Arguments:
		-c	: Normalizing constant, calculated as:
			  # Sample 1 reads / # Sample 2 reads [Default=1.0]
		-n	: Only count non-ambiguous (N) positions in reference
		-d	: FASTA header delimiter [default=None]
		-m	: Minimum scaffold length to report [default=None]
		-M	: Maximum proportion of Ns to retain a contig [default=0.5]
		-o	: Output file prefix [default=out]
```

The meaning of these various options, and the format of the required inputs, are discussed below.

### ADratio inputs 

ADratio requires a very simple input file which can be parsed from a [bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) file, or produced using the [bedtools](https://bedtools.readthedocs.io/en/latest/) [genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html) tool. All that is required is a tab-delimited field (for each of 2 samples) formatted like so:
```
AB1011.1  0	10	1
AB1011.1  10	14	2
AB1011.1  14	18	3
...
...
...
```

The first field is the scaffold ID (which should match a corresponding header in the multi-fasta for the genome). The second field is the physical base position (starting from the left) on the scaffold, using *0-based indexing* (as is output by bedtools genomecov using the -bg option). The third field is the final position of the interval in *half-open* format (meaning the last position of a chromosome of length N would be N-1). The final field is the number of reads which piled up on that exact coordinate (i.e. the per-base read depth). 

So, reading this format, the coverage per-base for scaffold "AB1011.1" shown in the above bedGraph is:
```
11111111122223333....
```

Below I provide an example pipeline using standard bioinformatics tools for generating these inputs, starting from raw reads formatted as fastq. 

### Optional arguments



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
Next, I mapped the trimmed reads to the available [black bear genome assembly](https://www.ncbi.nlm.nih.gov/assembly/GCA_003344425.1/), which at the time of writing was scaffold-level with N=111,495 scaffolds and a scaffold N50 of ~190k. We'll be removing scaffolds below a length threshold later. I performed the mapping using [bwa mem](http://bio-bwa.sourceforge.net/):

```
#indexing the reference genome 
bwa index GCA_003344425.1_ASM334442v1_genomic.fna

#bwa mapping 
bwa mem -t 32 GCA_003344425.1_ASM334442v1_genomic.fna male_R1.fq.gz male_R2.fq.gz > male_align.sam
bwa mem -t 32 GCA_003344425.1_ASM334442v1_genomic.fna female_R1.fq.gz female_R2.fq.gz > female_align.sam
```

In my case, I was running on an HPC cluster with 32 cores per node (-t 32), so each of these was actually submitted separately as individual jobs, each on one 32-core node. You should change the "-t" flag to match the number of cores available to you.

The next step did some file conversions using [samtools](http://www.htslib.org/), namely converting the map output (.sam) to a more compressed (binary) format (.bam), and sorting that output to improve the compression. These steps will need to be replicated on both male and female sequences (note here only shown for female):
```
samtools view -h -b -S female_align.sam > female_align.bam
samtools view -b -F 4 female_align.bam > female_mapped.bam 
samtools sort female_mapped.bam > female_mapped_sorted.bam
```

Next, I used the [picard](https://broadinstitute.github.io/picard/) MarkDuplicates tool (read about it's purpose [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)) to remove PCR artefacts which can present as optical duplicates on the sequencer:
```
java -jar /share/apps/bioinformatics/picard/picard-tools-2.17.10/picard.jar MarkDuplicates I=female_mapped_sorted.bam O=female_removeDups.bam M=marked_dup_metrics.txt REMOVE_DUPLICATES=true
java -jar /share/apps/bioinformatics/picard/picard-tools-2.17.10/picard.jar MarkDuplicates I=male_mapped_sorted.bam O=male_removeDups.bam M=marked_dup_metrics.txt REMOVE_DUPLICATES=true
```

Finally, I produced the exact output needed for ADratio.py using the [genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html) tool from [bedtools](https://bedtools.readthedocs.io/en/latest/), using the "-d" option to produce a per-base depth report (using *1-based* indexing):
```
/share/apps/bioinformatics/bedtools2/2.25.0/bin/bedtools genomecov -d -ibam female_removeDups.bam > female_coverage.txt
/share/apps/bioinformatics/bedtools2/2.25.0/bin/bedtools genomecov -d -ibam male_removeDups.bam > male_coverage.txt
```

That's it! After completing these steps you can run ADratio. Note that the runtimes can be *considerable* depending on the amount of sequence data, your genome size, and the specifics of your machine (e.g. # available cores), so running on an HPC is recommended. For reference, here are rough runtimes for the example above: Trimmomatic: ~3 hours per individual; BWA+samtools: ~12 hours (with 32 cores; CPU time ~70 hours); Picard MarkDuplicates: ~4 hours (CPU time 9 hours); bedtools genomecov: 12 hours. Benchmarking ADratio for these files is also presented below. 

## ADratio runtimes and benchmarking

Coming soon
