# ADratio
Python program for calculating normalized average depth ratio ("AD ratio") per-scaffold, given depth information for a pair of samples, and probabilistically classifying scaffolds based on ADratio values. The utility of this method is identifying candidate scaffolds or contigs being sex-linked, when chromosome-level assemblies for closely related organisms are not available to map your non-model genome to. 

The idea for this approach comes from [Bidon et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4524476/), which is the first place I've seen it. If you are aware of an earlier reference for this method, please let me know, but for now if you use it, please cite:

Bidon T, Schreck N, Hailer F, Nilsson MA, & Janke A. 2015. Genome-wide search identifies 1.9 Mb from the Polar Bear Y chromosome for evolutionary analysis. Genome Biology and Evolution. 7(7): 2010-2022.

Below I provide a description of options for running ADratio, as well as an example pipeline for preparing the necessary inputs using BWA, Samtools, Picard, and Bedtools. If you use those programs, please cite the original authors. 

## Installation 

ADratio is a Python3 program with the following dependencies:
- numpy
- seaborn (only for plotting functions)
- pandas

These can be installed via [anaconda](https://www.anaconda.com/products/individual), or [pip](https://pip.pypa.io/en/stable/):
```
#conda installation
conda install -c anaconda -c conda-forge numpy seaborn pandas

#or, pip3 installation
pip3 install numpy pandas seaborn
```

## Allele depth ratios
The 'allele depth ratio' can be used to identify sex-linked chromosomes in [heterogametic](https://en.wikipedia.org/wiki/Heterogametic_sex) species. For example, in a species in which males have XY and females XX, the results of whole-genome shotgun sequencing should produce roughly twice as many X chromosome reads for the female than the male (2:1 ratio), and essentially no coverage of the Y chromosome in females (0:1). Autosomes, however, are expected to be ~1:1 between the two sexes. This expectation allows us to potentially classify scaffolds or contigs as being sex-linked without having access to a chromosome-level assembly. To do so, we simply calculate for *each scaffold* the ratio female mean depth / male mean depth, multiplied by a normalizing constant (=total male reads/total female reads). The normalizing constant is required because biased sequencing effort towards one of the two samples will skew the mean depths. In an XY system, the X chromosome is expected to have an ADratio ~= 2.0; the Y chromosome ~= 0.0; and an autosome ~= 2.0. 

## Classification in ADratio

In addition to computing allele-depth ratios, ADratio also uses a companion library (nbClassifier) to build a [Naive Bayes classifier](https://en.wikipedia.org/wiki/Naive_Bayes_classifier) to use a probabilistic model to assign class labels (e.g. X, Y, auto) to input scaffolds. The classifier in ADratio assumes Gaussian priors to compute posterior probabilities of assignment to each of any arbitrary number of input classes. By default, ADratio has priors for X, Y, and autosomal with means or 1.0, 0.0, and 2.0, respectively, all with standard deviations of 0.1. However, these are all customizable by the user to fit any potential scheme (e.g. Z/W systems), and any Gaussian parameters. Options are also included for the user to easily fit the classifier to an input dataset (which could represent simulated or empirical values). More on that below. 

After computing posterior probabilities, ADratio then computes a relative evidence measurement, [Jayne's (2002) evidence](https://books.google.com/books?hl=en&lr=&id=UjsgAwAAQBAJ&oi=fnd&pg=PR17&ots=OVLtaP_kaC&sig=aVU96BIUz5jq6crltDQrBvBdlRE#v=onepage&q&f=false), of membership to each class as P(class) / P(not_class); where P(not_class) is the product of assignment probabilities to all other classes. For example, with 3 classes (A, B, C), this would be equal to P(A|x) / P(B|x) * P(C|x), and so on. For interpretability, the result is transformed as 10 times the base-10 logarithm. This classification scheme was inspired by [Tyrell 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6526653/).

# Running ADratio 

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
		
	Optional arguments:
		-d	: FASTA header delimiter [default=None]
		-o	: Output file prefix [default=out]
		-x	: Toggle to turn OFF plotting
		-b	: Binwidth for plotting [default=0.1]
		
	AD ratio arguments:
		-c	: Normalizing constant, calculated as:
			  # Sample 2 reads / # Sample 1 reads [Default=1.0]
		-n	: Only count non-ambiguous (N) positions in reference
		-m	: Minimum scaffold length to report [default=None]
		-M	: Maximum proportion of Ns to retain a contig [default=0.5]
		
	Classifier arguments:
		-N	: Toggle to classify scaffolds to chromosome type (e.g. X, Y, autosome)
		-p	: (Optional) Params file to customize chr type priors
			   See documentation. By default, we assume three Gaussian 
			   priors representing how we expect ADratio to vary by chr type:
			   Class	AD_mean	AD_sd	Prob
			   X	2.0	0.1	1.0
			   Y	0.0	0.1	1.0
			   auto	1.0	0.1	1.0
			   NOTE: Here we assume <-1> female and <-2> male.
		-F	: (Optional) Fit classifier to a tab-delimited data file
			   See documentation. Format should be like so: 
			   Class AD
			   X	2.09
			   X	1.99
			   Y	0.001
			   auto	1.1
			   auto	0.977
			   ...
		-f	: Toggle on to set class probabilities equal when using -F
		-P	: Maximum a posteriori (MAP) threshold to keep a classification 
		-J	: Toggle on to calculate Jayne's 'evidence' for each class
		-j	: Jayne's evidence (db) threshold [default=30]
```

The meaning of these various options, and the format of the required inputs, are discussed below. Also note that this assumes your Python3 interpreter is installed at /usr/bin/python... If this isn't the case, you can specify python like so:
```
#if python3 in your path:
python3 ADratio.py -h

#if python3 not in your path (for example here, at ~/miniconda3/bin
~/miniconda3/bin/python3 ADratio.py -h
```

### Input file

ADratio requires a very simple input file which can be parsed from a [bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) file, or produced using the [bedtools](https://bedtools.readthedocs.io/en/latest/) [genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html) tool. Note that bedGraph uses a *zero-relative, half-open coordinate system*. More on what that means below. The bedGraph (for each of 2 samples) should be a tab-delimited text file with four fields:
```
AB1011.1  0	10	1
AB1011.1  10	15	3
AB1011.1  15	16	4
AB1011.1	17	20	5
...
...
...
```

The first field is the scaffold ID (which should match a corresponding header in the multi-fasta for the genome). The second field is the physical base position (starting from the left) on the scaffold, using *0-based indexing* (as is output by bedtools genomecov using the -bg option). The third field is the final position of the interval in *half-open* format (meaning the last position of a chromosome of length N would be N-1). The final field is the number of reads which piled up on that exact coordinate (i.e. the per-base read depth). 

So, reading this format, the coverage per-base for scaffold "AB1011.1" shown in the above bedGraph is:
```
111111111133333440555....
```

Below I provide an example pipeline using standard bioinformatics tools for generating these inputs, starting from raw reads formatted as fastq. 

### ADratio running options

ADratio includes a few basic filtering options. Using the <-n> flag, all N (ambiguous) positions in the reference scaffold sequences will be un-counted; meaning that they are both removed from the input bedGraphs, AND subtracted from the total scaffold length. This is to prevent bias caused by mapping to ambiguous positions. The user can also remove scaffolds below a certain length threshold <-m>, and specify a maximum allowable N proportion to retain a scaffold <-N>. 

For computation of the allele-depth ratio, the user can also specify a normalizing constant using <-c>. This is left to the user because the user may have certain criteria for which reads were used for mapping (e.g. removing duplicates), and can more easily calculate read counts using other tools such as samtools. 

ADratio comes with a companion library called nbClassifier which provides a Gaussian Naive Bayes classifier. This can be used to assign scaffolds to chromosome classes (X-linked, etc) using a probabilistic model that allows a user-defined degree of variation around the expected values of X=2.0; Y=0.0; auto=1.0, including options for the user to customize however they want, and including options to weight the probability of each class (i.e. unconditional on AD value). By default, ADratio generates trhee Gaussian priors, each with a standard deviation of 0.1. This is an arbitrarily chosen value. The user can also specify their own Gaussian priors using the <-p> option to point to a parameter file, formatted like so:
```
Class	AD_mean	AD_sd	Prob
X	2.0	0.1	1.0
Y	0.0	0.05	1.0
auto	1.0	0.2	1.0
```
Note that in this file, the "classes" can be anything (e.g. Z, W, X), and the values should all be non-negative. 

If you have empirical data, ADratio can fit a classifier to those observed values using the <-F> option to point to the requisite file:
```
Name    AD
X       2.00
auto    1.07
Y       0.006
Y       0.01
Y       0.23
auto    1.43
X       2.11
X       2.65
auto    0.999
auto    0.8004
auto    0.6
auto    0.99
auto    1.07
auto    1.01
auto    1.18
...
...
...
```

Here, the class probabilities P(Class) will be taken as the proportion of observations (e.g. if 1% of scaffolds are Y; P(Y)=0.01). To make them all equal across classes, you can use the <-f> flag. 

### Output files

The first outputs will be mean depths for each sample, in two files called $out_ind1_cov.txt and $out_ind2_cov.txt. These will be tab-delimited files with two columns:
```
Scaffold        MeanDepth
LZNR01000001.1  21.240256760698635
LZNR01000002.1  19.008165059060865
LZNR01000003.1  4.991910840857328
LZNR01000004.1  19.50037574791044
LZNR01000005.1  17.31112428558967
...
...
...
```

Raw AD ratios (without classification) will be output immediately after calculation in a file called $out_AD.txt. This file will be formatted in a similar fashion to the individual coverage files, but with the 2nd column containing the normalized AD ratio values for each scaffold:
```
Scaffold	AD
LZNR01000001.1  1.995234
LZNR01000002.1  0.006007
...
...
...
```

If Naive Bayes classification <-N> is turned on, results will be outputted into a table called $out_classify.txt. If Jaynes evidence calculation (-J) is turned on, that will also be present as a column. The full file in that case would have the following columns:
```
Scaffold	AD	X	Y	auto	MAP	MAP_value	X_J	Y_J	auto_J	JAYNE	JAYNE_value
```

The additional columns are the classification probability values for each class (X, Y, auto); the chosen class according to the maximum probability estimate (MAP); the probability value for the MAP classification (MAP_value); and 5 additional columns only present when <-J>. These are: the relative evidence values for each (X_J, etc); the selected (highest evidence) classification (JAYNE); and the evidence value for the selected class (JAYNE_value). 

Finally, unless plotting is turned off (using <-x>), you will receive a histogram of the raw AD values ($out_hist.pdf): 
![alt text](https://github.com/tkchafin/ADratio/blob/main/images/example_noClass.png)

Additionally, if using classification, you will get a version of this file with values colored by their classifications ($out_MAP_hist.pdf and $out_JAYNE_hist.pdf):
![alt text](https://github.com/tkchafin/ADratio/blob/main/images/example_jayne.png)

You can control the binwidth of these plots using <-b>.

### Resuming previous runs 
The script supports a few different options for resuming previous runs. This could be needed for example if something went wrong and a run didn't finish, or if you already have coverage information from a previous source and you just want to calculate allele-depth ratios. If the latter case, you should have two 'coverage' files formatted as above, and names like $out_ind1_cov.txt and $out_ind2_cov.txt (where $out = the output prefix for the run; -o). You can then skip the coverage calculation by calling ADratio with <-R 1>, which will skip straight to AD ratio calculations using the input coverage data. 

You can also just run the classification model on already calculated AD ratios, by calling with <-R 2>. In this case, there should be a file called $out_AD.txt, formatted as above. 

### Important considerations when running ADratio
*Heterochromosome scaffolds can only be assigned if they are logically present in the reference*. By this, I mean that, for example, in order to use ADratio to suggest that a scaffold is Y-linked, the individual from which the reference was sequenced should have actually had a Y-chromosome...

*Individual 1 should almost always be the female*. For most use cases, the XX individual should be <-1> and the XY individual <-2>. This is because, under the expected 0:1 ratio for Y DNA between the XX and XY individual, if they were reversed it would create a divide-by-zero error. In those cases, ADratio will *SKIP* the scaffold. So, there are actually two important notes here: 1) To identify a Y chromosome, <-1> must be XX and <-2> XY; and 2) Any scaffolds exclusively present in the <-1> individual will be excluded.

## Examples 

### Tutorial using large American black bear dataset
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
Next, I mapped the trimmed reads to the available [black bear genome assembly](https://www.ncbi.nlm.nih.gov/assembly/GCA_003344425.1/), which at the time of writing was scaffold-level with N=111,495 scaffolds and a scaffold N50 of ~190k. It goes without saying that in order to assign Y-chromosome scaffolds, they must be present in the assembly... In this case, the individual from which raw reads were generated for our reference assembly *was male*, so we're good to go! We'll be removing scaffolds below a length threshold later. I performed the mapping using [bwa mem](http://bio-bwa.sourceforge.net/):

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

The normalizing constant was then calculated using samtools applied to each of the filtered bam files, then dividing the total number of (passing) male reads by the total female reads:
```
samtools view -c -F 260 XXX_removeDups.bam
```
In this case, the normalizing constant was 1.003, so very close to the default used by ADratio. 

That's it! After completing these steps you can run ADratio. Note that the runtimes can be *considerable* depending on the amount of sequence data, your genome size, and the specifics of your machine (e.g. # available cores), so running on an HPC is recommended. For reference, here are rough runtimes for the example above: Trimmomatic: ~3 hours per individual; BWA+samtools: ~12 hours (with 32 cores; CPU time ~70 hours); Picard MarkDuplicates: ~4 hours (CPU time 9 hours); bedtools genomecov: 12 hours. Benchmarking ADratio for these files is also presented below. 

You can then run ADratio, using the default classification priors, excluding N positions, and excluding scaffolds having >10% N content or being shorter than 1000bp, like so:
```
./ADratio.py -r .GCA_003344425.1_ASM334442v1_genomic.fna -1 female_coverage.txt -2 male_coverage.txt -m 1000 -M 0.1 -N -n -J -o "bba" -d " " -c 0.997
```

### Test cases using example files

ADratio using the above large dataset (>300 million reads per sample and XXX scaffolds), ADratio took 1 hour and 13 minutes to parse both bedGraph files (~75000 scaffolds each after removing short and high N-content scaffolds) and calculate mean depth of coverage for all scaffolds. This step is a bit faster when not skipping N bases (<-n>), but I still recommend using that option as mapping depths can be misleading in scaffolds regions with high ambiguity. Calculating allele-depth ratios and performing classification using default priors took an additional XXX

#### Full run; classification with default priors

Because datasets from organisms with large genome sizes take such considerable time to parse, I've included some test files to help verify the functionality of ADratio. To do a basic test of all steps, with default classification, you can run:
```
./ADratio.py -r example/ref.fa -1 example/sample1.bedgraph -2 example/sample2.bedgraph -m 10 -M 0.5 -n -N -J
```

This will run ADratio using a very small dataset of 5 scaffolds:
```
$ cat example/ref.fa
>contig1
AAAAAAAAAAAAAAAAANAAAANNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>contig2
NANAAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>contig3
AAAAAAAAAAAAAAAAANAAAATTTAAAAAAAAAAAAAAAAAAAAAAAAAAA
>contig4
AGAAANNNNNNNNNNAA
>contig5
AAAAA
```

The <-m 10> filter removes scaffolds below 10 bases (contig 5) and <-M 0.5> removes scaffolds with >=50% N content (contig 4). You can verify these are removed in the first lines of the ADratio output: 
```
tyler:ADratio $ python3 ./ADratio.py -r example/ref.fa -1 example/sample1.bedgraph -2ple/sample2.bedgraph -m 10 -M 0.5 -n -N -J

Reading reference genome example/ref.fa

Total contigs read: 5
Contigs skipped below min length: 1
Contigs skipped above max N proportion: 1
Kept 3 contigs.
```

The <-n> option excludes N bases from retained scaffolds in the depth calculations. The next lines of ADratio output tell you the output files from parsing the bedGraphs and calculating mean coverage and AD ratios for each scaffold:
```
Parsing bedgraph for individual 1: example/sample1.bedgraph
Outputting mean coverages to: out_ind1_cov.txt

Parsing bedgraph for individual 2: example/sample2.bedgraph
Outputting mean coverages to: out_ind2_cov.txt

Computing ADratio for each scaffold...
...Using the normalizing constant: 1.0
Outputting AD results to: out_AD.txt
```

Finally, the <-N> option engages the classification model (with default priors), and the <-J> calculates Jaynes evidence:
```
Classifying scaffolds using the following priors:
  Class  AD_mean  AD_sd  classProb
0     X      2.0    0.1        1.0
1     Y      0.0    0.1        1.0
2  auto      1.0    0.1        1.0

Outputting classified results to: out_classify.txt

Done!
```
In the full output out_classify.txt, you can see the results for all of our example scaffolds:

| Scaffold | AD                    | X                      | Y                     | auto                   | MAP_value          | MAP  | X_J                | Y_J                | auto_J             | JAYNE_value        | JAYNE |
|----------|-----------------------|------------------------|-----------------------|------------------------|--------------------|------|--------------------|--------------------|--------------------|--------------------|-------|
| contig2  | 1.0756972111553786    | 1.1200753506219021e-18 | 2.980375741341931e-25 | 2.9955958827339573     | 2.9955958827339573 | auto | 60.98492995045966  | -70.51459445052521 | 429.52964965321894 | 429.52964965321894 | auto  |
| contig1  | 2.0                   | 3.989422804014327      | 5.520948362159921e-87 | 7.694598626706474e-22  | 3.989422804014327  | X    | 1079.72710409992   | -657.4508235130871 | 645.4326221966683  | 1079.72710409992   | X     |
| contig3  | 0.0062499999999999995 | 1.923240324000599e-86  | 3.9816385668688663    | 1.4347353220919023e-21 | 3.9816385668688663 | Y    | -654.7280006870919 | 1071.5925648783336 | 642.7267639988725  | 1071.5925648783336 | Y     |

#### Classification using custom priors

You can also specify your own priors. An example prior config file is provided in example/nb_config.txt:
```
Class   AD_mean AD_sd   Prob
X       2.0     0.2     0.15
Y       0.0     0.1     0.05
auto    1.0     0.2     0.8
```

You can run the test case again, but this time with these custom priors, like so:
```
python3 ./ADratio.py -r example/ref.fa -1 example/sample1.bedgraph -2 example/sample2.bedgraph -m 10 -M 0.5 -n -N -J -p example/nb_config.txt 
```

You should now see that the priors in the command-line output have changed:
```
Classifying scaffolds using the following priors:
  Class  AD_mean  AD_sd  classProb
0     X      2.0    0.2       0.15
1     Y      0.0    0.1       0.05
2  auto      1.0    0.2       0.80
```

You should also see that, although the classification didn't change, the probabilities and relative evidence values have by looking in the new output (out_classify.txt):
| Scaffold | AD                    | X                     | Y                      | auto                  | MAP_value           | MAP  | X_J                | Y_J                 | auto_J             | JAYNE_value        | JAYNE |
|----------|-----------------------|-----------------------|------------------------|-----------------------|---------------------|------|--------------------|---------------------|--------------------|--------------------|-------|
| contig2  | 1.0756972111553786    | 6.887405018617033e-06 | 1.4901878706709654e-26 | 1.4854681583307487    | 1.4854681583307487  | auto | 204.9295125162716  | -208.36677945516135 | 311.60566700068665 | 311.60566700068665 | auto  |
| contig1  | 2.0                   | 0.2992067103010745    | 2.7604741810799605e-88 | 5.946878058937202e-06 | 0.2992067103010745  | X    | 922.6069860634777  | -818.0927667306429  | 828.5733401463901  | 922.6069860634777  | X     |
| contig3  | 0.0062499999999999995 | 7.884101409912027e-23 | 0.19908192834344332    | 6.949210837924389e-06 | 0.19908192834344332 | Y    | -162.4421512525781 | 265.6034414865151   | 176.46151447763825 | 265.6034414865151  | Y     |

#### Classification using empirical allele-depth ratios
Finally, ADratio also allows setting priors using empirical data. An example empirical dataset is provided in example/nb_testfit.txt. There are a few ways you can run this: 1) Using the empirical frequencies (in addition to empirical AD values); or 2) Treating class frequencies as equal (<-f>). Here, we'll do the latter:
```
python3 ./ADratio.py -r example/ref.fa -1 example/sample1.bedgraph -2 example/sample2.bedgraph -m 10 -M 0.5 -n -N -J -F example/nb_testfit.txt -f
```

This time, ADratio should report reading the 'nb_testfit.txt' dataset, as well as the new empirical priors:
```
Fitting classifier using provided dataset: example/nb_testfit.txt

Classifying scaffolds using the following priors:
  Class   AD_mean     AD_sd  classProb
0     X  1.975294  0.245575        1.0
1     Y  0.073963  0.111899        1.0
2  auto  0.964522  0.188521        1.0
```

And, once again, the probabilities change in response to the new priors:
| Scaffold | AD                    | X                      | Y                      | auto                  | MAP_value          | MAP  | X_J                | Y_J                 | auto_J            | JAYNE_value        | JAYNE |
|----------|-----------------------|------------------------|------------------------|-----------------------|--------------------|------|--------------------|---------------------|-------------------|--------------------|-------|
| contig2  | 1.0756972111553786    | 0.0019805482214041224  | 1.412183808284889e-17  | 1.778411584521212     | 1.778411584521212  | auto | 138.96861914488574 | -143.96926471852058 | 198.0335563014142 | 198.0335563014142  | auto  |
| contig1  | 2.0                   | 1.6163204137968048     | 1.6593872092108116e-64 | 5.948821904879058e-07 | 1.6163204137968048 | X    | 702.1414875258806  | -577.6301068711439  | 573.4595577108098 | 702.1414875258806  | X     |
| contig3  | 0.0062499999999999995 | 1.7799141492436244e-14 | 2.9687218463039926     | 5.187030146205668e-06 | 2.9687218463039926 | Y    | -89.37089224985118 | 195.07251680880455  | 79.91950208219647 | 195.07251680880455 | Y     |
