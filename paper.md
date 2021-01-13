---
title: 'ADratio: Probabilistic classification of sex-linked scaffolds from pileup data'
tags:
  - Python
  - genomics
  - sex chromosomes
  - Next-generation sequencing
  - Read pileups
authors:
  - name: Tyler K. Chafin^[Corresponding author]
    orcid: 0000-0001-8687-5905
    affiliation: "1" 
affiliations:
 - name: Department of Biological Science, University of Arkansas
   index: 1
date: 12 January 2021
bibliography: paper.bib

---



# Statement of need
ADratio provides a framework for classifying scaffolds in existing genome assemblies as belonging to sex-chromosomes or autosomes. Contrary to previous implementations, it does so within a probabilistic classification framework, wherein posterior probabilities of assignment are computed given either theoretical priors, or priors fitted to empirical data. This is needed because it allows the 'salvaging' of sex-chromosome sequences in genomes post hoc, without the need for sequencing large numbers of individuals, or expensive targeted-enrichment assays. Code is available open-source at github.com/tkchafin/ADratio. 


# Summary

Sex-limited chromosomes (e.g., the Y-chromosome in mammals) are notoriously difficult to sequence due to their highly repetitive natures [@Tomaszkiewicz2017]. However, identification of sex-chromosome sequences or sex-linked scaffolds in genomic datasets is useful for a varierty of reasons, such as the generation of sex-linked markers for sex determination [@Krueger-Hadfield2020], contrasting sex-chromosome and autosomal evolution [@Pennell2018], or study of functional genomic traits underlying reproductive isolation [@Liu2018]. Historically, approaches to do so relied on bacterial artificial chromosome (BAC) or fluorescent in situ hybridization (FISH) approaches [@Tomaszkiewicz2017]. However, these methods are time-consuming and labor-intensive. 

Modern approaches instead attempt to identify sex-linked sequences after sequencing, often via sanalysis of large numbers of loci from individuals of known sex [@Gamble2014], and subsequently seeking statistical associations therein [@Muyle2016]. Other methods use enrichment methods to specifically target sex-limited chromosomes [@Cruz-Dávalos2018]. However, both strategies require either large numbers of samples or changes to study design prior to sequencing. Methods to extract sex-linked scaffolds from existing assemblies, with minimal numbers of individuals for which data has been acquired, are thus necessary.

A more accessible method, particularly for non-model organisms, relies on sequencing coverage under the expectation that the ratio of coverage between the heterogametic (e.g., XY) and homogametic (e.g., XX) sexes should differ between scaffolds originating from sex-chromosomes versus autosomes. Such an approach, using either 'chromosome quotient' [@Hall2013] or 'allele-depth ratio' [@Bidon2015] formulations, has proven successful in numerous organisms [@Vicoso2015, @Chen2012]. However, available implementations either lack open-source code bases or lack a statistical classification scheme. 

I here present an option, ADratio, for calcuating normalized allele-depth ratios (AD-ratios) rapidly from read pileup data from male and female sequences. It does so in a statistically meaningful manner, using a Gaussian Naive Bayes classifier (with provided complementary Python library, nbClassifier), with options for users to input custom priors if desired. This provides a robust probabilistic framework from which prior information on expected allele-depth ratios or expected occurence frequencies can be used to compute posterior probabilities of chromosome class assignment (e.g., X, Y, autosome) for candidate scaffolds or contigs. The approach may be applied using either theoretical expectations, for example that allele-depths in an XY system should exhibit 1:1 (female:male) coverage for autosomes, 0:1 for Y chromosome, and 2:1 for X-chromosomes, or by fitting the classifier to input empirical data. Specifically, the relative likelihood is first computed as the likelihood for membership of a scaffold under each class, divided by the sum of likelihoods for all other possible assignments (i.e., p(AD|X)/(p(AD|Y)+p(AD|autosome))[@Tyrrell2019] . Posterior probabilities may then be calculated using Bayes' theorem, where prior probabilties may either be computed from user-defined Gaussian priors (provided using -p) or fitted from empirical data (-F), and with the marginal probability, being invariant across classes in this case [@Webb2005], being ignored. Users may additionally specify either posterior probability thresholds to 'keep' a classification, or a threshold with respect to the Jaynes evidence criterion [@jaynes2003probability], which is defined as 10 times the base-10 logarithm of the ratio of the maximum a posteriori (MAP) probability of class membership divided by the probability of membership to other classes [@Tyrrell2019]. This quantity may be automatically computed using the -J option.

# Performance and example
I tested ADratio using whole-genome sequencing data for two American black bear individuals (<i>Ursus americanus</i>; Male: NCBI accession SRR7813601; Female: SRR830685). Sequences were first filtered and trimmed using Trimmomatic [@Bolger2014], and subsequently mapped against the (male) reference genome (GCA_003344425) using bwa [@Li2009b]. Aligned reads were sorted and compressed using samtools [@Li2009], and PCR duplicates removed using the picard MarkDuplicates tool [@Picard2019toolkit]. These steps cumulatively took 18 hours running across 32 cores. Finally, read pileups were computed using the genomecov utility in the bedtools package [@Quinlan2010], which took an additional 12 hours.

ADratio was run directly on the bedtools output, and additionally removing scaffolds shorter than 1kb (-m 1000), or having >10% ambiguous base or gap content in the reference genome (-M 0.1). A normalizing constant (see [@Bidon2015]) was computed as the ratio of total analyzed reads for the female, divided by those of the male (and provided to ADratio.py using -c), and skipping ambiguous (N) bases (-n). With this dataset, comprising >200 million rads per sample, took 1 hour and 13 minutes to process in ADratio, which generated results for >75,000 scaffolds after filtering. Of note, ADratio may be re-run from pre-processed data (via -R), which leads to dramatic speed improvements, e.g., when re-running using alternate priors. For example, re-running this dataset after only the first step (-R 1) with alternate priors completed in only 30 seconds, including plotting and writing outputs. ADratio thus provides a flexible and statistically robust method for rapidly classifying sex-chromosomes in genome assemblies from pileup data. Tutorials to run ADratio using the above pipeline, as well as to validate its behavior from example data, are available with the complete code at: github.com/tkchafin/ADratio. 


# Acknowledgements

This research is supported by the Arkansas High Performance Computing Center (AHPCC) which is funded through multiple National Science Foundation grants and the Arkansas Economic Development Commission.

# References
