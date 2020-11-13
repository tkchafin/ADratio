---
title: 'ADratio: Classifying sex-linked and autosomal scaffolds from pileup data'
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
date: 12 November 2020
bibliography: paper.bib

---

# Summary

Sex-limited chromosomes (e.g., the Y-chromosome in mammals) are notoriously difficult to sequence due to their highly repetitive natures (https://www.sciencedirect.com/science/article/pii/S0168952517300197). However, identification of sex-chromosome sequences in next-generation datasets is useful for a variety of reasons, such as the generation of sex-linked markers for sex determination (<CITE>), or the study of sex-biased evolutionar processes such as dispersal, reproductive success, and hybridization (<CITE>; <CITE>). Reference discussion in https://www.sciencedirect.com/science/article/pii/S0168952517300197#bib0280

Historically, approached relied on bacterial artificial chromosome (BAC) or fluorescent in situ hybridization (FISH) approaches. However, these methods are time-consuming and labor-intensive. Modern approaches instead attempt to identify sex-chromosome sequences in genome assemblies. 

Sex-detector (https://pubmed.ncbi.nlm.nih.gov/27492231/); synteny/ identity search; 
Using large numbers of samples/ markers to identify sex-specific sequences: https://onlinelibrary.wiley.com/doi/pdf/10.1111/1755-0998.12237
Complex evolutionary trajectories of sex chromosomes across bird taxa

Other approaches use enrichment methods to specifically target sex-limited chromosomes (<CITE>), but these methods require changes to the study design prior to sequencing. Here, we focus on adapting methods for the identification of sex-linked scaffolds post hoc. 


A more accessible method for non-model organisms relies on sequencing coverage per-scaffold under the expectation that the ratio of coverage between the heterogametic (e.g., XY) and homogametic (e.g., XX) sexes should differ between sex-chromosomes and autosomes.

'allele-depth ratio' Bidon et al https://academic.oup.com/gbe/article/7/7/2010/630631
Six novel Y chromosome genes in Anopheles mosquitoes discovered by independently sequencing males and females ('Chromosome quotient' method)
  Applied in:
    Evolutionary analysis of the female-specific avian W chromosome
    Genome-wide wearch identifies 1.9 Mb from the polar bear Y chromosome for evolutionary analyses
    Complex evolutionary trajectories of sex chromosomes across bird taxa
    Identification of avian W-linked contigs by short-read sequencing
    Comparative sex chromosome genomics in snakes: differentiation, evolutionary strata, and lack of global dosage compensation
    Radical remodeling of the Y chromosome in a recent radiation of malaria mosquitoes
    Numerous transitions of sex chromosomes in diptera
   

Efficient identification of Y chromosome sequences in the human and Drosophila genomes
https://onlinelibrary.wiley.com/doi/full/10.1111/imb.12602
https://www.genetics.org/content/166/3/1291.short

# Statement of need
ADratio is XXX; useful because YYY



# Performance and testing



# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

# Acknowledgements

I acknowledge

# References
