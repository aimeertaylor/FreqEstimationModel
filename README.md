# FreqEstimationModel

An R package developed to estimate population-level frequencies of genetic markers of antimalarial resistance where

- genetic markers of antimalarial resistance are either alleles at single biallelic SNPs or sequences of alleles (haplotypes) over a small number of biallelic SNPs (at most seven) in the malaria parasite genome;
- population-level parasite frequencies (vs within-infection parasite frequencies) are estimated by jointly modelling parasite genetic data on many infections, including those
that contain genetically distinct parasites (i.e., those that are polyclonal);
- parasite genetic data are prevalence data (i.e., they describe the presence of alleles). A non-default version of the model does exploit information on read-depths (see note on chapter 6 below), but it has not been peer-reviewed. 

## Current state

Beware: this readme is seven years late. I wrote this code during my PhD. It is the first code I ever made public. It is far from perfect, but better public than private (I hope). Please email me (aimee.taylor@pasteur.fr) if you find any bugs. Package documentation is bare-bones. In addition to this readme, I've added a couple of vignettes, which are based on inst/Test_Run.R and inst/visualise_results.R

## More information 

The model used to estimate frequencies was first described here:
[Taylor et al. 2014](https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-13-102)

The code was first made available in conjunction with this publication: 
[Taylor et al. 2016](https://academic.oup.com/ofid/article/4/1/ofw229/2282866)

A comprehensive description of the model can be found here:
[methodological chapter of PhD thesis](https://github.com/aimeertaylor/FreqEstimationModel/blob/master/inst/Thesis_methods_chapter.pdf)

Various aspects of the code are specific to unpublished chapters of my PhD ([full thesis](https://ora.ox.ac.uk/objects/uuid:c192e7cb-b6e0-4e23-a880-de46d668ef07)). For example, accounting for inter-child variability when data come from children who experience multiple episodes (chapter 5) and using read-depth information from short-read sequencing data (chapter 6). 

## Installation and use

Follow the code below to install the model. 

```r
install.packages("devtools")

devtools::install_github("aimeertaylor/FreqEstimationModel", build_vignettes = TRUE, dependencies = TRUE)

library(FreqEstimationModel)

help(package = "FreqEstimationModel")

browseVignettes("FreqEstimationModel")

# A quick example of results that can be generated using a simulated data set: 
vignette("Quick_example", "FreqEstimationModel") 

# Loads an pdf example of the long-form results that can be generated following
the quick example:
vignette("Long_results", "FreqEstimationModel") 
```

## License
[MIT](https://choosealicense.com/licenses/mit/)

## Note-to-self
- Ignoring roxygen2 Encoding: "UTF-8" warning for now as no special characters in documention (see https://stackoverflow.com/questions/51694929/warning-about-utf-8-with-roxygen2). 
