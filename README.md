# FreqEstimationModel

An R package developed to estimate population-level frequencies of genetic markers of antimalarial resistance.  

Here, 
- genetic markers of antimalarial resistance are alleles at single biallelic SNPs or sequences of alleles (haplotypes) over a small number of biallelic SNPs (at most seven);
- population-level parasite frequencies (vs within-infection parasite frequencies) are estimated by jointly modelling parasite genetic data on many infections, including those
that contain genetically distinct parasites (i.e., those that are polyclonal);
- parasite genetic data are prevalence data (i.e., they describe the presence of alleles - the model does not exploit information on within-infection allele frequencies). 

Beware: this readme is seven years late and hasitily added. Package documentation is bare-bones. I'm in the process of checking for bugs (I need to check expand.grid()).  

The model used to estimate frequencies was first described here:
[Taylor et al. 2014](https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-13-102)

A comprehensive discription of the model can be found here:
[Thesis_methods_chapter](https://github.com/aimeertaylor/FreqEstimationModel/blob/master/inst/Thesis_methods_chapter.pdf)

Follow the code below to install the model. Thereafter, download the files in the Inst folder and run Test_run.R from within the downloaded folder (I need to make this into a vignette). 

```r
install.packages("devtools")

devtools::install_github("aimeertaylor/FreqEstimationModel", build_vignettes = TRUE)

library(FreqEstimationModel)

help(package = "FreqEstimationModel")
```
