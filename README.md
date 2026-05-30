# The Roadmap of the Spleen: A Meta-Analysis of Morphometric and Vascular Anatomy

**Ilias Miltiadis¹,\*, Edgar S. Kafarov², Ali S. Dadashev², Pavel Burko³, Daria A. Sukmanova⁴, Nadezhda A. Stashevskaya⁵, Oleg K. Zenin⁴**

¹ Department of Biomedicine, Neuroscience and Advanced Diagnostics (BiND), University of Palermo, Palermo, Italy.

² Department of Normal and Topographic Anatomy with Operative Surgery, Chechen State University, Grozny.

³ Department of Biomedical Sciences, Moscow University “Synergy”, Moscow.

⁴ Department of Human Anatomy, Penza State University, Penza.

⁵ RUDN University, Moscow.

*\*Correspondence:* Ilias Miltiadis (ORCID: [0000-0002-9130-3255](https://orcid.org/0000-0002-9130-3255))

## Overview
This repository contains the dataset, statistical scripts, and generated figures for the systematic review and meta-analysis of splenic morphometric parameters and vascular anatomy.

The primary objectives of the meta-analysis are to establish summary estimates for:
1. **Splenic Dimensions:** Length, Width, Thickness, and Volume (including male vs. female volume comparisons).
2. **Splenic Vasculature:** Splenic Artery Diameter and the prevalence of Primary Splenic Branches (2 vs. 3 branches).

## Repository Structure
* `/data/` - Contains the primary dataset (`spleen-database.csv`) detailing study characteristics, sample sizes, means, standard errors, and quality assessments.
* `/figures/` - Output directory for all generated Forest plots and Trim-and-Fill Funnel plots.
* `analysis.R` - The primary R script used to execute all random-effects meta-analyses, subgroup analyses (by Modality and Quality), and plot generation.

## Prerequisites & Installation
To reproduce the analyses, you will need [R](https://www.r-project.org/) (and preferably [RStudio](https://posit.co/)). The code relies on the `meta` package for statistical pooling and bias assessment.

Install the required packages by running the following in your R console:
```R
install.packages("meta")
install.packages("readr")
install.packages("grid")