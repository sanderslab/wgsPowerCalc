# wgsPowerTest
This R package runs power calculations for the discovery of variants in whole genome sequencing data

Install the wgsPowerCalc project in R:
```{r}
devtools::install_github("stephansanders/wgsPowerTest")
library(wgsPowerCalc)
```

Set the working directory for where you want the PDFs to be printed
```{r}
setwd("/path/to/directory/")
```

Run an example to estimate the power to detect an excess of de novo protein truncating variants (PTVs) in autism spectrum disorder as the relative risk is varied:
* Maximum relative risk (R) = 25
* Sample size (N) = 5000 cases 
* Case to control ratio (r) = 1
* P-value threshold (p_thres_denovo_burden) = 0.05
* Number of de novo PTVs per person (q) = 0.0996
* Proportion of de novo PTVs that mediate risk (f) = 0.05
* Prefix for PDF of power (name) = Test1
* Color of line in PDF (col) = "#4daf4a"
```{r}
plotDnBurdenByRelativeRisk(R=25, N=5000, r=1, q=0.0996, f=0.05, p_thres_denovo_burden=0.05, name="Test1", col="#4daf4a") 
```

Run an example to estimate the power to detect an excess of rare protein truncating variants (PTVs) in autism spectrum disorder as the relative risk is varied:
* Maximum relative risk (R) = 2.5
* Sample size (N) = 20000 cases 
* Case to control ratio (r) = 1
* P-value threshold (p_thres_cc_burden) = 0.05
* Number of rare PTVs per person (q) = 5.04
* Proportion of rare PTVs that mediate risk (f) = 0.05
* Prefix for PDF of power (name) = Test2
* Color of line in PDF (col) = "#4daf4a"
```{r}
plotCcBurdenByRelativeRisk(R=2.5, N=20000, r=1, q=5.04, f=0.05, p_thres_cc_burden=0.05, name="Test2", col="#4daf4a")
```

Run an example to estimate the power to detect a specific locus through an excess of de novo protein truncating variants (PTVs) in autism spectrum disorder as the relative risk is varied:
* Maximum relative risk (R) = 25
* Sample size (N) = 5000 cases 
* Case to control ratio (r) = 1
* P-value threshold (p_thres_denovo_locus) = 2.5e-06 
* Number of de novo PTVs per person (q) = 0.0996
* Proportion of de novo PTVs that mediate risk (f) = 0.05
* Prefix for PDF of power (name) = Test3
* Color of line in PDF (col) = "#4daf4a"
```{r}
plotDnLocusBySampleSize(R=25, N=5000, r=1, q=0.0996, f=0.05, p_thres_denovo_locus=2.5e-06, name="Test3", col="#4daf4a") 
```

Run an example to estimate the power to detect a specific locus through an excess of rare protein truncating variants (PTVs) in autism spectrum disorder as the relative risk is varied:
* Maximum relative risk (R) = 2.5
* Sample size (N) = 20000 cases 
* Case to control ratio (r) = 1
* Number of possible rare variants per gene (s) = 615
* Number of possible rare PTVs per gene (f_gene) = 123
* Proportion of rare PTVs that mediate risk (f) = 0.05
* Prevalence of autism (K) = 0.01
* Mean allele frequency of rare variants (AF_bar) = 0.001
* Number of replicate simulations (N_rep) = 50
* P-value threshold (p_thres_cc_locus_single) = 1.7e-11
* Prefix for PDF of power (name) = Test4
* Color of line in PDF (col) = "#4daf4a"
```{r}
plotCcLocusSingleByRelativeRisk(R=2.5, N=20000, r=1, s=615, f_gene=123, f=0.05, K=0.01, AF_bar=0.001, N_rep=10, p_thres_cc_locus_single=1.7e-11, name="Test4", col="#4daf4a") 
```

For the analyses used to make the figures in the WGSPD manuscript, please see th vignettes.