## DmBTC-package

The DmBTC package provides functions to conduct meta-analyses on the diagnostic accuracy of biomarkers for disease detection using a Bayesian approach. It allows users to perform Bayesian modeling of diagnostic accuracy, compare various models, generate comprehensive visualizations, and conduct diagnostics on convergence. This package is intended for researchers and statisticians analyzing the diagnostic power of biomarkers in a clinical setting.

 ## Package Structure
The DmBTC package consists of multiple R scripts and documentation resources to facilitate diagnostic meta-analyses:

```
├── R/
│ ├── DmBTCmeta.R
│ ├── SCSmeta.R
│ ├── Bamdit.R
│ ├── Simulation.R
│ ├── Compare_DmBTC.R
│ ├── Compare_Simulation.R
│ ├── Traceplot.R
│ ├── Forest_plot.R
│ ├── SROC_curve.R
├── DESCRIPTION
├── NAMESPACE
├── man/
├── inst/
│ └── extdata/DmBTC_model.txt
├── Data/ 
│ ├── Study1.csv
│ ├── Study2.csv
│ ├── Study3.csv
│ ├── Study4.csv
│ ├── Study5.csv
.
.
│ ├── Study45.csv
```
 
## Functions
- DmBTCmeta: Conducts a Bayesian meta-analysis on diagnostic accuracy data, pooling results from multiple studies.
- SCSmeta: Conducts a meta-analysis based on SCS model
- Bamdit: Conducts a meta-analysis based on bamdit model
- Simulation: Conducts a simulation for diagnostic accuracy data
- Compare_DmBTC: Allows comparison of multiple diagnostic accuracy models with data from included in DmBTC package.
- Compare_Simulation: Allows comparison of multiple diagnostic accuracy models with data from simulation.
- Plot_Forest: Visualizes results as forest plots for sensitivity and specificity.
- Traceplot: Creates trace plots for convergence diagnostics in Bayesian models.
- SROCplot: Visualizes results as SROC plots for sensitivity and specificity.

## Datasets
- sepsis_data: Sample dataset of diagnostic accuracy values for evaluating biomarkers in sepsis detection.

## Updates
######
######
- 2024-11-07. **v0.0.0.** A very preliminary version.

## References
- Verde PE: Meta-analysis of diagnostic test data: a bivariate Bayesian modeling approach. Stat Med 2010, 29(30):3088-3102.
- https://cran.r-project.org/web/packages/bamdit/index.html
- Furuya-Kanamori L, Kostoulas P, Doi SAR: A new method for synthesizing test accuracy data outperformed the bivariate method. J Clin Epidemiol 2021, 132:51-58.
