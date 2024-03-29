# plasmaContamination

A package to get contamination score for blood proteomic samples or differential anaylysis results

Contamination signatures for coagulation, platelet and erithrocytes obtained from: https://www.embopress.org/doi/full/10.15252/emmm.201910427

# Installation
```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("saezlab/plasmaContamination")
```

# Tutorial
```r
library(plasmaContamination)
data(test_data)

conta_scores <- plasma_contamination(test_data)
```

# Citation

Pre-analytical processing of plasma and serum samples for combined proteome and metabolome analysis
https://www.frontiersin.org/articles/10.3389/fmolb.2022.961448/full

