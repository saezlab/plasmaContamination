# plasmaContamination

A package to get contamination score for blood proteomic samples or differential anaylysis results

Contamination signatures obtained from: https://www.sciencedirect.com/science/article/pii/S2405471216300722

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
