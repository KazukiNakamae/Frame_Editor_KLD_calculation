# Frame_Editor_KLD_calculation

`MapKLValue.R` performed KL divergence calculation and created heatmaps.

```R
# create the enviroment
install.packages("exactRankTests")
install.packages("ggalluvial")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("stringr")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("CrispRVariants")

# run MapKLValue.R script
source("./MapKLValue.R");
```

The output is saved into the MapKLValue_output directory.







