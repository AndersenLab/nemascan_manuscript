# Evaluating the power and limitations of genome-wide association mapping in _C. elegans_

## Abstract
A central goal of evolutionary genetics in Caenorhabditis elegans is to understand the genetic basis of traits that contribute to adaptation and fitness. Genome-wide association (GWA) mappings scan the genome for individual genetic variants that are significantly correlated with phenotypic variation in a population, or quantitative trait loci (QTL). GWA mappings are a popular choice for quantitative genetic analyses because the QTL that are discovered segregate in natural populations. Despite numerous successful mapping experiments, the empirical performance of GWA mappings has not, to date, been formally evaluated for this species. We developed an open-source GWA mapping pipeline called NemaScan and used a simulation-based approach to provide benchmarks of mapping performance among wild C. elegans strains. Simulated trait heritability and complexity determined the spectrum of QTL detected by GWA mappings. Power to detect smaller-effect QTL increased with the number of strains sampled from the C. elegans Natural Diversity Resource (CeNDR). Population structure was a major driver of variation in GWA mapping performance, with populations shaped by recent selection exhibiting significantly lower false discovery rates than populations composed of more divergent strains. We also recapitulated previous GWA mappings of experimentally validated quantitative trait variants. Our simulation-based evaluation of GWA performance provides the community with critical context for pursuing quantitative genetic studies using CeNDR to elucidate the genetic basis of complex traits in C. elegans natural populations.

## Code
### Required software:
1. [R-v3.6.0](https://www.r-project.org/)
2. [R/tidyverse](https://www.tidyverse.org/)
3. [R/magrittr](https://magrittr.tidyverse.org/)
4. [R/cowplot](https://github.com/wilkelab/cowplot)
5. [R/rstatix](https://github.com/kassambara/rstatix)
6. [R/wesanderson](https://github.com/karthik/wesanderson)
7. [R/ggbeeswarm](https://github.com/eclarke/ggbeeswarm)
8. [R/ggrepel](https://ggrepel.slowkow.com/)
9. [R/RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)
10. [R/ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html#:~:text='ggtree'%20extends%20the%20'ggplot2,structures%20with%20their%20annotation%20data.)
11. [R/ggplotify](https://github.com/GuangchuangYu/ggplotify)
12. [R/data.table](https://rdatatable.gitlab.io/data.table/)

#### A. `algorithm.comp.R`
  **Inputs:**
  1. `File_S1.RData`

  **Outputs:**
  1. Figure 1
  2. Supplemental Figure 1
  3. Supplemental Table 2
  4. Supplemental Table 3

#### B. `architecture.R`
  **Inputs:**
  1. `File_S2.RData`

  **Outputs:**
  1. Figure 2
  2. Supplemental Figure 2

#### C. `sample.size.R`
  **Inputs:**
  1. `File_S3.RData`

  **Outputs:**
  1. Figure 3
  2. Table 1
  3. Supplemental Figure 3
  4. Supplemental Table 4

#### D. `pop.demo.144.R`
  **Inputs:**
  1. `File_S4.RData`

  **Outputs:**
  1. Figure 4
  2. Supplemental Table 5

#### E. `local.performance.R`
  **Inputs:**
  1. `File_S5.RData`
  2. `File_S6.RData`

  **Outputs:**
  1. Figure 5
  2. Supplemental Table 6
  3. Supplemental Table 7
  4. Supplemental Table 8
  5. Supplemental Table 9
  6. Supplemental Table 10

#### F. `recapitulating.mappings.R`
  **Inputs:**
  1. `File_S7.RData`
  2. `File_S8.RData`

  **Outputs:**
  1. Figure 6
  2. Supplemental Figure 4
  3. Supplemental Figure 5
