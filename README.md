# EpiSC pluripotency manuscript

This repository contains an Rmarkdown document (episc-manuscript.Rmd) and the
required code and data to reproduce all the analysis, figures and some of the
tables related to the manuscript *Dissection of pluripotency and self-reneweal
master regulators for the mouse epiblast stem cells* by Li et.al.


## Software requirements

The [R-system](https://cran.r-project.org/) software package and the following
packages are required to execute the code in the *episc-manuscript.Rmd* file:

- [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html)
- [knitr](https://cran.r-project.org/web/packages/rmarkdown/index.html)
- [viper](https://www.bioconductor.org/packages/release/bioc/html/viper.html)
- [umap](https://cran.r-project.org/web/packages/umap/index.html)
- [mixtools](https://cran.r-project.org/web/packages/mixtools/index.html)
- [DESeq](https://bioconductor.org/packages/3.10/bioc/html/DESeq.html) This package has been deprecated by Bioconductor. For convenience, the source code for version 1.38 is available under the `code` directory in this repository.
- [limma](https://www.bioconductor.org/packages/release/bioc/html/limma.html)
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html)
- [gridExtra](https://cran.r-project.org/web/packages/gridExtra/)
- [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html)
- [igraph](https://igraph.org/r/)
- [GO.db](https://bioconductor.org/packages/release/data/annotation/html/GO.db.html)
- [tibble](https://cran.r-project.org/web/packages/tibble/index.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
- [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
- [openxlsx](https://cran.r-project.org/web/packages/openxlsx/index.html)
- [fitdistrplus](https://cran.r-project.org/web/packages/fitdistrplus/index.html)



The following packages are required to run the Gene Ontology analysis from scratch:

- [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html)
- [ViSEAGO](https://bioconductor.org/packages/release/bioc/html/ViSEAGO.html)


While not strictly required, the following software components are
recommended:

- [Rstudio](https://rstudio.com/)
- [Rmarkdown](https://rmarkdown.rstudio.com/)


## Updating the configuration to your system

Before running the code and rendering the Rmarkdown document, the right paths
must be defined in the configuration code chunk called `configuration`.
The value of the following variable should be properly defined:

- `repository_dir`: Directory where the repository was locally cloned.

## How to render the Rmarkdown document

The easiest way to render this document is with the *Knitr* bottom in **Rstudio**.
Alternatively, the `episc-manuscript.Rmd` file can be rendered from an R-session:

```r
library(rmarkdown)
rmarkdown::render("episc-manuscript.Rmd")
```

**NOTE:** Prior to render the document, please, edit the `episc-manuscript.Rmd`
file and update the variable `repository_dir` to the right directory path to be
used as *repository directory*.


## Data files included in this repository

Data files are included in the `/data/` directory, including:

1. `perturbation-optimization.rda`: EpiSC perturbed with retinoic acid and 10 small
molecule compounds for 24h and 36h (n=36). Also available from GEO:
[GSE197632](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197632).

2. `episc-perturbation.rda`: EpiSC cells perturbed with 6 morphogens for 24h and
small molecule compounds for 36h used to reverse engineer the EpiSC interactome (n=276).
Also available from GEO:
[GSE197414](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197414).

3. `EpiSC_interactome.rds`: EpiSC interactome as a regulon-class object.

4. `episc-differentiation_counts.rds`: Raw-counts for the EpiSC differentiation
time-course experiment (n=144). Also available from GEO:
[GSE199114](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199114).

5. `episc-differentiation-vipermat.rds`: VIPER matrix for the EpiSC differentiation
time-course, including the VIPER-inferred relative protein activity for 1,393
transcriptional regulators (rows) and 38 differentiation conditions (columns).

6. `episc-plateseq-silencing.rda`: Raw-counts for the shRNA-mediated silencing
of candidate MRs expression (n=798). Also available from GEO:
[GSE199855](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199855).

7. `incell-primary-screen.rds`: Raw-data for the Incell-based screen for Pou5f1
protein levels in EpiSC in response to shRNA-mediated MR silencing.
This dataset is a list of 22 plates, each one represented by a list of samples
with Pou5f1 protein levels for individual cells (n=23,317,327).

8. `episc-mi.rds`: Multual information between 210 potential transcriptional regulators
and 20,310 transcripts estimated by the adaptive partitioning method implemented in
[ARACNE-AP](https://github.com/califano-lab/ARACNe-AP).
This data is formatted as a list of vectors.

9. `entrez2gene.rds`: Vector of Entrez Gene ID to Gene symbol conversions.

10. `mm10.rda`: Provides the variables `genePosition` and `geneLength` with the
genomic location and average transcript length for known genes based on EntrezID.

11. `incell.txt`: Entrez GeneIDs for the 172 candidate MR genes whose shNRA-mediated
silencing was evaluated by Pou5f1 protein levels in the InCell assay.

12. `plateseq.txt`: Entrez GeneIDs for the 154 candidate MR genes whose
shRNA-mediated silencing was evaluated by PLATE-seq expression profile.

13. `nonmr.txt`: Entrez GeneIDs for the 15 known regulators of pluripotency that
are not considered candidate MRs.

14. `candidates.tsv`: Entrez GeneIDs for 220 genes selected for shRNA assays
and annotation of whether they are candidate MRs or positive controls (known
  regulators of pluripotency).

15. `gobp_results.rds`: Matrix of GO-BP enrichment analysis results, using shadow
analysis, for each of the 120 MRs represented in the pluripotency network.

16. `gobp_mm.rda`: List with GO-BP categories associated to each mouse known gene.

17. `msigdb_results.rds`: Matrix of MSigDB enrichment analysis results, using shadow
analysis, for each of the 120 MRs represented in the pluripotency network.

18. `msig_mouse_entrez_v7.rds`: MSigDB v7 database.


-- END --
