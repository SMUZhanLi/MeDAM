# MeDAM
extrapolation of potential Metabolite-Disease Associations by Mining biomedical knowledge

## Installation
``` r
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
remotes::install_github("rstudio/pool")
remotes::install_github("SMUZhanLi/MeDAM", ref = "dev", force = TRUE)

library(pool)
library(MeDAM)

## Create a pool object connected MeDAM.db
medamdb <- dbPool(drv = RSQLite::SQLite(), dbname = "MeDAM.db")

## pre-eclampsia metabolomics
data(PE.metabolomics)
abun <- PE.metabolomics$abundance
meta <- PE.metabolomics$metadata
pe_batch <- medam_batch(medamdb,
                        disease = "pre-eclampsia",
                        abundance = abun, resp = meta$Statue2,
                        compared = c("preeclampsia", "normotension"),
                        ortho = TRUE)
## IBD metabolomics
## Discovery cohort (PRISM)
data(IBD.metabolomics)
abun2 <- IBD.metabolomics$discovery$abundance
meta2 <- IBD.metabolomics$discovery$metadata

## CD v.s. Control
abun_cd_ctrl <- abun2[meta2$Diagnosis %in% c("CD", "Control"), ]
meta_cd_ctrl <- meta2[meta2$Diagnosis %in% c("CD", "Control"), ]
cd_vs_ctrl_batch <- medam_batch(medamdb,
                                disease = "Crohn's disease",
                                abundance = abun_cd_ctrl,
                                resp = meta_cd_ctrl$Diagnosis,
                                compared = c("CD", "Control"),
                                ortho = TRUE,
                                islog10 = FALSE)

## UC v.s. Control
abun_uc_ctrl <- abun2[meta2$Diagnosis %in% c("UC", "Control"), ]
meta_uc_ctrl <- meta2[meta2$Diagnosis %in% c("UC", "Control"), ]
uc_vs_ctrl_batch <- medam_batch(medamdb,
                                disease = "ulcerative colitis",
                                abundance = abun_uc_ctrl,
                                resp = meta_uc_ctrl$Diagnosis,
                                compared = c("UC", "Control"),
                                ortho = TRUE,
                                islog10 = FALSE)
```
