#' @title (Data) Non-targeted metabolomics for pre-eclampsia
#' @description Contained two datasets, abundance and  metadata.
#' * \code{abundance} contained 60 samples (nrow) and 8407 metabolites (ncol),
#' values have been log10-transformed.
#' * \code{metadata} contained 30 patients with pre-eclampsia (PE) and
#' 30 age-matched normotensive pregnant women (Ctrl) in their third trimester.
#' Ref: fig 1.C in Huimin Zheng, et al. (2022).
#' @name PE.metabolomics
#' @aliases PE.metabolomics
#' @docType data
#' @keywords data
#' @references Huimin Zheng, et al. (2022). In silico method to maximise the
#' biological potential of understudied metabolomic biomarkers: a study in
#' pre-eclampsia.
#' @examples
#' data(PE.metabolomics)
#' abundance <- PE.metabolomics$abundance
#' metadata <- PE.metabolomics$metadata
#' unique(metadata$Statue2)
#' ## [1] "normotension" "preeclampsia"
NA


#' @title (Data) Untargeted LC–MS metabolomics for IBD cohort
#' @description A 155-member discovery cohort and a 65-member validation cohort,
#' each containing a crosssectional sampling of UC, CD and control patients.
#' * \code{discovery} Discovery cohort (PRISM)
#' * \code{validation} Validation cohort (LifeLines DEEP and NLIBD)
#' Ref: Eric A. Franzosa, et al. (2019).
#' @name IBD.metabolomics
#' @aliases IBD.metabolomics
#' @docType data
#' @keywords data
#' @references Eric A. Franzosa, et al. (2019). Gut microbiome structure and
#' metabolic activity in inflammatory bowel disease.
#' @examples
#' data(IBD.metabolomics)
#' ## Discovery cohort (PRISM)
#' discovery <- IBD.metabolomics$discovery
#' ## Validation cohort (LifeLines DEEP and NLIBD)
#' validation <- IBD.metabolomics$validation
NA