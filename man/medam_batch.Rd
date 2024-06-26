% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/medam-batch.R
\name{medam_batch}
\alias{medam_batch}
\alias{medam_batch_manual}
\title{MeDAM batch search}
\usage{
medam_batch(
  medam,
  abundance,
  disease = NULL,
  resp,
  compared,
  ortho = FALSE,
  log10L = FALSE,
  scaling = "standard",
  vip_cutoff = 1,
  test = "t.test",
  padj_method = "none",
  padj_cutoff = 0.01,
  auc_cutoff = 0.5,
  log2fc_cutoff = 1.5,
  score = 700,
  ecpi_score = 800,
  universe = 21306,
  wgcna = NULL,
  wgcna_params = list(abundance = abundance)
)

medam_batch_manual(
  medam,
  metabolites,
  disease = NULL,
  score = 700,
  ecpi_score = 800,
  universe = 21306,
  wgcna = NULL,
  wgcna_params = list(abundance = NULL)
)
}
\arguments{
\item{medam}{A Pool object connected MeDAM.db.}

\item{abundance}{Abundance table, a matrix or data frame in which columns are
metabolites and rows ar samples.}

\item{disease}{Aliases of disease (ignored case),  e.g. "pre-eclampsia".}

\item{resp}{Response (i.e, group) to be modelled, a factor (same length as
\code{abundance} row number).}

\item{compared}{Define the condition A and condition B in "response" for
calculating fold change.
For example, set "compared" to c("case", "control"), then the fold change is
"case/control".}

\item{ortho}{FALSE (for PLS-DA case) or TRUE (for OPLS-DA).}

\item{log10L}{Whether the 'Abudnace Table' table needs to be
log10-transformed.}

\item{scaling}{Scaling methods for "abundance", one of "none", "center",
"pareto", and "standard" (default).}

\item{vip_cutoff}{Threshold of VIP score of the (O)PLS-DA model.}

\item{test}{Use t.test or wilcox.test for differential abundance analysis.}

\item{padj_method}{P-values adjusted methods, one of "holm", "hochberg",
"hommel", "bonferroni", "BH", "BY", "fdr" and "none".}

\item{padj_cutoff}{Threshold of P-value adjusted for t.test (or wilcox.test).}

\item{auc_cutoff}{Threshold of auc of ROC curve.}

\item{log2fc_cutoff}{Threshold of absolute value of log2 fold change.}

\item{score}{Threshold of significance to include a interaction, a number
between 0 and 1000 (default: 700).}

\item{ecpi_score}{Threshold of significance for metabilite-protein
experimental interactions (from 0 to 1000, default: 800). See
\code{stitch_expt_cpi()}..}

\item{universe}{Count of universal background protein-coding genes.
For example, background of homo sapiens is 21306.}

\item{wgcna}{Result of WGCNA analysis for metabolite abundance. A dataframe
contains metabolite and module columns.}

\item{wgcna_params}{Define the parameters in WGCNA analysis when the
\code{wgcna} is NULL and use \code{wgcna_analysis()} in medam batch.}

\item{metabolites}{A biomarker metabolites list.}
}
\value{
a list contained:
\itemize{
\item \code{daa} Result of the differential abundance analysis for metabolomics.
See \code{\link{diff_metabolites}} and \code{\link{wgcna_analysis}}.
\item \code{network} Interaction network. See \code{\link{string_network}} and
\code{\link{stitch_network}}.
\item \code{drgora} Disease-related genes ORA. See \code{\link{drgene_ora}}
\item \code{disease} Common name, DOID and related (entrez) genes of disease.
\item \code{enrichment} NULL by default. Use \code{\link{medam_add_enrichment}}
to excute enriment analysis and add result into \verb{$enrichment}
}
}
\description{
Batch search for potential metabolite-disease associations using
MeDAM pipeline. Ref: fig 1.C in Huimin Zheng, et al. (2022).
}
\references{
Huimin Zheng, et al. (2022). In silico method to maximise the
biological potential of understudied metabolomic biomarkers: a study in
pre-eclampsia.
}
