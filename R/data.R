#' A sample RNA-seq dataset containing 4,000 genes and 8 individuals (4M, 4F),
#'  each with 15 immune-specific cell types present.
#'
#' @source The La Jolla Institute of Immunology - Database of Immune Cell eQTLs,
#'  Expression, and Epigenomics. \url{https://dice-database.org/}
#' @format A data frame with columns:
#' \describe{
#'  \item{Male1_MONOCYTES}{Expression levels from Classical Monocytes.}
#'  \item{Male1_B_CELL_NAIVE}{Expression levels from naive B cells.}
#'  \item{Male1_CD8_NAIVE}{Expression levels from naive CD8 T cells}
#'  \item{Male1_CD8_STIM}{Expression levels from stimulated CD8 T cells.}
#'  \item{Male1_TH2}{Expression levels from Th2 T cells.}
#'  \item{Male1_THSTAR}{Expression levels from Th1/17 T cells.}
#'  \item{Male1_TH17}{Expression levels from Th17 T cells.}
#'  \item{Male1_TH1}{Expression levels from Th1 T cells.}
#'  \item{Male1_TFH}{Expression levels from T follicular helper T cells.}
#'  \item{Male1_TREG_NAIVE}{Expression levels from naive regulatory T cells.}
#'  \item{Male1_TREG_MEM}{Expression levels from memory regulatory T cells.}
#'  \item{Male1_CD4_NAIVE}{Expression levels from naive CD4 T cells.}
#'  \item{Male1_CD4_STIM}{Expression levels from stimulated CD4 T cells.}
#'  \item{Male1_M2}{Expression levels from Non-classical Monocytes.}
#'  \item{Male1_NK}{Expression levels from Natural Killer cells.}
#'  \item{Male2_MONOCYTES}{Expression levels from Classical Monocytes.}
#'  \item{Male2_B_CELL_NAIVE}{Expression levels from naive B cells.}
#'  \item{Male2_CD8_NAIVE}{Expression levels from naive CD8 T cells}
#'  \item{Male2_CD8_STIM}{Expression levels from stimulated CD8 T cells.}
#'  \item{Male2_TH2}{Expression levels from Th2 T cells.}
#'  \item{Male2_THSTAR}{Expression levels from Th1/17 T cells.}
#'  \item{Male2_TH17}{Expression levels from Th17 T cells.}
#'  \item{Male2_TH1}{Expression levels from Th1 T cells.}
#'  \item{Male2_TFH}{Expression levels from T follicular helper T cells.}
#'  \item{Male2_TREG_NAIVE}{Expression levels from naive regulatory T cells.}
#'  \item{Male2_TREG_MEM}{Expression levels from memory regulatory T cells.}
#'  \item{Male2_CD4_NAIVE}{Expression levels from naive CD4 T cells.}
#'  \item{Male2_CD4_STIM}{Expression levels from stimulated CD4 T cells.}
#'  \item{Male2_M2}{Expression levels from Non-classical Monocytes.}
#'  \item{Male2_NK}{Expression levels from Natural Killer cells.}
#'  \item{Male3_MONOCYTES}{Expression levels from Classical Monocytes.}
#'  \item{Male3_B_CELL_NAIVE}{Expression levels from naive B cells.}
#'  \item{Male3_CD8_NAIVE}{Expression levels from naive CD8 T cells}
#'  \item{Male3_CD8_STIM}{Expression levels from stimulated CD8 T cells.}
#'  \item{Male3_TH2}{Expression levels from Th2 T cells.}
#'  \item{Male3_THSTAR}{Expression levels from Th1/17 T cells.}
#'  \item{Male3_TH17}{Expression levels from Th17 T cells.}
#'  \item{Male3_TH1}{Expression levels from Th1 T cells.}
#'  \item{Male3_TFH}{Expression levels from T follicular helper T cells.}
#'  \item{Male3_TREG_NAIVE}{Expression levels from naive regulatory T cells.}
#'  \item{Male3_TREG_MEM}{Expression levels from memory regulatory T cells.}
#'  \item{Male3_CD4_NAIVE}{Expression levels from naive CD4 T cells.}
#'  \item{Male3_CD4_STIM}{Expression levels from stimulated CD4 T cells.}
#'  \item{Male3_M2}{Expression levels from Non-classical Monocytes.}
#'  \item{Male3_NK}{Expression levels from Natural Killer cells.}
#'  \item{Male4_MONOCYTES}{Expression levels from Classical Monocytes.}
#'  \item{Male4_B_CELL_NAIVE}{Expression levels from naive B cells.}
#'  \item{Male4_CD8_NAIVE}{Expression levels from naive CD8 T cells}
#'  \item{Male4_CD8_STIM}{Expression levels from stimulated CD8 T cells.}
#'  \item{Male4_TH2}{Expression levels from Th2 T cells.}
#'  \item{Male4_THSTAR}{Expression levels from Th1/17 T cells.}
#'  \item{Male4_TH17}{Expression levels from Th17 T cells.}
#'  \item{Male4_TH1}{Expression levels from Th1 T cells.}
#'  \item{Male4_TFH}{Expression levels from T follicular helper T cells.}
#'  \item{Male4_TREG_NAIVE}{Expression levels from naive regulatory T cells.}
#'  \item{Male4_TREG_MEM}{Expression levels from memory regulatory T cells.}
#'  \item{Male4_CD4_NAIVE}{Expression levels from naive CD4 T cells.}
#'  \item{Male4_CD4_STIM}{Expression levels from stimulated CD4 T cells.}
#'  \item{Male4_M2}{Expression levels from Non-classical Monocytes.}
#'  \item{Male4_NK}{Expression levels from Natural Killer cells.}
#'  \item{Female1_MONOCYTES}{Expression levels from Classical Monocytes.}
#'  \item{Female1_B_CELL_NAIVE}{Expression levels from naive B cells.}
#'  \item{Female1_CD8_NAIVE}{Expression levels from naive CD8 T cells}
#'  \item{Female1_CD8_STIM}{Expression levels from stimulated CD8 T cells.}
#'  \item{Female1_TH2}{Expression levels from Th2 T cells.}
#'  \item{Female1_THSTAR}{Expression levels from Th1/17 T cells.}
#'  \item{Female1_TH17}{Expression levels from Th17 T cells.}
#'  \item{Female1_TH1}{Expression levels from Th1 T cells.}
#'  \item{Female1_TFH}{Expression levels from T follicular helper T cells.}
#'  \item{Female1_TREG_NAIVE}{Expression levels from naive regulatory T cells.}
#'  \item{Female1_TREG_MEM}{Expression levels from memory regulatory T cells.}
#'  \item{Female1_CD4_NAIVE}{Expression levels from naive CD4 T cells.}
#'  \item{Female1_CD4_STIM}{Expression levels from stimulated CD4 T cells.}
#'  \item{Female1_M2}{Expression levels from Non-classical Monocytes.}
#'  \item{Female1_NK}{Expression levels from Natural Killer cells.}
#'  \item{Female2_MONOCYTES}{Expression levels from Classical Monocytes.}
#'  \item{Female2_B_CELL_NAIVE}{Expression levels from naive B cells.}
#'  \item{Female2_CD8_NAIVE}{Expression levels from naive CD8 T cells}
#'  \item{Female2_CD8_STIM}{Expression levels from stimulated CD8 T cells.}
#'  \item{Female2_TH2}{Expression levels from Th2 T cells.}
#'  \item{Female2_THSTAR}{Expression levels from Th1/17 T cells.}
#'  \item{Female2_TH17}{Expression levels from Th17 T cells.}
#'  \item{Female2_TH1}{Expression levels from Th1 T cells.}
#'  \item{Female2_TFH}{Expression levels from T follicular helper T cells.}
#'  \item{Female2_TREG_NAIVE}{Expression levels from naive regulatory T cells.}
#'  \item{Female2_TREG_MEM}{Expression levels from memory regulatory T cells.}
#'  \item{Female2_CD4_NAIVE}{Expression levels from naive CD4 T cells.}
#'  \item{Female2_CD4_STIM}{Expression levels from stimulated CD4 T cells.}
#'  \item{Female2_M2}{Expression levels from Non-classical Monocytes.}
#'  \item{Female2_NK}{Expression levels from Natural Killer cells.}
#'  \item{Female3_MONOCYTES}{Expression levels from Classical Monocytes.}
#'  \item{Female3_B_CELL_NAIVE}{Expression levels from naive B cells.}
#'  \item{Female3_CD8_NAIVE}{Expression levels from naive CD8 T cells}
#'  \item{Female3_CD8_STIM}{Expression levels from stimulated CD8 T cells.}
#'  \item{Female3_TH2}{Expression levels from Th2 T cells.}
#'  \item{Female3_THSTAR}{Expression levels from Th1/17 T cells.}
#'  \item{Female3_TH17}{Expression levels from Th17 T cells.}
#'  \item{Female3_TH1}{Expression levels from Th1 T cells.}
#'  \item{Female3_TFH}{Expression levels from T follicular helper T cells.}
#'  \item{Female3_TREG_NAIVE}{Expression levels from naive regulatory T cells.}
#'  \item{Female3_TREG_MEM}{Expression levels from memory regulatory T cells.}
#'  \item{Female3_CD4_NAIVE}{Expression levels from naive CD4 T cells.}
#'  \item{Female3_CD4_STIM}{Expression levels from stimulated CD4 T cells.}
#'  \item{Female3_M2}{Expression levels from Non-classical Monocytes.}
#'  \item{Female3_NK}{Expression levels from Natural Killer cells.}
#'  \item{Female4_MONOCYTES}{Expression levels from Classical Monocytes.}
#'  \item{Female4_B_CELL_NAIVE}{Expression levels from naive B cells.}
#'  \item{Female4_CD8_NAIVE}{Expression levels from naive CD8 T cells}
#'  \item{Female4_CD8_STIM}{Expression levels from stimulated CD8 T cells.}
#'  \item{Female4_TH2}{Expression levels from Th2 T cells.}
#'  \item{Female4_THSTAR}{Expression levels from Th1/17 T cells.}
#'  \item{Female4_TH17}{Expression levels from Th17 T cells.}
#'  \item{Female4_TH1}{Expression levels from Th1 T cells.}
#'  \item{Female4_TFH}{Expression levels from T follicular helper T cells.}
#'  \item{Female4_TREG_NAIVE}{Expression levels from naive regulatory T cells.}
#'  \item{Female4_TREG_MEM}{Expression levels from memory regulatory T cells.}
#'  \item{Female4_CD4_NAIVE}{Expression levels from naive CD4 T cells.}
#'  \item{Female4_CD4_STIM}{Expression levels from stimulated CD4 T cells.}
#'  \item{Female4_M2}{Expression levels from Non-classical Monocytes.}
#'  \item{Female4_NK}{Expression levels from Natural Killer cells.}
#' }
"rnaSeq"

#' Protein coding genes found on the Y-chromosome. These genes are significantly
#'  expressed in the 'rnaSeq' dataset provided. They are to be used as the sample
#'  ground truth set in exploring the GECO package.
#'
"yChromSet"
