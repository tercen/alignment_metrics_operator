library(tercen)
library(dplyr)
library(seqinr)
library(tidyr)
library(msa)

options("tercen.workflowId" = "a77770c3923fad0ca99b77fa8905471d")
options("tercen.stepId"     = "7db7e18c-2668-49a6-b99d-91993d0d9942")

(ctx = tercenCtx())

sequence_type <- "protein"
if(!is.null(ctx$op.value('sequence_type'))) { sequence_type <- ctx$op.value('sequence_type') } 
substitutionMatrix <- "BLOSUM62"
if(!is.null(ctx$op.value('substitution_matrix'))) { substitutionMatrix <- ctx$op.value('substitution_matrix') }
gapvsgap <- NULL
if(!is.null(ctx$op.value('gap_vs_gap'))) { 
  gapvsgap <- as.numeric(ctx$op.value('gap_vs_gap'))
  if(is.na(gapvsgap)) gapvsgap <- NULL
}

letters <- ctx$select(ctx$colors[[1]])[[1]]

id_na_pos <- which(is.na(ctx$cselect(ctx$cnames[[1]]))) - 1
if(length(id_na_pos) == 0) id_na_pos <- -1

df <- ctx %>% select(.ri, .ci) %>%
  mutate(letter = letters) %>%
  filter(.ci != id_na_pos) %>% 
  spread(.ci, letter)

df[is.na(df)] <- "-"

aln <- apply(df[,-1], 1, function(x) paste0(x[!is.na(x)], collapse = ""))

if(sequence_type == "protein") AAset <- AAStringSet(aln)
if(sequence_type == "dna") AAset <- DNAStringSet(aln)
if(sequence_type == "rna") AAset <- RNAStringSet(aln)

conMat <- consensusMatrix(AAset)
conMat <- conMat[rownames(conMat) != "*", ]

## compute consensus scores using the BLOSUM62 matrix
data(list=substitutionMatrix)
submat <- eval(parse(text = substitutionMatrix))
ngaps <- conMat["-", ]
propgaps <- ngaps / length(aln)
cscore <- msaConservationScore(conMat, submat, gapVsGap = gapvsgap)

df_out <- data.frame(
  .ci = 1:length(cscore) - 1,
  conservation_score = cscore,
  n_gaps = ngaps,
  gap_proportion = propgaps
) %>%
  ctx$addNamespace() %>%
  ctx$save()



