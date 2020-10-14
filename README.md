# Alignment metrics operator

##### Description

`alignment_metrics` operator returns some metrics summarising a multiple sequence alignment.

##### Usage

Input projection|.
---|---
`row`        | numeric, input data, per cell 
`col`        | numeric, input data, per cell 
`colors`        | numeric, input data, per cell 
`y-axis`        | numeric, input data, per cell 

Properties|.
---|---
`sequence_type`        | factor, sequence type (protein, dna, rna)
`substitution_matrix`        | factor, substitution amtrix (default: BLOSUM62)
`gap_vs_gap`        | factor, gap vs gap score (default: NA)

Output relations|.
---|---
`conservation_score`        | numeric, conservation score at each position (column)
`n_gaps`        | numeric, number of gaps at each position (column)
`gap_proportion`        | numeric, conservation score at each position (column)

##### See also

[msa_operator](https://github.com/tercen/msa_operator)

[dist_alignment_operator](https://github.com/tercen/dist_alignment_operator)

[readfasta_operator](https://github.com/tercen/readfasta_operator)
