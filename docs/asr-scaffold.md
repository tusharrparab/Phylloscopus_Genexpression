# ASR Scaffold

`run_asr_scaffold.py` is an optional downstream module. It is only meaningful once upstream loci are real and orthology has been vetted.

## Implemented Now

- per-locus FASTA parsing
- taxa-count checks
- optional `mafft` alignments
- optional `IQ-TREE` execution with `-asr`
- plan-script generation per locus

## Still Missing

- codon-aware alignment for coding loci
- locus occupancy thresholds
- tree conflict checks against a species tree
- partitioning strategy
- masking or filtering of poor alignments

## Honest Interpretation

- ASR on stub-generated sequences is not biologically interpretable
- ASR on poorly validated candidate-gene loci is risky
- ASR becomes reasonable only after real homologous sequences, robust alignments, and explicit orthology checks are in place
