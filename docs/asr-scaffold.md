# ASR Scaffold

The ASR scaffold is implemented in [run_asr_scaffold.py](/Users/sam/Documents/New%20project/bin/run_asr_scaffold.py).

It expects per-locus FASTA files, such as the merged ortholog FASTA bundles under `results/merged/ortholog_sequences/`, and then:

1. counts taxa and locus lengths
2. aligns loci with `mafft --auto` when `mafft` is installed
3. falls back to the input FASTA when alignment is unavailable or fails
4. plans or runs ancestral reconstruction with `IQ-TREE` using `-asr`
5. writes per-locus status tables and runnable shell scripts

The `IQ-TREE` command syntax follows the official command reference, where `-asr` writes ancestral states and `-te` can be used with a user-supplied tree ([IQ-TREE command reference](https://iqtree.github.io/doc/Command-Reference)).

## Typical Usage

```bash
python3 bin/run_asr_scaffold.py \
  --sequence-dir results/merged/ortholog_sequences \
  --species-tree examples/asr/phylloscopus_example_tree.nwk \
  --outdir results/asr \
  --mode scaffold
```

## Outputs

- `asr/locus_manifest.tsv`
- `asr/asr_plan.tsv`
- `asr/asr_summary.tsv`
- `asr/plans/<locus>/run_asr.sh`

When `IQ-TREE` is not installed, the scaffold still emits `run_asr.sh` plans for each eligible locus instead of failing the entire downstream stage.
