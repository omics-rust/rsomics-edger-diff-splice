# rsomics-edger-diff-splice

Differential exon usage from an exon-level negative-binomial GLM — a Rust port of
edgeR's `diffSpliceDGE()` + `topSpliceDGE()`. This is the count-based counterpart
to limma's `diffSplice` (which works on log-expression).

Given an exon-level count matrix (rows = exons, each annotated with a gene id), it
tests whether each exon's log-fold-change for the chosen coefficient differs from
its gene's overall log-fold-change. The gene log-fold-change is the slope of the
gene's summed-exon counts; the per-exon statistic is the drop in NB deviance when
that exon's tested coefficient is held at the gene slope.

| column (per-exon, `--test exon`) | meaning |
|---|---|
| `logFC` | natural-log fold change of the exon relative to its gene |
| `exon.LR` | likelihood-ratio statistic, 1 d.f. |
| `P.Value` | chi-square upper tail of `exon.LR` |
| `FDR` | Benjamini-Hochberg across exons (with `--fdr`) |

`--test gene` reports a gene-level LR test (`gene.LR` = sum of exon LRs, d.f. =
number of exons). `--test Simes` reports the Simes-combined exon p-value per gene.
Single-exon genes are not tested.

## Usage

```
rsomics-edger-diff-splice counts.tsv --design design.tsv --genes genes.tsv \
    [--coef N | --contrast c.tsv] [--test exon|Simes|gene] \
    [--dispersion D | --dispersion-file f.tsv] [--norm-factors f.tsv] [--fdr] [-o out.tsv]
```

- `counts.tsv` — header `exon<TAB>sample…`, one integer-count row per exon.
- `genes.tsv` — gene id per exon, one per line in the same row order (a one-token
  header is allowed).
- `design.tsv` — header of coefficient names, then one numeric row per sample.
- `--coef N` — 1-based design column to test (default: the last). `--contrast`
  takes a per-coefficient weight vector instead.
- `--dispersion` / `--dispersion-file` — common or per-exon NB dispersion (the
  exon dispersions from `estimateDisp`).

```
rsomics-edger-diff-splice counts.tsv --design design.tsv --genes genes.tsv --coef 2 --test Simes --fdr -o gene.tsv
```

## Origin

This crate is an independent, clean-room Rust reimplementation of edgeR's
`diffSpliceDGE` / `topSpliceDGE` based on:

- The published method: McCarthy DJ, Chen Y, Smyth GK, "Differential expression
  analysis of multifactor RNA-Seq experiments with respect to biological
  variation", *Nucleic Acids Research* 40(10):4288-4297, 2012.
  DOI: 10.1093/nar/gks042. Robinson MD, McCarthy DJ, Smyth GK, "edgeR: a
  Bioconductor package for differential expression analysis of digital gene
  expression data", *Bioinformatics* 26(1):139-140, 2010.
  DOI: 10.1093/bioinformatics/btp616.
- The public edgeR R-level documentation of `diffSpliceDGE` / `topSpliceDGE`
  (the per-exon-vs-gene log-fold-change construction, the gene LR = sum of exon
  LRs, the Simes combination), observed black-box.
- Black-box behaviour testing against the edgeR binary.

edgeR is GPL-licensed. **No edgeR source code was read or used as a reference
during implementation** — the algorithm was reconstructed from the published
method and the documented R-level behaviour, then validated value-exact against
the edgeR oracle. Test fixtures are independently generated.

License: MIT OR Apache-2.0.
Upstream credit: edgeR <https://bioconductor.org/packages/edgeR> (GPL).
