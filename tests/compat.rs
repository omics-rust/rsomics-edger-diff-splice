//! Differential compat against edgeR's diffSpliceDGE() + topSpliceDGE().
//!
//! `golden_matches_committed` always runs: our binary against the committed
//! R-captured goldens (exon / Simes / gene) on two fixtures — a well-conditioned
//! large-count set and a small-count set (with zeros) that exercises edgeR's
//! prior.count=0.125 shrinkage, the nexons-1 gene d.f., and the floored Simes.
//! `live_matches_edger` runs the real R upstream when an r-bioc/rs-edger Rscript
//! is found (RSOMICS_RSCRIPT overrides), else loud-skips.
//!
//! Every value matches edgeR to rel < 1e-6: coefficients come from the
//! prior-augmented fit, the likelihood ratios from the raw-count deviance, exactly
//! as diffSpliceDGE constructs them.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-edger-diff-splice"))
}

fn parse(text: &str) -> HashMap<String, Vec<f64>> {
    let mut lines = text.lines();
    lines.next();
    let mut rows = HashMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let key = f.next().unwrap().to_string();
        let vals: Vec<f64> = f.filter_map(|v| v.parse().ok()).collect();
        rows.insert(key, vals);
    }
    rows
}

/// `tol` is `(column, abs_tol, rel_tol)`; a value passes if it is within abs OR rel.
fn compare(expected: &str, got: &str, tol: &[(usize, f64, f64)]) {
    let e = parse(expected);
    let g = parse(got);
    assert_eq!(e.len(), g.len(), "row count mismatch");
    for (key, ev) in &e {
        let gv = g.get(key).unwrap_or_else(|| panic!("missing {key}"));
        for &(col, atol, rtol) in tol {
            let (a, b) = (ev[col], gv[col]);
            let ad = (a - b).abs();
            let rd = if a.abs() > 1e-12 { ad / a.abs() } else { ad };
            assert!(
                ad <= atol || rd <= rtol,
                "{key} col{col}: expected {a} got {b} (abs {ad:.2e} rel {rd:.2e})"
            );
        }
    }
}

fn run(prefix: &str, mode: &str, golden: &Path) -> String {
    let out = Command::new(bin())
        .arg(golden.join(format!("{prefix}counts.tsv")))
        .arg("--design")
        .arg(golden.join(format!("{prefix}design.tsv")))
        .arg("--genes")
        .arg(golden.join(format!("{prefix}genes.tsv")))
        .args(["--coef", "2", "--test", mode, "--fdr"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

/// NExons is an exact integer; every other column matches edgeR to rel < 1e-6,
/// with a small absolute floor for the near-zero (no-splicing) rows whose
/// likelihood ratio is a difference of two near-equal deviances.
const EXON_TOL: &[(usize, f64, f64)] = &[
    (0, 2e-6, 1e-6),
    (1, 2e-6, 1e-6),
    (2, 2e-6, 1e-6),
    (3, 2e-6, 1e-6),
];
const GENE_TOL: &[(usize, f64, f64)] = &[
    (0, 0.0, 0.0),
    (1, 2e-6, 1e-6),
    (2, 2e-6, 1e-6),
    (3, 2e-6, 1e-6),
];
const SIMES_TOL: &[(usize, f64, f64)] = &[(0, 0.0, 0.0), (1, 2e-6, 1e-6), (2, 2e-6, 1e-6)];

fn check_fixture(prefix: &str, golden: &Path) {
    compare(
        &run(prefix, "exon", golden),
        &std::fs::read_to_string(golden.join(format!("{prefix}expected_exon.tsv"))).unwrap(),
        EXON_TOL,
    );
    compare(
        &run(prefix, "Simes", golden),
        &std::fs::read_to_string(golden.join(format!("{prefix}expected_simes.tsv"))).unwrap(),
        SIMES_TOL,
    );
    compare(
        &run(prefix, "gene", golden),
        &std::fs::read_to_string(golden.join(format!("{prefix}expected_gene.tsv"))).unwrap(),
        GENE_TOL,
    );
}

#[test]
fn golden_matches_committed() {
    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    check_fixture("", &golden);
    check_fixture("small_", &golden);
}

/// edgeR errors on negative or non-finite counts ("Negative counts not allowed");
/// so does the parser, before any fitting. A garbage row must exit non-zero, never
/// silently produce a NaN VCF-equivalent.
#[test]
fn rejects_bad_counts() {
    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let dir = tempfile::tempdir().unwrap();
    for (name, bad_row) in [
        ("neg", "G1.2\t-5\t3\t2\t4\t1\t2"),
        ("inf", "G1.2\tinf\t3\t2\t4\t1\t2"),
        ("nan", "G1.2\tnan\t3\t2\t4\t1\t2"),
    ] {
        let counts = dir.path().join(format!("{name}.tsv"));
        std::fs::write(
            &counts,
            format!("exon\tS1\tS2\tS3\tS4\tS5\tS6\nG1.1\t1\t2\t3\t4\t5\t6\n{bad_row}\nG1.3\t2\t2\t2\t2\t2\t2\n"),
        )
        .unwrap();
        let out = Command::new(bin())
            .arg(&counts)
            .arg("--design")
            .arg(golden.join("small_design.tsv"))
            .arg("--genes")
            .arg(golden.join("small_genes.tsv"))
            .output()
            .unwrap();
        assert!(
            !out.status.success(),
            "{name}: expected non-zero exit on a bad count"
        );
    }
}

fn find_rscript() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("RSOMICS_RSCRIPT") {
        let p = PathBuf::from(p);
        if p.exists() {
            return Some(p);
        }
    }
    let home = std::env::var("HOME").ok()?;
    for cand in [
        "miniconda3/envs/r-bioc/bin/Rscript",
        "miniforge3/envs/rs-edger/bin/Rscript",
    ] {
        let p = PathBuf::from(&home).join(cand);
        if p.exists() {
            return Some(p);
        }
    }
    None
}

#[test]
fn live_matches_edger() {
    let Some(rscript) = find_rscript() else {
        eprintln!("SKIP live_matches_edger: no r-bioc/rs-edger Rscript (set RSOMICS_RSCRIPT)");
        return;
    };
    let has_edger = Command::new(&rscript)
        .args(["-e", "suppressMessages(library(edgeR))"])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    if !has_edger {
        eprintln!("SKIP live_matches_edger: Rscript lacks edgeR");
        return;
    }

    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    for prefix in ["", "small_"] {
        let counts = golden.join(format!("{prefix}counts.tsv"));
        let design = golden.join(format!("{prefix}design.tsv"));
        let genes = golden.join(format!("{prefix}genes.tsv"));
        let scratch =
            std::env::temp_dir().join(format!("rsomics-edger-diff-splice-compat-{prefix}x"));
        std::fs::create_dir_all(&scratch).unwrap();
        let r_exon = scratch.join("r_exon.tsv");
        let r_simes = scratch.join("r_simes.tsv");
        let r_gene = scratch.join("r_gene.tsv");

        // gene.df.test holds nexons-1; topSpliceDGE reports the actual exon count,
        // so NExons = gene.df.test + 1.
        let script = format!(
            r#"
suppressMessages(library(edgeR))
options(digits=17)
counts <- as.matrix(read.delim("{c}", row.names=1))
design <- as.matrix(read.delim("{d}", check.names=FALSE))
geneid <- as.character(read.delim("{g}")[[1]])
exonid <- ave(geneid, geneid, FUN=seq_along)
fit <- glmFit(DGEList(counts=counts), design, dispersion=0.05)
ds <- diffSpliceDGE(fit, geneid=geneid, exonid=exonid, coef=2, verbose=FALSE)
ek <- paste0(ds$genes$GeneID, ".", ds$genes$ExonID)
ex <- data.frame(ExonID=ek, GeneID=ds$genes$GeneID, logFC=ds$coefficients,
  exon.LR=ds$exon.LR, P.Value=ds$exon.p.value, FDR=p.adjust(ds$exon.p.value,"BH"))
write.table(ex, "{oe}", sep="\t", quote=FALSE, row.names=FALSE)
gg <- rownames(ds$gene.df.test)
nex <- as.integer(ds$gene.df.test[,1]) + 1L
si <- data.frame(GeneID=gg, NExons=nex,
  P.Value=ds$gene.Simes.p.value, FDR=p.adjust(ds$gene.Simes.p.value,"BH"))
write.table(si, "{os}", sep="\t", quote=FALSE, row.names=FALSE)
ge <- data.frame(GeneID=gg, NExons=nex,
  gene.LR=as.numeric(ds$gene.LR), P.Value=as.numeric(ds$gene.p.value),
  FDR=p.adjust(as.numeric(ds$gene.p.value),"BH"))
write.table(ge, "{og}", sep="\t", quote=FALSE, row.names=FALSE)
"#,
            c = counts.display(),
            d = design.display(),
            g = genes.display(),
            oe = r_exon.display(),
            os = r_simes.display(),
            og = r_gene.display(),
        );
        let st = Command::new(&rscript)
            .args(["-e", &script])
            .status()
            .unwrap();
        assert!(st.success(), "R edgeR run failed");

        compare(
            &run(prefix, "exon", &golden),
            &std::fs::read_to_string(&r_exon).unwrap(),
            EXON_TOL,
        );
        compare(
            &run(prefix, "Simes", &golden),
            &std::fs::read_to_string(&r_simes).unwrap(),
            SIMES_TOL,
        );
        compare(
            &run(prefix, "gene", &golden),
            &std::fs::read_to_string(&r_gene).unwrap(),
            GENE_TOL,
        );
    }
}
