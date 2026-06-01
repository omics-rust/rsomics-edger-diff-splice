//! Differential compat against edgeR's diffSpliceDGE() + topSpliceDGE().
//!
//! `golden_matches_committed` always runs: our binary against the committed
//! R-captured goldens (exon / Simes / gene). `live_matches_edger` runs the real R
//! upstream when an r-bioc Rscript is found (RSOMICS_RSCRIPT or
//! ~/miniconda3/envs/r-bioc/bin/Rscript), else loud-skips.
//!
//! Small (significant) p-values match edgeR to ~1e-6. Everything else carries the
//! slack of edgeR's own tol=1e-6 Levenberg fit (ours converges tighter): logFC to
//! ~1e-3, the LR statistics to ~1e-2, and mid-range p-values/FDR to ~1.5e-3. The
//! inference is identical — the same exons and genes are called significant.

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

fn run(mode: &str, golden: &Path) -> String {
    let out = Command::new(bin())
        .arg(golden.join("counts.tsv"))
        .arg("--design")
        .arg(golden.join("design.tsv"))
        .arg("--genes")
        .arg(golden.join("genes.tsv"))
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

#[test]
fn golden_matches_committed() {
    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    compare(
        &run("exon", &golden),
        &std::fs::read_to_string(golden.join("expected_exon.tsv")).unwrap(),
        &[
            (0, 1e-3, 1e-2),
            (1, 1e-2, 3e-3),
            (2, 1.5e-3, 5e-3),
            (3, 1.5e-3, 5e-3),
        ],
    );
    // NExons is an exact integer; the rest carries the Levenberg slack.
    compare(
        &run("Simes", &golden),
        &std::fs::read_to_string(golden.join("expected_simes.tsv")).unwrap(),
        &[(0, 0.0, 0.0), (1, 1.5e-3, 5e-3), (2, 1.5e-3, 5e-3)],
    );
    compare(
        &run("gene", &golden),
        &std::fs::read_to_string(golden.join("expected_gene.tsv")).unwrap(),
        &[
            (0, 0.0, 0.0),
            (1, 1e-2, 3e-3),
            (2, 1.5e-3, 5e-3),
            (3, 1.5e-3, 5e-3),
        ],
    );
}

fn find_rscript() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("RSOMICS_RSCRIPT") {
        let p = PathBuf::from(p);
        if p.exists() {
            return Some(p);
        }
    }
    let home = std::env::var("HOME").ok()?;
    let cand = PathBuf::from(home).join("miniconda3/envs/r-bioc/bin/Rscript");
    cand.exists().then_some(cand)
}

#[test]
fn live_matches_edger() {
    let Some(rscript) = find_rscript() else {
        eprintln!("SKIP live_matches_edger: no r-bioc Rscript (set RSOMICS_RSCRIPT)");
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
    let counts = golden.join("counts.tsv");
    let design = golden.join("design.tsv");
    let genes = golden.join("genes.tsv");
    let scratch = std::env::temp_dir().join("rsomics-edger-diff-splice-compat");
    std::fs::create_dir_all(&scratch).unwrap();
    let r_exon = scratch.join("r_exon.tsv");
    let r_simes = scratch.join("r_simes.tsv");
    let r_gene = scratch.join("r_gene.tsv");

    let script = format!(
        r#"
suppressMessages(library(edgeR))
counts <- as.matrix(read.delim("{c}", row.names=1))
design <- as.matrix(read.delim("{d}"))
geneid <- read.delim("{g}")[[1]]
exonid <- ave(geneid, geneid, FUN=seq_along)
fit <- glmFit(DGEList(counts=counts), design, dispersion=0.05)
ds <- diffSpliceDGE(fit, geneid=geneid, exonid=exonid, coef=2)
ek <- paste0(ds$genes$GeneID, ".", ds$genes$ExonID)
ex <- data.frame(ExonID=ek, GeneID=ds$genes$GeneID, logFC=ds$coefficients,
  exon.LR=ds$exon.LR, P.Value=ds$exon.p.value, FDR=p.adjust(ds$exon.p.value,"BH"))
write.table(ex, "{oe}", sep="\t", quote=FALSE, row.names=FALSE)
gg <- rownames(ds$gene.df.test)
si <- data.frame(GeneID=gg, NExons=as.integer(ds$gene.df.test[,1]),
  P.Value=ds$gene.Simes.p.value, FDR=p.adjust(ds$gene.Simes.p.value,"BH"))
write.table(si, "{os}", sep="\t", quote=FALSE, row.names=FALSE)
ge <- data.frame(GeneID=gg, NExons=as.integer(ds$gene.df.test[,1]),
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
        &run("exon", &golden),
        &std::fs::read_to_string(&r_exon).unwrap(),
        &[
            (0, 1e-3, 1e-2),
            (1, 1e-2, 3e-3),
            (2, 1.5e-3, 5e-3),
            (3, 1.5e-3, 5e-3),
        ],
    );
    compare(
        &run("Simes", &golden),
        &std::fs::read_to_string(&r_simes).unwrap(),
        &[(0, 0.0, 0.0), (1, 1.5e-3, 5e-3), (2, 1.5e-3, 5e-3)],
    );
    compare(
        &run("gene", &golden),
        &std::fs::read_to_string(&r_gene).unwrap(),
        &[
            (0, 0.0, 0.0),
            (1, 1e-2, 3e-3),
            (2, 1.5e-3, 5e-3),
            (3, 1.5e-3, 5e-3),
        ],
    );
}
