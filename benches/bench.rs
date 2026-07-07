use std::hint::black_box;
use std::io::Write;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_edger_diff_splice::{DiffSpliceArgs, GeneTest, diff_splice};

fn make_inputs(
    n_genes: usize,
    exons_per_gene: usize,
    n_samples: usize,
) -> (
    tempfile::TempDir,
    std::path::PathBuf,
    std::path::PathBuf,
    std::path::PathBuf,
) {
    let dir = tempfile::tempdir().unwrap();
    let counts = dir.path().join("counts.tsv");
    let design = dir.path().join("design.tsv");
    let genes = dir.path().join("genes.tsv");
    let mut c = std::fs::File::create(&counts).unwrap();
    let mut gf = std::fs::File::create(&genes).unwrap();
    write!(c, "exon").unwrap();
    for s in 0..n_samples {
        write!(c, "\tS{s}").unwrap();
    }
    writeln!(c).unwrap();
    writeln!(gf, "geneid").unwrap();
    let mut seed = 0x1234_5678u64;
    let mut rng = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    for g in 0..n_genes {
        for x in 0..exons_per_gene {
            write!(c, "G{g}.{x}").unwrap();
            for _ in 0..n_samples {
                write!(c, "\t{}", 50 + rng() % 400).unwrap();
            }
            writeln!(c).unwrap();
            writeln!(gf, "G{g}").unwrap();
        }
    }
    let mut d = std::fs::File::create(&design).unwrap();
    writeln!(d, "Intercept\tgroup").unwrap();
    for s in 0..n_samples {
        writeln!(d, "1\t{}", if s < n_samples / 2 { 0 } else { 1 }).unwrap();
    }
    (dir, counts, design, genes)
}

fn bench_diff_splice(c: &mut Criterion) {
    let (_dir, counts, design, genes) = make_inputs(4000, 7, 6);
    c.bench_function("diff_splice_28000x6", |b| {
        b.iter(|| {
            let mut out = Vec::new();
            diff_splice(
                &DiffSpliceArgs {
                    counts: black_box(&counts),
                    design: black_box(&design),
                    genes: black_box(&genes),
                    coef: Some(2),
                    contrast: None,
                    dispersion: 0.05,
                    dispersion_file: None,
                    norm_factors: None,
                    prior_count: 0.125,
                    test: GeneTest::Exon,
                    fdr: true,
                },
                &mut out,
            )
            .unwrap();
            black_box(out);
        })
    });
}

criterion_group!(benches, bench_diff_splice);
criterion_main!(benches);
