//! edgeR diffSpliceDGE + topSpliceDGE: differential exon usage from an exon-level
//! NB-GLM fit. For the tested coefficient, each exon's log-fold-change is compared
//! to its gene's overall log-fold-change (the slope of the gene's summed counts).
//! Method: McCarthy, Chen & Smyth (2012) NAR 40:4288-4297 (NB GLM), with the
//! diffSpliceDGE construction reconstructed clean-room from the documented method
//! and black-box behaviour. Per-exon LR = drop in NB deviance when the tested
//! coefficient is held at the gene slope; gene LR = sum over exons (df = nExons − 1);
//! Simes combines exon p-values per gene. Counts carry edgeR's prior.count=0.125.

mod special;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

const MAXIT: usize = 200;
const TOL: f64 = 1e-11;

/// edgeR's `diffSpliceDGE` fits the gene-level betabar with `glmFit`'s default
/// dispersion (0.05), never the exon dispersions passed in — replicated here so
/// the reported log-fold-changes are value-exact even when exon dispersions differ.
const GENE_FIT_DISPERSION: f64 = 0.05;

pub struct Matrix {
    pub header: String,
    pub exons: Vec<String>,
    pub counts: Vec<f64>,
    pub n_samples: usize,
}

impl Matrix {
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
        let mut lines = BufReader::new(file).lines();
        let header = lines
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("empty count matrix".into()))?
            .map_err(RsomicsError::Io)?;
        let n_samples = header.split('\t').count() - 1;
        if n_samples == 0 {
            return Err(RsomicsError::InvalidInput(
                "count matrix has no sample columns".into(),
            ));
        }
        let mut exons = Vec::new();
        let mut counts = Vec::new();
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.is_empty() {
                continue;
            }
            let mut fields = line.split('\t');
            let exon = fields
                .next()
                .ok_or_else(|| RsomicsError::InvalidInput("row without an exon id".into()))?;
            exons.push(exon.to_string());
            let before = counts.len();
            for f in fields {
                let c = f.parse::<f64>().map_err(|_| {
                    RsomicsError::InvalidInput(format!("non-numeric count '{f}' for exon {exon}"))
                })?;
                if !c.is_finite() {
                    return Err(RsomicsError::InvalidInput(format!(
                        "non-finite count '{f}' for exon {exon}"
                    )));
                }
                if c < 0.0 {
                    return Err(RsomicsError::InvalidInput(format!(
                        "negative count {c} for exon {exon}"
                    )));
                }
                counts.push(c);
            }
            if counts.len() - before != n_samples {
                return Err(RsomicsError::InvalidInput(format!(
                    "exon {exon}: {} values, header has {n_samples} samples",
                    counts.len() - before
                )));
            }
        }
        Ok(Self {
            header,
            exons,
            counts,
            n_samples,
        })
    }

    pub fn n_exons(&self) -> usize {
        self.exons.len()
    }
    fn row(&self, e: usize) -> &[f64] {
        &self.counts[e * self.n_samples..(e + 1) * self.n_samples]
    }
}

pub struct Design {
    pub data: Vec<f64>,
    pub n_samples: usize,
    pub n_coef: usize,
    pub coef_names: Vec<String>,
}

impl Design {
    fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
        let mut lines = BufReader::new(file).lines();
        let header = lines
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("empty design matrix".into()))?
            .map_err(RsomicsError::Io)?;
        let coef_names: Vec<String> = header.split('\t').map(str::to_string).collect();
        let n_coef = coef_names.len();
        let mut data = Vec::new();
        let mut n_samples = 0;
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.is_empty() {
                continue;
            }
            let before = data.len();
            for f in line.split('\t') {
                data.push(f.parse::<f64>().map_err(|_| {
                    RsomicsError::InvalidInput(format!("non-numeric design value '{f}'"))
                })?);
            }
            if data.len() - before != n_coef {
                return Err(RsomicsError::InvalidInput(format!(
                    "design row {n_samples}: {} values, header has {n_coef} columns",
                    data.len() - before
                )));
            }
            n_samples += 1;
        }
        Ok(Self {
            data,
            n_samples,
            n_coef,
            coef_names,
        })
    }

    fn row(&self, s: usize) -> &[f64] {
        &self.data[s * self.n_coef..(s + 1) * self.n_coef]
    }
}

fn load_genes(path: &Path, n_exons: usize) -> Result<Vec<String>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut lines = BufReader::new(file).lines();
    let first = lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("empty genes file".into()))?
        .map_err(RsomicsError::Io)?;
    let mut genes = Vec::with_capacity(n_exons);
    // A single-token header ("geneid") is consumed; a value-looking first line is kept.
    let header_like = first.split('\t').count() == 1 && first.parse::<f64>().is_err();
    let take = |g: &str| g.split('\t').next().unwrap_or(g).to_string();
    if !header_like {
        genes.push(take(&first));
    }
    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        genes.push(take(&line));
    }
    if genes.len() != n_exons {
        return Err(RsomicsError::InvalidInput(format!(
            "{} gene ids for {n_exons} exons",
            genes.len()
        )));
    }
    Ok(genes)
}

fn load_one_per_line(path: &Path, n: usize, what: &str) -> Result<Vec<f64>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut v = Vec::with_capacity(n);
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let val = line.rsplit('\t').next().unwrap_or(line);
        if let Ok(x) = val.parse::<f64>() {
            v.push(x);
        }
    }
    if v.len() != n {
        return Err(RsomicsError::InvalidInput(format!(
            "{} {what} for {n} expected",
            v.len()
        )));
    }
    Ok(v)
}

/// NB deviance, edgeR nbinomDeviance.
fn nb_deviance(y: &[f64], mu: &[f64], dispersion: f64) -> f64 {
    let r = 1.0 / dispersion;
    let mut dev = 0.0;
    for (&yi, &mui) in y.iter().zip(mu) {
        let term_y = if yi > 0.0 { yi * (yi / mui).ln() } else { 0.0 };
        dev += term_y - (yi + r) * ((yi + r) / (mui + r)).ln();
    }
    2.0 * dev
}

struct Fit {
    beta: Vec<f64>,
    deviance: f64,
}

/// Per-exon precompute: the prior-shrunk tested coefficient (for the reported
/// logFC) alongside the raw-count fit (its beta seeds the null refit, its deviance
/// is the unshrunk full deviance for the LR).
struct ExonFit {
    beta_aug_tested: f64,
    beta_raw: Vec<f64>,
    dev_raw: f64,
}

/// IRLS NB GLM, log link, per-sample offset, Levenberg ridge — edgeR mglmLevenberg.
/// `fixed` pins coefficient indices to their starting value (their score rows are
/// zeroed so they never move), used for the constrained per-exon null fit.
fn fit_nb_glm(
    y: &[f64],
    x: &[f64],
    offset: &[f64],
    dispersion: f64,
    start: &[f64],
    fixed: Option<usize>,
) -> Fit {
    let n = offset.len();
    let p = start.len();
    let mut beta = start.to_vec();
    let mut mu = vec![0.0f64; n];
    eta_mu(x, offset, &beta, n, p, &mut mu);
    let mut dev = nb_deviance(y, &mu, dispersion);
    let mut lambda = 0.0f64;
    let mut xtwx = vec![0.0f64; p * p];
    let mut xtr = vec![0.0f64; p];
    let mut a = vec![0.0f64; p * p];
    let mut rhs = vec![0.0f64; p];
    let mut step = vec![0.0f64; p];
    let mut trial = vec![0.0f64; p];
    let mut mu_t = vec![0.0f64; n];

    for _ in 0..MAXIT {
        xtwx.iter_mut().for_each(|v| *v = 0.0);
        xtr.iter_mut().for_each(|v| *v = 0.0);
        for s in 0..n {
            let xr = &x[s * p..s * p + p];
            let mui = mu[s];
            let denom = 1.0 + dispersion * mui;
            let w = mui / denom;
            let resid = (y[s] - mui) / denom;
            for j in 0..p {
                xtr[j] += xr[j] * resid;
                let xjw = xr[j] * w;
                for k in 0..p {
                    xtwx[j * p + k] += xjw * xr[k];
                }
            }
        }
        if let Some(f) = fixed {
            xtr[f] = 0.0;
            for k in 0..p {
                xtwx[f * p + k] = 0.0;
                xtwx[k * p + f] = 0.0;
            }
            xtwx[f * p + f] = 1.0;
        }
        let mut accepted = false;
        for _ in 0..30 {
            a.copy_from_slice(&xtwx);
            for d in 0..p {
                if Some(d) != fixed {
                    a[d * p + d] += lambda * xtwx[d * p + d].max(1e-6);
                }
            }
            rhs.copy_from_slice(&xtr);
            if !solve(&mut a, &mut rhs, &mut step, p) {
                lambda = if lambda == 0.0 { 1.0 } else { lambda * 2.0 };
                continue;
            }
            for j in 0..p {
                trial[j] = beta[j] + step[j];
            }
            eta_mu(x, offset, &trial, n, p, &mut mu_t);
            let dev_t = nb_deviance(y, &mu_t, dispersion);
            if dev_t <= dev + 1e-8 * (1.0 + dev.abs()) {
                let max_step = step.iter().fold(0.0f64, |m, s| m.max(s.abs()));
                beta.copy_from_slice(&trial);
                mu.copy_from_slice(&mu_t);
                dev = dev_t;
                lambda *= 0.5;
                accepted = true;
                if max_step < TOL {
                    return Fit {
                        beta,
                        deviance: dev,
                    };
                }
                break;
            }
            lambda = if lambda == 0.0 { 1.0 } else { lambda * 4.0 };
        }
        if !accepted {
            break;
        }
    }
    Fit {
        beta,
        deviance: dev,
    }
}

fn eta_mu(x: &[f64], offset: &[f64], beta: &[f64], n: usize, p: usize, mu: &mut [f64]) {
    for s in 0..n {
        let xr = &x[s * p..s * p + p];
        let mut eta = offset[s];
        for (&xv, &b) in xr.iter().zip(beta) {
            eta += xv * b;
        }
        mu[s] = eta.exp();
    }
}

/// One-group NB fit: the intercept-only starting mean, in natural-log.
fn one_group_start(y: &[f64], offset: &[f64]) -> f64 {
    let total: f64 = y.iter().sum();
    let n = y.len() as f64;
    let mean_off = offset.iter().sum::<f64>() / n;
    if total <= 0.0 {
        return -1e8;
    }
    (total / n).ln() - mean_off
}

/// edgeR start.method="null": fit a common mean b0, project onto the design as
/// b0·(XᵀX)⁻¹Xᵀ1 so every gene/exon starts from a sound intercept-only point.
fn start_direction(x: &[f64], n: usize, p: usize) -> Vec<f64> {
    let mut xtx = vec![0.0f64; p * p];
    let mut xt1 = vec![0.0f64; p];
    for s in 0..n {
        let xr = &x[s * p..s * p + p];
        for j in 0..p {
            xt1[j] += xr[j];
            for k in 0..p {
                xtx[j * p + k] += xr[j] * xr[k];
            }
        }
    }
    let mut d = vec![0.0f64; p];
    solve(&mut xtx, &mut xt1, &mut d, p);
    d
}

/// Solve A x = b (row-major p×p), Gaussian elimination with partial pivoting.
fn solve(a: &mut [f64], rhs: &mut [f64], x: &mut [f64], p: usize) -> bool {
    for col in 0..p {
        let mut piv = col;
        let mut best = a[col * p + col].abs();
        for r in (col + 1)..p {
            let v = a[r * p + col].abs();
            if v > best {
                best = v;
                piv = r;
            }
        }
        if best < 1e-13 {
            return false;
        }
        if piv != col {
            for k in 0..p {
                a.swap(col * p + k, piv * p + k);
            }
            rhs.swap(col, piv);
        }
        let d = a[col * p + col];
        for r in (col + 1)..p {
            let f = a[r * p + col] / d;
            if f == 0.0 {
                continue;
            }
            for k in col..p {
                a[r * p + k] -= f * a[col * p + k];
            }
            rhs[r] -= f * rhs[col];
        }
    }
    for col in (0..p).rev() {
        let mut s = rhs[col];
        for k in (col + 1)..p {
            s -= a[col * p + k] * x[k];
        }
        x[col] = s / a[col * p + col];
    }
    true
}

fn bh_fdr(pvals: &[f64]) -> Vec<f64> {
    let n = pvals.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| pvals[b].partial_cmp(&pvals[a]).unwrap());
    let mut adj = vec![0.0f64; n];
    let mut cummin = f64::INFINITY;
    for (rank, &i) in order.iter().enumerate() {
        let m = n - rank;
        let v = (pvals[i] * n as f64 / m as f64).min(1.0);
        cummin = cummin.min(v);
        adj[i] = cummin;
    }
    adj
}

/// edgeR's Simes combination of a gene's exon p-values: the minimum over sorted
/// p of p_(r)·max((n-1)/r, 1). Differs from textbook Simes (n·p_(r)/r) — edgeR
/// uses n-1 and floors the multiplier at 1 so the largest p-value enters as-is.
fn simes(pvals: &[f64]) -> f64 {
    let mut p: Vec<f64> = pvals.to_vec();
    p.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = p.len() as f64;
    p.iter()
        .enumerate()
        .map(|(i, &pi)| {
            let r = i as f64 + 1.0;
            pi * ((n - 1.0) / r).max(1.0)
        })
        .fold(f64::INFINITY, f64::min)
        .min(1.0)
}

#[derive(Clone, Copy)]
pub enum GeneTest {
    Exon,
    Simes,
    Gene,
}

pub struct DiffSpliceArgs<'a> {
    pub counts: &'a Path,
    pub design: &'a Path,
    pub genes: &'a Path,
    pub coef: Option<usize>,
    pub contrast: Option<&'a Path>,
    pub dispersion: f64,
    pub dispersion_file: Option<&'a Path>,
    pub norm_factors: Option<&'a Path>,
    pub prior_count: f64,
    pub test: GeneTest,
    pub fdr: bool,
}

struct ExonResult {
    coef: f64,
    lr: f64,
    pval: f64,
}

pub fn diff_splice(args: &DiffSpliceArgs, output: &mut dyn Write) -> Result<u64> {
    let m = Matrix::load(args.counts)?;
    let design = Design::load(args.design)?;
    if design.n_samples != m.n_samples {
        return Err(RsomicsError::InvalidInput(format!(
            "design has {} rows but matrix has {} samples",
            design.n_samples, m.n_samples
        )));
    }
    let genes = load_genes(args.genes, m.n_exons())?;
    let p = design.n_coef;

    let norm_factors = match args.norm_factors {
        Some(path) => load_one_per_line(path, m.n_samples, "norm factors")?,
        None => vec![1.0; m.n_samples],
    };

    let mut lib = vec![0.0f64; m.n_samples];
    for row in m.counts.chunks_exact(m.n_samples) {
        for (s, &c) in lib.iter_mut().zip(row) {
            *s += c;
        }
    }
    let offset: Vec<f64> = lib
        .iter()
        .zip(&norm_factors)
        .map(|(&l, &f)| (l * f).ln())
        .collect();

    // A contrast is tested by rotating the design so the contrast becomes the
    // first coefficient (Householder QR); a coef is tested in place.
    let (xdata, tested_col, contrast_scale) = match (args.coef, args.contrast) {
        (Some(_), Some(_)) => {
            return Err(RsomicsError::InvalidInput(
                "give --coef or --contrast, not both".into(),
            ));
        }
        (Some(c), None) => {
            if c == 0 || c > p {
                return Err(RsomicsError::InvalidInput(format!(
                    "--coef {c} out of range 1..={p}"
                )));
            }
            (design.data.clone(), c - 1, 1.0)
        }
        (None, Some(cpath)) => {
            let contrast = load_one_per_line(cpath, p, "contrast weights")?;
            let q = householder_basis(&contrast);
            let mut xr = vec![0.0f64; design.n_samples * p];
            for s in 0..design.n_samples {
                let row = design.row(s);
                for j in 0..p {
                    let mut v = 0.0;
                    for (k, &rk) in row.iter().enumerate() {
                        v += rk * q[k * p + j];
                    }
                    xr[s * p + j] = v;
                }
            }
            let norm = contrast.iter().map(|x| x * x).sum::<f64>().sqrt();
            (xr, 0usize, norm)
        }
        (None, None) => (design.data.clone(), p - 1, 1.0),
    };

    let dispersions = match args.dispersion_file {
        Some(path) => load_one_per_line(path, m.n_exons(), "dispersions")?,
        None => vec![args.dispersion; m.n_exons()],
    };

    let n = m.n_samples;

    // edgeR addPriorCount: a per-sample prior proportional to library size is added
    // to every fitted response and the offset grown to match. prior_count defaults
    // to edgeR's 0.125; without it small-count exon/gene fits diverge from edgeR.
    let libsize: Vec<f64> = offset.iter().map(|&o| o.exp()).collect();
    let mean_ls = libsize.iter().sum::<f64>() / n as f64;
    let prior_scaled: Vec<f64> = libsize
        .iter()
        .map(|&l| args.prior_count * l / mean_ls)
        .collect();
    let offset_adj: Vec<f64> = libsize
        .iter()
        .zip(&prior_scaled)
        .map(|(&l, &pc)| (l + 2.0 * pc).ln())
        .collect();

    // Group exons by gene in first-seen order; single-exon genes are not tested.
    let mut gene_order: Vec<String> = Vec::new();
    let mut gene_exons: HashMap<String, Vec<usize>> = HashMap::new();
    for (e, g) in genes.iter().enumerate() {
        gene_exons
            .entry(g.clone())
            .or_insert_with(|| {
                gene_order.push(g.clone());
                Vec::new()
            })
            .push(e);
    }

    let dir = start_direction(&xdata, n, p);
    // edgeR reports a prior-shrunk log-fold-change but an unshrunk (raw) likelihood
    // ratio: coefficients come from the prior.count-augmented fit, the LR from the
    // raw-count deviance. Each exon therefore needs both fits.
    let per_exon = |e: usize| -> ExonFit {
        let raw = m.row(e);
        let y_aug: Vec<f64> = raw
            .iter()
            .zip(&prior_scaled)
            .map(|(&c, &pc)| c + pc)
            .collect();
        let start_a: Vec<f64> = {
            let b0 = one_group_start(&y_aug, &offset_adj);
            dir.iter().map(|&d| b0 * d).collect()
        };
        let aug = fit_nb_glm(&y_aug, &xdata, &offset_adj, dispersions[e], &start_a, None);

        let y_raw = raw.to_vec();
        let start_r: Vec<f64> = {
            let b0 = one_group_start(&y_raw, &offset);
            dir.iter().map(|&d| b0 * d).collect()
        };
        let raw_fit = fit_nb_glm(&y_raw, &xdata, &offset, dispersions[e], &start_r, None);

        ExonFit {
            beta_aug_tested: aug.beta[tested_col],
            beta_raw: raw_fit.beta,
            dev_raw: raw_fit.deviance,
        }
    };

    let exon_fits: Vec<ExonFit> = if rayon::current_num_threads() > 1 {
        use rayon::prelude::*;
        (0..m.n_exons()).into_par_iter().map(per_exon).collect()
    } else {
        (0..m.n_exons()).map(per_exon).collect()
    };

    let mut results: HashMap<usize, ExonResult> = HashMap::new();
    let mut gene_lr: HashMap<String, (f64, usize)> = HashMap::new();
    let mut gene_simes: HashMap<String, f64> = HashMap::new();
    let mut n_tested_exons = 0usize;

    for g in &gene_order {
        let exons = &gene_exons[g];
        if exons.len() < 2 {
            continue;
        }
        // Gene log-fold-change: slope of the gene's summed exon counts.
        let mut total = vec![0.0f64; n];
        for &e in exons {
            for (s, &c) in total.iter_mut().zip(m.row(e)) {
                *s += c;
            }
        }
        let total_adj: Vec<f64> = total
            .iter()
            .zip(&prior_scaled)
            .map(|(&t, &pc)| t + pc)
            .collect();
        let b0 = one_group_start(&total_adj, &offset_adj);
        let start: Vec<f64> = dir.iter().map(|&d| b0 * d).collect();
        let gene_fit = fit_nb_glm(
            &total_adj,
            &xdata,
            &offset_adj,
            GENE_FIT_DISPERSION,
            &start,
            None,
        );
        let beta_gene = gene_fit.beta[tested_col];

        let mut g_lr = 0.0;
        let mut exon_pvals = Vec::with_capacity(exons.len());
        for &e in exons {
            let y = m.row(e);
            let disp = dispersions[e];
            let coef = exon_fits[e].beta_aug_tested - beta_gene;
            // Null deviance: the raw-count exon refit with the tested coef pinned to
            // the (shrunk) gene slope. The LR is unshrunk, so this fit sees no prior.
            let mut start = exon_fits[e].beta_raw.clone();
            start[tested_col] = beta_gene;
            let null_fit = fit_nb_glm(y, &xdata, &offset, disp, &start, Some(tested_col));
            let lr = (null_fit.deviance - exon_fits[e].dev_raw).max(0.0);
            let pval = special::pchisq_upper(lr, 1.0);
            g_lr += lr;
            exon_pvals.push(pval);
            results.insert(
                e,
                ExonResult {
                    coef: coef * contrast_scale,
                    lr,
                    pval,
                },
            );
            n_tested_exons += 1;
        }
        gene_lr.insert(g.clone(), (g_lr, exons.len()));
        gene_simes.insert(g.clone(), simes(&exon_pvals));
    }

    write_output(
        output,
        args,
        &m,
        &genes,
        &gene_order,
        &results,
        &gene_lr,
        &gene_simes,
        n_tested_exons,
    )
}

#[allow(clippy::too_many_arguments)]
fn write_output(
    output: &mut dyn Write,
    args: &DiffSpliceArgs,
    m: &Matrix,
    genes: &[String],
    gene_order: &[String],
    results: &HashMap<usize, ExonResult>,
    gene_lr: &HashMap<String, (f64, usize)>,
    gene_simes: &HashMap<String, f64>,
    n_tested_exons: usize,
) -> Result<u64> {
    let exon_col = m.header.split('\t').next().unwrap_or("ExonID");
    match args.test {
        GeneTest::Exon => {
            let mut idx: Vec<usize> = (0..m.n_exons())
                .filter(|e| results.contains_key(e))
                .collect();
            let pvals: Vec<f64> = idx.iter().map(|&e| results[&e].pval).collect();
            let fdr = bh_fdr(&pvals);
            let fdr_map: HashMap<usize, f64> = idx.iter().copied().zip(fdr).collect();
            let mut header = format!("{exon_col}\tGeneID\tlogFC\texon.LR\tP.Value");
            if args.fdr {
                header.push_str("\tFDR");
            }
            writeln!(output, "{header}").map_err(RsomicsError::Io)?;
            idx.sort_unstable();
            for e in idx {
                let r = &results[&e];
                write!(
                    output,
                    "{}\t{}\t{:.7e}\t{:.7e}\t{:.7e}",
                    m.exons[e], genes[e], r.coef, r.lr, r.pval
                )
                .map_err(RsomicsError::Io)?;
                if args.fdr {
                    write!(output, "\t{:.7e}", fdr_map[&e]).map_err(RsomicsError::Io)?;
                }
                writeln!(output).map_err(RsomicsError::Io)?;
            }
            Ok(n_tested_exons as u64)
        }
        GeneTest::Gene => {
            let tested: Vec<&String> = gene_order
                .iter()
                .filter(|g| gene_lr.contains_key(*g))
                .collect();
            let pvals: Vec<f64> = tested
                .iter()
                .map(|g| {
                    let (lr, nexons) = gene_lr[*g];
                    special::pchisq_upper(lr, (nexons - 1) as f64)
                })
                .collect();
            let fdr = bh_fdr(&pvals);
            let mut header = "GeneID\tNExons\tgene.LR\tP.Value".to_string();
            if args.fdr {
                header.push_str("\tFDR");
            }
            writeln!(output, "{header}").map_err(RsomicsError::Io)?;
            for (i, g) in tested.iter().enumerate() {
                let (lr, nexons) = gene_lr[*g];
                write!(output, "{}\t{}\t{:.7e}\t{:.7e}", g, nexons, lr, pvals[i])
                    .map_err(RsomicsError::Io)?;
                if args.fdr {
                    write!(output, "\t{:.7e}", fdr[i]).map_err(RsomicsError::Io)?;
                }
                writeln!(output).map_err(RsomicsError::Io)?;
            }
            Ok(tested.len() as u64)
        }
        GeneTest::Simes => {
            let tested: Vec<&String> = gene_order
                .iter()
                .filter(|g| gene_simes.contains_key(*g))
                .collect();
            let pvals: Vec<f64> = tested.iter().map(|g| gene_simes[*g]).collect();
            let fdr = bh_fdr(&pvals);
            let mut header = "GeneID\tNExons\tP.Value".to_string();
            if args.fdr {
                header.push_str("\tFDR");
            }
            writeln!(output, "{header}").map_err(RsomicsError::Io)?;
            for (i, g) in tested.iter().enumerate() {
                let (_, nexons) = gene_lr[*g];
                write!(output, "{}\t{}\t{:.7e}", g, nexons, pvals[i]).map_err(RsomicsError::Io)?;
                if args.fdr {
                    write!(output, "\t{:.7e}", fdr[i]).map_err(RsomicsError::Io)?;
                }
                writeln!(output).map_err(RsomicsError::Io)?;
            }
            Ok(tested.len() as u64)
        }
    }
}

/// Orthonormal p×p basis (row-major, columns are basis vectors) whose first
/// column is `v` normalized; the rest complete it via Gram-Schmidt.
fn householder_basis(v: &[f64]) -> Vec<f64> {
    let p = v.len();
    let mut cols: Vec<Vec<f64>> = Vec::with_capacity(p);
    let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    cols.push(v.iter().map(|x| x / norm).collect());
    for e in 0..p {
        let mut cand = vec![0.0f64; p];
        cand[e] = 1.0;
        for c in &cols {
            let d: f64 = cand.iter().zip(c).map(|(a, b)| a * b).sum();
            for i in 0..p {
                cand[i] -= d * c[i];
            }
        }
        let nrm = cand.iter().map(|x| x * x).sum::<f64>().sqrt();
        if nrm > 1e-9 {
            for x in &mut cand {
                *x /= nrm;
            }
            cols.push(cand);
            if cols.len() == p {
                break;
            }
        }
    }
    let mut q = vec![0.0f64; p * p];
    for (j, c) in cols.iter().enumerate() {
        for (i, &val) in c.iter().enumerate() {
            q[i * p + j] = val;
        }
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deviance_zero_at_mle() {
        let y = [10.0, 12.0, 8.0];
        assert!(nb_deviance(&y, &y, 0.1).abs() < 1e-12);
    }

    #[test]
    fn simes_single_small() {
        // edgeR variant: sorted [0.001, 0.4, 0.5], n=3; min of 0.001·2, 0.4·1, 0.5·1.
        let p = [0.001, 0.5, 0.4];
        assert!((simes(&p) - 0.002).abs() < 1e-12);
    }

    #[test]
    fn solve_diag() {
        let mut a = vec![2.0, 0.0, 0.0, 3.0];
        let mut b = [4.0, 9.0];
        let mut x = vec![0.0; 2];
        assert!(solve(&mut a, &mut b, &mut x, 2));
        assert!((x[0] - 2.0).abs() < 1e-12 && (x[1] - 3.0).abs() < 1e-12);
    }
}
