use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Section};

use rsomics_edger_diff_splice::{DiffSpliceArgs, GeneTest, diff_splice};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Copy, Clone, Debug, ValueEnum)]
enum TestKind {
    Exon,
    #[value(alias = "Simes")]
    Simes,
    Gene,
}

#[derive(Parser, Debug)]
#[command(name = "rsomics-edger-diff-splice", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    pub counts: PathBuf,
    #[arg(long, value_name = "PATH")]
    design: PathBuf,
    #[arg(long, value_name = "PATH")]
    genes: PathBuf,
    #[arg(long, value_name = "N")]
    coef: Option<usize>,
    #[arg(long, value_name = "PATH")]
    contrast: Option<PathBuf>,
    #[arg(long, default_value_t = 0.05)]
    dispersion: f64,
    #[arg(long, value_name = "PATH")]
    dispersion_file: Option<PathBuf>,
    #[arg(long, value_name = "PATH")]
    norm_factors: Option<PathBuf>,
    #[arg(long, value_enum, default_value_t = TestKind::Exon)]
    test: TestKind,
    #[arg(long)]
    fdr: bool,
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.common.install_rayon_pool()?;
        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        let test = match self.test {
            TestKind::Exon => GeneTest::Exon,
            TestKind::Simes => GeneTest::Simes,
            TestKind::Gene => GeneTest::Gene,
        };
        let n = diff_splice(
            &DiffSpliceArgs {
                counts: &self.counts,
                design: &self.design,
                genes: &self.genes,
                coef: self.coef,
                contrast: self.contrast.as_deref(),
                dispersion: self.dispersion,
                dispersion_file: self.dispersion_file.as_deref(),
                norm_factors: self.norm_factors.as_deref(),
                test,
                fdr: self.fdr,
            },
            &mut out,
        )?;
        if !self.common.quiet {
            eprintln!("{n} rows reported");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Differential exon usage from an exon-level NB-GLM (edgeR diffSpliceDGE/topSpliceDGE).",
    origin: None,
    usage_lines: &[
        "<counts.tsv> --design design.tsv --genes genes.tsv [--coef N | --contrast c.tsv] [--test exon|Simes|gene] [--dispersion D | --dispersion-file f.tsv] [--norm-factors f.tsv] [--fdr] [-o out.tsv]",
    ],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "design",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Design matrix TSV: a header of coefficient names then one numeric row per sample.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "genes",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Gene id per exon, one per line in count-matrix row order (header optional).",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "coef",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("last coefficient"),
                description: "1-based design column tested for differential exon usage.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "contrast",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Contrast vector (one weight per design coefficient) to test instead of --coef.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "test",
                aliases: &[],
                value: Some("<kind>"),
                type_hint: Some("exon|Simes|gene"),
                required: false,
                default: Some("exon"),
                description: "Report per-exon stats, per-gene Simes p-values, or the gene-level LR test.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "dispersion",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0.05"),
                description: "Common negative-binomial dispersion shared across exons.",
                why_default: Some("edgeR's fallback when no dispersion is estimated."),
            },
            FlagSpec {
                short: None,
                long: "dispersion-file",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Per-exon dispersions (one per row, exon order), overriding --dispersion.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "norm-factors",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Per-sample normalization factors; multiplied into library sizes.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "fdr",
                aliases: &[],
                value: None,
                type_hint: Some("flag"),
                required: false,
                default: None,
                description: "Append a Benjamini-Hochberg FDR column.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Per-exon differential usage for the last coefficient",
            command: "rsomics-edger-diff-splice counts.tsv --design design.tsv --genes genes.tsv -o exon.tsv",
        },
        Example {
            description: "Gene-level Simes calls testing coefficient 2",
            command: "rsomics-edger-diff-splice counts.tsv --design design.tsv --genes genes.tsv --coef 2 --test Simes --fdr -o gene.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
