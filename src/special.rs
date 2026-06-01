//! chi-square upper-tail probability matching R's pchisq(lower.tail=FALSE),
//! built on a port of the public-domain pgamma (series + Lentz continued
//! fraction, DLMF 8.11). Accurate to ~1e-12, the precision the LRT p-values need.

fn lgamma(x: f64) -> f64 {
    libm::lgamma(x)
}
fn log1p(x: f64) -> f64 {
    libm::log1p(x)
}
fn expm1(x: f64) -> f64 {
    libm::expm1(x)
}

const M_LN2: f64 = std::f64::consts::LN_2;

fn log1mexp(x: f64) -> f64 {
    if x > -M_LN2 {
        (-expm1(x)).ln()
    } else {
        log1p(-x.exp())
    }
}

fn log_gammp_series(a: f64, x: f64) -> f64 {
    let mut ap = a;
    let mut del = 1.0 / a;
    let mut sum = del;
    loop {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * 1e-16 {
            break;
        }
    }
    -x + a * x.ln() - lgamma(a) + sum.ln()
}

fn log_gammq_cf(a: f64, x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / TINY;
    let mut d = 1.0 / b;
    let mut h = d;
    let mut i = 1.0;
    loop {
        let an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < TINY {
            d = TINY;
        }
        c = b + an / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-16 {
            break;
        }
        i += 1.0;
    }
    -x + a * x.ln() - lgamma(a) + h.ln()
}

fn pgamma_upper_log(x: f64, a: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        log1mexp(log_gammp_series(a, x))
    } else {
        log_gammq_cf(a, x)
    }
}

/// Upper-tail chi-square probability P(X > stat) on `df` degrees of freedom,
/// i.e. R's pchisq(stat, df, lower.tail=FALSE). The LRT p-value.
pub fn pchisq_upper(stat: f64, df: f64) -> f64 {
    if stat <= 0.0 {
        return 1.0;
    }
    pgamma_upper_log(stat / 2.0, df / 2.0).exp().min(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol * (1.0 + b.abs()), "{a} vs {b}");
    }

    #[test]
    fn pchisq_matches_r() {
        // R: pchisq(c(3.84, 10, 0.5), df=c(1,2,3), lower.tail=FALSE)
        close(pchisq_upper(3.84, 1.0), 0.0500435212, 1e-9);
        close(pchisq_upper(10.0, 2.0), 6.737946999e-3, 1e-9);
        close(pchisq_upper(0.5, 3.0), 0.9188914117, 1e-9);
    }
}
