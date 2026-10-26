//! Small numeric helpers that reproduce numpy / scipy semantics exactly where
//! scVelo depends on them (percentiles with linear interpolation, ddof=0 std,
//! Scott's-rule Gaussian KDE, ...).

/// numpy.percentile with the default "linear" method. NaNs propagate (returns NaN).
pub fn percentile(values: &[f64], q: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    if values.iter().any(|v| v.is_nan()) {
        return f64::NAN;
    }
    let mut v = values.to_vec();
    v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    percentile_sorted(&v, q)
}

/// numpy.nanpercentile: NaNs are ignored.
pub fn nanpercentile(values: &[f64], q: f64) -> f64 {
    let mut v: Vec<f64> = values.iter().copied().filter(|x| !x.is_nan()).collect();
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    percentile_sorted(&v, q)
}

/// Linear-interpolated percentile on an already sorted slice.
pub fn percentile_sorted(sorted: &[f64], q: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    if n == 1 {
        return sorted[0];
    }
    let pos = q / 100.0 * (n as f64 - 1.0);
    let lo = pos.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    let frac = pos - lo as f64;
    let (a, b) = (sorted[lo], sorted[hi]);
    // numpy uses a lerp that is exact at the end points
    if frac == 0.0 {
        a
    } else {
        a + (b - a) * frac
    }
}

pub fn mean(v: &[f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.iter().sum::<f64>() / v.len() as f64
}

pub fn nanmean(v: &[f64]) -> f64 {
    let mut s = 0.0;
    let mut n = 0usize;
    for &x in v {
        if !x.is_nan() {
            s += x;
            n += 1;
        }
    }
    if n == 0 {
        f64::NAN
    } else {
        s / n as f64
    }
}

/// Population standard deviation (numpy default ddof=0).
pub fn std(v: &[f64]) -> f64 {
    var(v).sqrt()
}

pub fn var(v: &[f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    let m = mean(v);
    v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / v.len() as f64
}

pub fn max(v: &[f64]) -> f64 {
    v.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

pub fn min(v: &[f64]) -> f64 {
    v.iter().copied().fold(f64::INFINITY, f64::min)
}

/// Index of the minimum (first occurrence), like numpy.argmin (NaN wins like numpy).
pub fn argmin(v: &[f64]) -> usize {
    let mut best = 0;
    for (i, &x) in v.iter().enumerate() {
        if x.is_nan() {
            return i;
        }
        if x < v[best] {
            best = i;
        }
    }
    best
}

pub fn argmax(v: &[f64]) -> usize {
    let mut best = 0;
    for (i, &x) in v.iter().enumerate() {
        if x.is_nan() {
            return i;
        }
        if x > v[best] {
            best = i;
        }
    }
    best
}

/// numpy.linspace(a, b, num) (endpoint included).
pub fn linspace(a: f64, b: f64, num: usize) -> Vec<f64> {
    if num == 0 {
        return vec![];
    }
    if num == 1 {
        return vec![a];
    }
    let step = (b - a) / (num as f64 - 1.0);
    let mut out: Vec<f64> = (0..num).map(|i| a + step * i as f64).collect();
    out[num - 1] = b;
    out
}

/// Survival function of the standard normal distribution.
pub fn norm_sf(x: f64) -> f64 {
    0.5 * libm::erfc(x / std::f64::consts::SQRT_2)
}

/// scipy.stats.gaussian_kde(x)(grid) for 1-d data with Scott's bandwidth.
/// Returns None if the data variance is zero (scipy raises in that case).
pub fn gaussian_kde_1d(x: &[f64], grid: &[f64]) -> Option<Vec<f64>> {
    let n = x.len();
    if n < 2 {
        return None;
    }
    let m = mean(x);
    let var_ddof1 = x.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / (n as f64 - 1.0);
    if !(var_ddof1 > 0.0) || !var_ddof1.is_finite() {
        return None;
    }
    let factor = (n as f64).powf(-1.0 / 5.0);
    let cov = var_ddof1 * factor * factor;
    let norm = 1.0 / ((2.0 * std::f64::consts::PI * cov).sqrt() * n as f64);
    Some(
        grid.iter()
            .map(|&g| {
                let s: f64 = x.iter().map(|&xi| (-(g - xi) * (g - xi) / (2.0 * cov)).exp()).sum();
                s * norm
            })
            .collect(),
    )
}

/// Min-max scaling into [0, 1] ignoring non-finite values (scvelo.tools.utils.scale).
pub fn scale01(x: &[f64]) -> Vec<f64> {
    let finite: Vec<f64> = x.iter().copied().filter(|v| v.is_finite()).collect();
    if finite.is_empty() {
        return x.to_vec();
    }
    let lo = min(&finite);
    let shifted: Vec<f64> = x.iter().map(|v| v - lo).collect();
    let hi = shifted.iter().copied().filter(|v| v.is_finite()).fold(f64::NEG_INFINITY, f64::max);
    if hi != 0.0 {
        shifted.iter().map(|v| v / hi).collect()
    } else {
        shifted
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn percentile_matches_numpy() {
        let v = [1.0, 2.0, 3.0, 4.0];
        assert!((percentile(&v, 50.0) - 2.5).abs() < 1e-12);
        assert!((percentile(&v, 95.0) - 3.85).abs() < 1e-12);
        assert_eq!(percentile(&v, 0.0), 1.0);
        assert_eq!(percentile(&v, 100.0), 4.0);
    }

    #[test]
    fn sf() {
        assert!((norm_sf(0.0) - 0.5).abs() < 1e-15);
        assert!((norm_sf(3.0902323061678132) - 1e-3).abs() < 1e-9);
    }
}
