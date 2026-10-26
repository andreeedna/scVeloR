//! Streamline integration on a regular grid, a port of matplotlib's
//! `streamplot` algorithm (adaptive Heun / RK12 integration, a density mask that
//! prevents streamlines from crowding, seeds spiralling in from the border).
//! Rendering is left to the host language (ggplot2, matplotlib, ...).

#[derive(Clone, Debug)]
pub struct StreamOptions {
    /// Density of streamlines (matplotlib/scVelo: 2 -> 60 x 60 mask).
    pub density: (f64, f64),
    /// Minimum / maximum streamline length in axes coordinates.
    pub min_length: f64,
    pub max_length: f64,
    /// "both", "forward" or "backward".
    pub direction: Direction,
    pub broken_streamlines: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Direction {
    Both,
    Forward,
    Backward,
}

impl Default for StreamOptions {
    fn default() -> Self {
        StreamOptions { density: (2.0, 2.0), min_length: 0.1, max_length: 4.0, direction: Direction::Both, broken_streamlines: true }
    }
}

/// One streamline in data coordinates, plus the arrow segment matplotlib
/// would draw at its midpoint and the interpolated speed along it.
#[derive(Clone, Debug)]
pub struct Streamline {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub speed: Vec<f64>,
    pub arrow_tail: (f64, f64),
    pub arrow_head: (f64, f64),
}

struct Grid {
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    x0: f64,
    y0: f64,
}

impl Grid {
    fn within(&self, xi: f64, yi: f64) -> bool {
        xi >= 0.0 && xi <= (self.nx - 1) as f64 && yi >= 0.0 && yi <= (self.ny - 1) as f64
    }
}

struct Mask {
    nx: usize,
    ny: usize,
    m: Vec<bool>,
    traj: Vec<(usize, usize)>,
    current: Option<(usize, usize)>,
}

impl Mask {
    fn get(&self, xm: usize, ym: usize) -> bool {
        self.m[ym * self.nx + xm]
    }
    /// Returns false if the cell is occupied (InvalidIndexError).
    fn update(&mut self, xm: usize, ym: usize, broken: bool) -> bool {
        if self.current != Some((xm, ym)) {
            if !self.get(xm, ym) {
                self.traj.push((ym, xm));
                self.m[ym * self.nx + xm] = true;
                self.current = Some((xm, ym));
            } else if broken {
                return false;
            }
        }
        true
    }
    fn undo(&mut self) {
        for &(ym, xm) in &self.traj {
            self.m[ym * self.nx + xm] = false;
        }
    }
}

struct DomainMap {
    x_grid2mask: f64,
    y_grid2mask: f64,
}

impl DomainMap {
    fn grid2mask(&self, xi: f64, yi: f64) -> (usize, usize) {
        (round_half_even(xi * self.x_grid2mask) as usize, round_half_even(yi * self.y_grid2mask) as usize)
    }
}

/// Python's round() (banker's rounding).
fn round_half_even(v: f64) -> f64 {
    let r = v.round();
    if (v - v.trunc()).abs() == 0.5 {
        2.0 * (v / 2.0).round()
    } else {
        r
    }
}

/// Bilinear interpolation on an integer grid; NaN if any corner is NaN
/// (matplotlib masked arrays terminate the trajectory there).
fn interpgrid(a: &[f64], nx: usize, ny: usize, xi: f64, yi: f64) -> f64 {
    let x = xi as usize;
    let y = yi as usize;
    let xn = if x == nx - 1 { x } else { x + 1 };
    let yn = if y == ny - 1 { y } else { y + 1 };
    let a00 = a[y * nx + x];
    let a01 = a[y * nx + xn];
    let a10 = a[yn * nx + x];
    let a11 = a[yn * nx + xn];
    let xt = xi - x as f64;
    let yt = yi - y as f64;
    let a0 = a00 * (1.0 - xt) + a01 * xt;
    let a1 = a10 * (1.0 - xt) + a11 * xt;
    a0 * (1.0 - yt) + a1 * yt
}

enum Step {
    Ok(f64, f64),
    OutOfBounds,
    Terminate,
}

struct Field<'a> {
    u: Vec<f64>,
    v: Vec<f64>,
    speed: Vec<f64>,
    nx: usize,
    ny: usize,
    grid: &'a Grid,
}

impl<'a> Field<'a> {
    fn forward(&self, xi: f64, yi: f64) -> Step {
        if !self.grid.within(xi, yi) {
            return Step::OutOfBounds;
        }
        let ds_dt = interpgrid(&self.speed, self.nx, self.ny, xi, yi);
        if ds_dt == 0.0 || ds_dt.is_nan() {
            return Step::Terminate;
        }
        let dt_ds = 1.0 / ds_dt;
        let ui = interpgrid(&self.u, self.nx, self.ny, xi, yi);
        let vi = interpgrid(&self.v, self.nx, self.ny, xi, yi);
        if ui.is_nan() || vi.is_nan() {
            return Step::Terminate;
        }
        Step::Ok(ui * dt_ds, vi * dt_ds)
    }
    fn eval(&self, xi: f64, yi: f64, backward: bool) -> Step {
        match self.forward(xi, yi) {
            Step::Ok(a, b) if backward => Step::Ok(-a, -b),
            s => s,
        }
    }
}

fn integrate_rk12(
    x0: f64,
    y0: f64,
    field: &Field,
    dmap: &DomainMap,
    mask: &mut Mask,
    backward: bool,
    maxlength: f64,
    broken: bool,
) -> (f64, Vec<(f64, f64)>) {
    let maxerror = 0.003;
    let maxds = (1.0 / mask.nx as f64).min(1.0 / mask.ny as f64).min(0.1);
    let mut ds = maxds;
    let mut stotal = 0.0;
    let (mut xi, mut yi) = (x0, y0);
    let mut traj: Vec<(f64, f64)> = Vec::new();
    let (nx, ny) = (field.nx as f64, field.ny as f64);
    loop {
        if field.grid.within(xi, yi) {
            traj.push((xi, yi));
        } else {
            // OutOfBounds raised before evaluating
            if !traj.is_empty() {
                stotal += euler_to_boundary(&mut traj, field, backward);
            }
            break;
        }
        let k1 = field.eval(xi, yi, backward);
        let (k1x, k1y) = match k1 {
            Step::Ok(a, b) => (a, b),
            Step::OutOfBounds => {
                if !traj.is_empty() {
                    stotal += euler_to_boundary(&mut traj, field, backward);
                }
                break;
            }
            Step::Terminate => break,
        };
        let (k2x, k2y) = match field.eval(xi + ds * k1x, yi + ds * k1y, backward) {
            Step::Ok(a, b) => (a, b),
            Step::OutOfBounds => {
                if !traj.is_empty() {
                    stotal += euler_to_boundary(&mut traj, field, backward);
                }
                break;
            }
            Step::Terminate => break,
        };
        let dx1 = ds * k1x;
        let dy1 = ds * k1y;
        let dx2 = ds * 0.5 * (k1x + k2x);
        let dy2 = ds * 0.5 * (k1y + k2y);
        let error = (((dx2 - dx1) / (nx - 1.0)).powi(2) + ((dy2 - dy1) / (ny - 1.0)).powi(2)).sqrt();
        if error < maxerror {
            xi += dx2;
            yi += dy2;
            if !field.grid.within(xi, yi) {
                // update_trajectory raises InvalidIndexError outside the grid
                break;
            }
            let (xm, ym) = dmap.grid2mask(xi, yi);
            if !mask.update(xm, ym, broken) {
                break;
            }
            if stotal + ds > maxlength {
                break;
            }
            stotal += ds;
        }
        if error == 0.0 {
            ds = maxds;
        } else {
            ds = maxds.min(0.85 * ds * (maxerror / error).sqrt());
        }
    }
    (stotal, traj)
}

fn euler_to_boundary(traj: &mut Vec<(f64, f64)>, field: &Field, backward: bool) -> f64 {
    let (nx, ny) = (field.nx as f64, field.ny as f64);
    let (xi, yi) = *traj.last().unwrap();
    let (cx, cy) = match field.eval(xi, yi, backward) {
        Step::Ok(a, b) => (a, b),
        _ => return 0.0,
    };
    let dsx = if cx == 0.0 { f64::INFINITY } else if cx < 0.0 { xi / -cx } else { (nx - 1.0 - xi) / cx };
    let dsy = if cy == 0.0 { f64::INFINITY } else if cy < 0.0 { yi / -cy } else { (ny - 1.0 - yi) / cy };
    let ds = dsx.min(dsy);
    traj.push((xi + cx * ds, yi + cy * ds));
    ds
}

/// Spiral of mask cells starting at a corner (matplotlib _gen_starting_points).
fn starting_points(nx: usize, ny: usize) -> Vec<(usize, usize)> {
    let (mut xfirst, mut yfirst, mut xlast, mut ylast) = (0i64, 1i64, nx as i64 - 1, ny as i64 - 1);
    let (mut x, mut y) = (0i64, 0i64);
    let mut dir = 0; // right, up, left, down
    let mut out = Vec::with_capacity(nx * ny);
    for _ in 0..(nx * ny) {
        out.push((x as usize, y as usize));
        match dir {
            0 => {
                x += 1;
                if x >= xlast {
                    xlast -= 1;
                    dir = 1;
                }
            }
            1 => {
                y += 1;
                if y >= ylast {
                    ylast -= 1;
                    dir = 2;
                }
            }
            2 => {
                x -= 1;
                if x <= xfirst {
                    xfirst += 1;
                    dir = 3;
                }
            }
            _ => {
                y -= 1;
                if y <= yfirst {
                    yfirst += 1;
                    dir = 0;
                }
            }
        }
    }
    out
}

/// Compute streamlines for the velocity field `(u, v)` given on the regular
/// grid `xs` x `ys` (`u[y * nx + x]`, NaN = masked).
pub fn streamlines(xs: &[f64], ys: &[f64], u: &[f64], v: &[f64], opts: &StreamOptions) -> Vec<Streamline> {
    let (nx, ny) = (xs.len(), ys.len());
    assert!(nx >= 2 && ny >= 2);
    assert_eq!(u.len(), nx * ny);
    let grid = Grid {
        nx,
        ny,
        dx: xs[1] - xs[0],
        dy: ys[1] - ys[0],
        x0: xs[0],
        y0: ys[0],
    };
    let mnx = (30.0 * opts.density.0) as usize;
    let mny = (30.0 * opts.density.1) as usize;
    let mut mask = Mask { nx: mnx, ny: mny, m: vec![false; mnx * mny], traj: Vec::new(), current: None };
    let dmap = DomainMap {
        x_grid2mask: (mnx as f64 - 1.0) / (nx as f64 - 1.0),
        y_grid2mask: (mny as f64 - 1.0) / (ny as f64 - 1.0),
    };
    // data -> grid units, then axes units for the speed
    let ug: Vec<f64> = u.iter().map(|a| a / grid.dx).collect();
    let vg: Vec<f64> = v.iter().map(|a| a / grid.dy).collect();
    // matplotlib masks u and v independently; a masked component masks speed
    let speed: Vec<f64> = ug
        .iter()
        .zip(&vg)
        .map(|(a, b)| ((a / (nx as f64 - 1.0)).powi(2) + (b / (ny as f64 - 1.0)).powi(2)).sqrt())
        .collect();
    let field = Field { u: ug, v: vg, speed, nx, ny, grid: &grid };
    let maxlength = if opts.direction == Direction::Both { opts.max_length / 2.0 } else { opts.max_length };

    let mut trajectories: Vec<Vec<(f64, f64)>> = Vec::new();
    for (xm, ym) in starting_points(mnx, mny) {
        if mask.get(xm, ym) {
            continue;
        }
        let xg = xm as f64 / dmap.x_grid2mask;
        let yg = ym as f64 / dmap.y_grid2mask;
        // start_trajectory
        mask.traj.clear();
        let (sxm, sym) = dmap.grid2mask(xg, yg);
        if !mask.update(sxm, sym, opts.broken_streamlines) {
            continue;
        }
        let mut stotal = 0.0;
        let mut traj: Vec<(f64, f64)> = Vec::new();
        if opts.direction != Direction::Forward {
            let (s, t) = integrate_rk12(xg, yg, &field, &dmap, &mut mask, true, maxlength, opts.broken_streamlines);
            stotal += s;
            traj.extend(t.into_iter().rev());
        }
        if opts.direction != Direction::Backward {
            mask.current = Some(dmap.grid2mask(xg, yg));
            let (s, t) = integrate_rk12(xg, yg, &field, &dmap, &mut mask, false, maxlength, opts.broken_streamlines);
            stotal += s;
            traj.extend(t.into_iter().skip(1));
        }
        if stotal > opts.min_length {
            trajectories.push(traj);
        } else {
            mask.undo();
        }
    }

    trajectories
        .into_iter()
        .filter(|t| t.len() >= 2)
        .map(|t| {
            let x: Vec<f64> = t.iter().map(|p| p.0 * grid.dx + grid.x0).collect();
            let y: Vec<f64> = t.iter().map(|p| p.1 * grid.dy + grid.y0).collect();
            let speed: Vec<f64> = t.iter().map(|p| interp_raw(u, v, nx, ny, p.0, p.1)).collect();
            let mut s = Vec::with_capacity(x.len() - 1);
            let mut acc = 0.0;
            for i in 1..x.len() {
                acc += ((x[i] - x[i - 1]).powi(2) + (y[i] - y[i - 1]).powi(2)).sqrt();
                s.push(acc);
            }
            let half = acc / 2.0;
            let nmid = s.partition_point(|&v| v < half);
            let tail = (x[nmid], y[nmid]);
            let n2 = (nmid + 2).min(x.len());
            let head = (
                x[nmid..n2].iter().sum::<f64>() / (n2 - nmid) as f64,
                y[nmid..n2].iter().sum::<f64>() / (n2 - nmid) as f64,
            );
            Streamline { x, y, speed, arrow_tail: tail, arrow_head: head }
        })
        .collect()
}

fn interp_raw(u: &[f64], v: &[f64], nx: usize, ny: usize, xi: f64, yi: f64) -> f64 {
    let xi = xi.clamp(0.0, (nx - 1) as f64);
    let yi = yi.clamp(0.0, (ny - 1) as f64);
    let a = interpgrid(u, nx, ny, xi, yi);
    let b = interpgrid(v, nx, ny, xi, yi);
    (a * a + b * b).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_flow_gives_horizontal_lines() {
        let xs: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let ys = xs.clone();
        let u = vec![1.0; 400];
        let v = vec![0.0; 400];
        let lines = streamlines(&xs, &ys, &u, &v, &StreamOptions::default());
        assert!(!lines.is_empty());
        for l in &lines {
            let y0 = l.y[0];
            assert!(l.y.iter().all(|&y| (y - y0).abs() < 1e-9));
            assert!(l.arrow_head.0 >= l.arrow_tail.0);
        }
    }
}
