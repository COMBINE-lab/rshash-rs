//! Build parameters and the derived window geometry.

use crate::io::{Reader, Writer, invalid};
use crate::shape::{MAX_WINDOW, ShapeError, Shapes};
use std::io;

/// Construction parameters (the `rshash build` options).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Params {
    /// k-mer length (ignored when shapes are given: the kernel length is used).
    pub k: u32,
    /// Number of minimiser levels (1..=3).
    pub levels: u32,
    /// Minimiser length per level.
    pub m: [u32; 3],
    /// Maximum bucket size per level (`--t1..--t3`).
    pub t: [u16; 3],
    /// Last level: k-mers (or gapped k-mers) occurring `>= threshold` times are
    /// dropped; 0 disables the filter (`-t`).
    pub threshold: u16,
    /// Build locate support (`--loc`).
    pub loc: bool,
    /// Last level as hash tables (default) or Elias–Fano (`--ht` in C++).
    pub use_ht: bool,
    /// Gapped k-mer shapes; empty = contiguous k-mers.
    pub shapes: Vec<u64>,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            k: 31,
            levels: 2,
            m: [18, 21, 23],
            t: [64, 64, 64],
            threshold: 0,
            loc: false,
            use_ht: true,
            shapes: vec![],
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum ParamError {
    #[error(transparent)]
    Shape(#[from] ShapeError),
    #[error("number of levels must be 1, 2 or 3 (got {0})")]
    Levels(u32),
    #[error("k must be in 1..={MAX_WINDOW} (got {0})")]
    K(u32),
    #[error("minimiser length m{0}={1} must be in 1..=min(31, k={2})")]
    M(usize, u32, u32),
    #[error("at most {max} shapes are supported (got {0})", max = crate::query::MAX_SHAPES)]
    TooManyShapes(usize),
}

/// Window geometry derived from the parameters.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Geometry {
    /// Kernel length `K`: the contiguous k-mer used for minimisers.
    pub kernel: u32,
    /// Bases on each side of the kernel (`O`; 0 without shapes).
    pub overlap: u32,
    /// Window length `L = K + 2·O` (= k without shapes).
    pub window: u32,
    pub shapes: Option<Shapes>,
}

impl Params {
    pub fn geometry(&self) -> Result<Geometry, ParamError> {
        if !(1..=3).contains(&self.levels) {
            return Err(ParamError::Levels(self.levels));
        }
        let g = if self.shapes.is_empty() {
            if self.k == 0 || self.k > MAX_WINDOW {
                return Err(ParamError::K(self.k));
            }
            Geometry { kernel: self.k, overlap: 0, window: self.k, shapes: None }
        } else {
            if self.shapes.len() > crate::query::MAX_SHAPES {
                return Err(ParamError::TooManyShapes(self.shapes.len()));
            }
            let s = Shapes::new(&self.shapes)?;
            Geometry { kernel: s.kernel_length, overlap: s.overlap, window: s.length, shapes: Some(s) }
        };
        for l in 0..self.levels as usize {
            let m = self.m[l];
            if m == 0 || m > 31 || m > g.kernel {
                return Err(ParamError::M(l + 1, m, g.kernel));
            }
        }
        Ok(g)
    }

    pub fn write(&self, w: &mut Writer) -> io::Result<()> {
        w.u32(self.k)?;
        w.u32(self.levels)?;
        for l in 0..3 {
            w.u32(self.m[l])?;
            w.u32(self.t[l] as u32)?;
        }
        w.u32(self.threshold as u32)?;
        w.bool(self.loc)?;
        w.bool(self.use_ht)?;
        w.vec_u64(&self.shapes)
    }

    pub fn read(r: &mut Reader) -> io::Result<Self> {
        let k = r.u32()?;
        let levels = r.u32()?;
        let mut m = [0; 3];
        let mut t = [0; 3];
        for l in 0..3 {
            m[l] = r.u32()?;
            t[l] = u16::try_from(r.u32()?).map_err(|_| invalid("bad threshold"))?;
        }
        let threshold = u16::try_from(r.u32()?).map_err(|_| invalid("bad threshold"))?;
        let loc = r.bool()?;
        let use_ht = r.bool()?;
        let shapes = r.vec_u64()?;
        Ok(Self { k, levels, m, t, threshold, loc, use_ht, shapes })
    }
}
