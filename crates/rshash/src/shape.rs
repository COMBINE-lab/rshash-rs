//! Gapped k-mer shapes (port of `shape.hpp`).
//!
//! A shape is an integer bit mask: bit `i` (from the least significant bit)
//! says whether base `i` of the window is a "care" position. Because windows
//! store their first base in the lowest bits, the binary literal reads
//! *reversed* with respect to sequence order (irrelevant for palindromic
//! shapes).
//!
//! Several shapes are aligned in one common window of `L = K + 2·O` bases
//! around a contiguous *kernel* of `K` bases (the shortest of the shapes'
//! longest runs of ones), which starts at base `O` of the window. Minimisers
//! are computed on the kernel only, so every shape's care positions include
//! the kernel and the kernel is centred (strand-symmetric).

use crate::word::{mask128, pdep64};
use std::fmt;

/// Longest run of ones (`find_long_run`): the first (lowest-bit) run wins ties.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Run {
    pub start: u32,
    pub len: u32,
}

impl Run {
    pub fn end(&self) -> u32 {
        self.start + self.len
    }
}

pub fn find_long_run(shape: u64) -> Run {
    let mut best = Run { start: 64, len: 0 };
    let mut pos = 0u32;
    let mut x = shape;
    while x != 0 {
        let zeros = x.trailing_zeros();
        pos += zeros;
        x >>= zeros;
        let len = (!x).trailing_zeros();
        if len > best.len {
            best = Run { start: pos, len };
        }
        pos += len;
        x = if len >= 64 { 0 } else { x >> len };
    }
    best
}

/// Bit-reverse the `bit_width(x)` significant bits of `x` (`reverse_shape`).
pub fn reverse_shape(x: u64) -> u64 {
    assert!(x != 0);
    let len = 64 - x.leading_zeros();
    x.reverse_bits() >> (64 - len)
}

/// Whether the shape is a palindrome (`canonical_shape` in C++).
pub fn is_palindrome(x: u64) -> bool {
    x == 0 || x == reverse_shape(x)
}

/// Base-level shape to 2-bit-per-base mask (`compute_shape_mask`).
pub fn shape_mask(shape: u64) -> u128 {
    let lo = pdep64(shape & 0xFFFF_FFFF, 0x5555_5555_5555_5555);
    let hi = pdep64(shape >> 32, 0x5555_5555_5555_5555);
    let lo = lo | (lo << 1);
    let hi = hi | (hi << 1);
    (lo as u128) | ((hi as u128) << 64)
}

/// Largest supported shape / window length in bases.
pub const MAX_WINDOW: u32 = 63;
/// Largest supported number of care positions per shape (keys must fit in 64 bits
/// for windows of up to 32 bases; kept for all windows as in C++).
pub const MAX_WEIGHT: u32 = 32;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Shape {
    pub value: u64,
    /// Span in bases (`bit_width(value)`).
    pub length: u32,
    /// Number of care positions.
    pub weight: u32,
    /// Longest run of ones: the shape's candidate kernel.
    pub run: Run,
    /// Bases before the run (C++: `overlap_right`).
    pub overlap_right: u32,
    /// Bases after the run (C++: `overlap_left`).
    pub overlap_left: u32,
    /// 2-bit-per-base mask over the shape's own window.
    pub mask: u128,
    pub is_palindrome: bool,
}

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum ShapeError {
    #[error("shape must be non-zero")]
    Zero,
    #[error("shape {0:#b}: length {1} exceeds the maximum of {MAX_WINDOW}")]
    TooLong(u64, u32),
    #[error("shape {0:#b}: weight {1} exceeds the maximum of {MAX_WEIGHT}")]
    TooHeavy(u64, u32),
    #[error("no shapes given")]
    Empty,
    #[error("aligned shapes window length {0} exceeds the maximum of {MAX_WINDOW}")]
    WindowTooLong(u32),
}

impl Shape {
    pub fn new(value: u64) -> Result<Self, ShapeError> {
        if value == 0 {
            return Err(ShapeError::Zero);
        }
        let length = 64 - value.leading_zeros();
        if length > MAX_WINDOW {
            return Err(ShapeError::TooLong(value, length));
        }
        let weight = value.count_ones();
        if weight > MAX_WEIGHT {
            return Err(ShapeError::TooHeavy(value, weight));
        }
        let run = find_long_run(value);
        Ok(Self {
            value,
            length,
            weight,
            run,
            overlap_right: run.start,
            overlap_left: length - run.end(),
            mask: shape_mask(value),
            is_palindrome: is_palindrome(value),
        })
    }

    /// Pretty string in sequence order (first base first), `1` = care, `0` = don't care.
    pub fn pattern(&self) -> String {
        (0..self.length).map(|i| if (self.value >> i) & 1 == 1 { '1' } else { '0' }).collect()
    }
}

/// A shape placed in the common window.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AlignedShape {
    pub shape: Shape,
    /// First base of the shape's span within the window.
    pub start: u32,
    /// One past the last base of the shape's span within the window.
    pub end: u32,
    /// 2-bit mask over the window (`w_mask`).
    pub w_mask: u128,
}

/// A set of shapes aligned around a common kernel (`Shapes32`/`Shapes64`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Shapes {
    pub shapes: Vec<AlignedShape>,
    /// Kernel length `K` (contiguous, used for minimisers).
    pub kernel_length: u32,
    /// Bases on each side of the kernel, `O`.
    pub overlap: u32,
    /// Window length `L = K + 2·O`.
    pub length: u32,
    /// 2-bit mask selecting the kernel within the window.
    pub kernel_mask: u128,
}

impl Shapes {
    pub fn new(values: &[u64]) -> Result<Self, ShapeError> {
        if values.is_empty() {
            return Err(ShapeError::Empty);
        }
        let shapes: Vec<Shape> = values.iter().map(|&v| Shape::new(v)).collect::<Result<_, _>>()?;
        let kernel_length = shapes.iter().map(|s| s.run.len).min().unwrap();
        // Every aligned shape must fit into the window: the run starts at O, so
        // O >= overlap_right and O + run + overlap_left <= K + 2·O.
        // (C++ uses max(overlap_left, overlap_right), which lets shapes whose
        // run is longer than K spill past the window; see DEVIATIONS.md.)
        let overlap =
            shapes.iter().map(|s| s.overlap_right.max(s.overlap_left + s.run.len - kernel_length)).max().unwrap();
        let length = kernel_length + 2 * overlap;
        if length > MAX_WINDOW {
            return Err(ShapeError::WindowTooLong(length));
        }
        let shapes = shapes
            .into_iter()
            .map(|shape| {
                let start = overlap - shape.overlap_right;
                AlignedShape { shape, start, end: start + shape.length, w_mask: shape.mask << (2 * start) }
            })
            .collect();
        Ok(Self { shapes, kernel_length, overlap, length, kernel_mask: mask128(2 * kernel_length) << (2 * overlap) })
    }

    pub fn values(&self) -> Vec<u64> {
        self.shapes.iter().map(|s| s.shape.value).collect()
    }
}

impl fmt::Display for Shapes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for s in &self.shapes {
            writeln!(
                f,
                "shape {} ({:#x}): weight {}, length {}, kernel run {}..{}, window span {}..{}{}",
                s.shape.pattern(),
                s.shape.value,
                s.shape.weight,
                s.shape.length,
                s.shape.run.start,
                s.shape.run.end(),
                s.start,
                s.end,
                if s.shape.is_palindrome { ", palindromic" } else { "" }
            )?;
        }
        write!(f, "window length {}, kernel length {}, overlap {}", self.length, self.kernel_length, self.overlap)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn benchmark_shapes() {
        // 0x5FFFFFFD = 1011111...1101 (palindrome): run bits 2..29 (27 ones)
        let s = Shape::new(0x5FFF_FFFD).unwrap();
        assert!(s.is_palindrome);
        assert_eq!(s.run, Run { start: 2, len: 27 });
        assert_eq!((s.overlap_right, s.overlap_left, s.length, s.weight), (2, 2, 31, 29));
        let w = Shapes::new(&[0x5FFF_FFFD]).unwrap();
        assert_eq!((w.kernel_length, w.overlap, w.length), (27, 2, 31));

        // 0x7FFFFFF9 (benchmark_rshash.sh): run bits 3..31, not a palindrome
        let s = Shape::new(0x7FFF_FFF9).unwrap();
        assert!(!s.is_palindrome);
        assert_eq!(s.run, Run { start: 3, len: 28 });
        assert_eq!((s.overlap_right, s.overlap_left), (3, 0));
        let w = Shapes::new(&[0x7FFF_FFF9]).unwrap();
        assert_eq!((w.kernel_length, w.overlap, w.length), (28, 3, 34));
        let a = w.shapes[0];
        assert_eq!((a.start, a.end), (0, 31));
    }

    #[test]
    fn multi_shape_alignment() {
        // runs of different lengths: the longer run must still fit
        let w = Shapes::new(&[0b1011111101, 0b111111111111]).unwrap();
        assert_eq!(w.kernel_length, 6);
        for s in &w.shapes {
            assert!(s.end <= w.length);
            // kernel inside care positions
            assert_eq!(s.w_mask & w.kernel_mask, w.kernel_mask);
        }
    }

    #[test]
    fn masks_and_reverse() {
        assert_eq!(shape_mask(0b101), 0b110011);
        assert_eq!(reverse_shape(0b1101), 0b1011);
        assert!(is_palindrome(0b1011101));
        let long = (1u64 << 62) | 0x3FFF_FFFF;
        let s = Shape::new(long).unwrap();
        assert_eq!(s.length, 63);
        assert_eq!(s.mask >> 124, 0b11);
        assert!(Shape::new(u64::MAX >> 1).is_err()); // weight 63
    }
}
