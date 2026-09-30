//! Minimal little-endian binary serialisation helpers.

use std::io::{self, Read, Write};

pub struct Writer<'a> {
    inner: &'a mut dyn Write,
}

impl<'a> Writer<'a> {
    pub fn new(inner: &'a mut dyn Write) -> Self {
        Self { inner }
    }
    pub fn bytes(&mut self, b: &[u8]) -> io::Result<()> {
        self.inner.write_all(b)
    }
    pub fn u8(&mut self, x: u8) -> io::Result<()> {
        self.bytes(&[x])
    }
    pub fn bool(&mut self, x: bool) -> io::Result<()> {
        self.u8(x as u8)
    }
    pub fn u32(&mut self, x: u32) -> io::Result<()> {
        self.bytes(&x.to_le_bytes())
    }
    pub fn u64(&mut self, x: u64) -> io::Result<()> {
        self.bytes(&x.to_le_bytes())
    }
    pub fn u128(&mut self, x: u128) -> io::Result<()> {
        self.bytes(&x.to_le_bytes())
    }
    pub fn vec_u64(&mut self, v: &[u64]) -> io::Result<()> {
        self.u64(v.len() as u64)?;
        // write in chunks to avoid per-element syscalls on unbuffered writers
        let mut buf = Vec::with_capacity(8 * v.len().min(1 << 16));
        for chunk in v.chunks(1 << 16) {
            buf.clear();
            for x in chunk {
                buf.extend_from_slice(&x.to_le_bytes());
            }
            self.bytes(&buf)?;
        }
        Ok(())
    }
}

pub struct Reader<'a> {
    inner: &'a mut dyn Read,
}

impl<'a> Reader<'a> {
    pub fn new(inner: &'a mut dyn Read) -> Self {
        Self { inner }
    }
    pub fn bytes(&mut self, b: &mut [u8]) -> io::Result<()> {
        self.inner.read_exact(b)
    }
    pub fn u8(&mut self) -> io::Result<u8> {
        let mut b = [0u8; 1];
        self.bytes(&mut b)?;
        Ok(b[0])
    }
    pub fn bool(&mut self) -> io::Result<bool> {
        match self.u8()? {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(invalid("invalid bool")),
        }
    }
    pub fn u32(&mut self) -> io::Result<u32> {
        let mut b = [0u8; 4];
        self.bytes(&mut b)?;
        Ok(u32::from_le_bytes(b))
    }
    pub fn u64(&mut self) -> io::Result<u64> {
        let mut b = [0u8; 8];
        self.bytes(&mut b)?;
        Ok(u64::from_le_bytes(b))
    }
    pub fn u128(&mut self) -> io::Result<u128> {
        let mut b = [0u8; 16];
        self.bytes(&mut b)?;
        Ok(u128::from_le_bytes(b))
    }
    pub fn vec_u64(&mut self) -> io::Result<Vec<u64>> {
        let n = self.u64()? as usize;
        let mut v = Vec::with_capacity(n.min(1 << 24));
        let mut buf = vec![0u8; 8 * n.min(1 << 16)];
        let mut left = n;
        while left > 0 {
            let c = left.min(1 << 16);
            self.bytes(&mut buf[..8 * c])?;
            v.extend(buf[..8 * c].chunks_exact(8).map(|b| u64::from_le_bytes(b.try_into().unwrap())));
            left -= c;
        }
        Ok(v)
    }
}

pub fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}
