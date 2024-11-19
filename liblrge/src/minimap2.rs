pub(crate) mod aligner;
pub(crate) mod mapping;
pub(crate) mod thread_buf;

use minimap2_sys::*;

pub(crate) use self::aligner::Aligner;

pub(crate) type MapOpt = mm_mapopt_t;
pub(crate) type IdxOpt = mm_idxopt_t;
pub(crate) const AVA_PB: &[u8] = b"ava-pb\0";
pub(crate) const AVA_ONT: &[u8] = b"ava-ont\0";
