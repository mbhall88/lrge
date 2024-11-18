use minimap2_sys::*;
use std::cell::RefCell;
use std::ffi::c_void;
use std::mem::MaybeUninit;
use std::num::NonZeroI32;
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

pub(crate) type MapOpt = mm_mapopt_t;
pub(crate) type IdxOpt = mm_idxopt_t;
pub(crate) const AVA_PB: &[u8] = b"ava-pb\0";
pub(crate) const AVA_ONT: &[u8] = b"ava-ont\0";

/// Strand enum
#[derive(Debug, PartialEq, Eq, Copy, Clone, Default)]
pub enum Strand {
    #[default]
    Forward,
    Reverse,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

/// Mapping result
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Mapping {
    // The query sequence name.
    pub query_name: Option<String>,
    pub query_len: Option<NonZeroI32>,
    pub query_start: i32,
    pub query_end: i32,
    pub strand: Strand,
    pub target_name: Option<String>,
    pub target_len: i32,
    pub target_start: i32,
    pub target_end: i32,
    pub match_len: i32,
    pub block_len: i32,
    pub mapq: u32,
    pub is_primary: bool,
    pub is_supplementary: bool,
}

// Thread local buffer (memory management) for minimap2
thread_local! {
    static BUF: RefCell<ThreadLocalBuffer> = RefCell::new(ThreadLocalBuffer::new());
}

/// ThreadLocalBuffer for minimap2 memory management
#[derive(Debug)]
struct ThreadLocalBuffer {
    buf: *mut mm_tbuf_t,
    max_uses: usize,
    uses: usize,
}

impl ThreadLocalBuffer {
    pub fn new() -> Self {
        let buf = unsafe { mm_tbuf_init() };
        Self {
            buf,
            max_uses: 1,
            uses: 0,
        }
    }
    /// Return the buffer, checking how many times it has been borrowed.
    /// Free the memory of the old buffer and reinitialise a new one If
    /// num_uses exceeds max_uses.
    pub fn get_buf(&mut self) -> *mut mm_tbuf_t {
        if self.uses > self.max_uses {
            self.free_buffer();
            let buf = unsafe { mm_tbuf_init() };
            self.buf = buf;
            self.uses = 0;
        }
        self.uses += 1;
        self.buf
    }

    fn free_buffer(&mut self) {
        unsafe { mm_tbuf_destroy(self.buf) };
    }
}

/// Handle destruction of thread local buffer properly.
impl Drop for ThreadLocalBuffer {
    fn drop(&mut self) {
        unsafe { mm_tbuf_destroy(self.buf) };
    }
}

impl Default for ThreadLocalBuffer {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone)]
pub struct Aligner {
    /// Index options passed to minimap2 (mm_idxopt_t)
    pub idxopt: IdxOpt,

    /// Mapping options passed to minimap2 (mm_mapopt_t)
    pub mapopt: MapOpt,

    /// Number of threads to create the index with
    pub threads: usize,

    /// Index created by minimap2
    pub idx: Option<*mut mm_idx_t>,

    /// Index reader created by minimap2
    pub idx_reader: Option<mm_idx_reader_t>,
}

/// Create a default aligner
impl Default for Aligner {
    fn default() -> Self {
        Self {
            idxopt: Default::default(),
            mapopt: Default::default(),
            threads: 1,
            idx: None,
            idx_reader: None,
        }
    }
}

mod send {
    use super::*;
    unsafe impl Sync for Aligner {}
    unsafe impl Send for Aligner {}
}

impl Aligner {
    /// Create a new aligner with default options
    pub fn builder() -> Self {
        // let mut aligner = Aligner {
        //     mapopt: MapOpt {
        //         seed: 11,
        //         best_n: 5,
        //         ..Default::default()
        //     },
        //     ..Default::default()
        // };
        let mut aligner = Aligner::default();

        let result = unsafe {
            let preset: i32 = 0;
            mm_set_opt(
                preset as *const i8,
                &mut aligner.idxopt,
                &mut aligner.mapopt,
            )
        };
        eprintln!("{:?}", aligner.mapopt);
        // if the result is 0, success, if -1 then issue with preset
        if result == -1 {
            panic!("Preset not found: {:?}", aligner.mapopt);
        }

        aligner
    }

    pub fn preset(mut self, preset: &[u8]) -> Self {
        let result = unsafe {
            mm_set_opt(
                preset.as_ptr() as *const i8,
                &mut self.idxopt,
                &mut self.mapopt,
            )
        };

        eprintln!("{:?}", self.mapopt);
        // if the result is 0, success, if -1 then the preset is not found
        if result == -1 {
            panic!("Preset not found: {:?}", String::from_utf8_lossy(preset));
        }

        self
    }

    pub fn dual(mut self, yes: bool) -> Self {
        if yes {
            // Set the `--dual=yes` flag. to do this, we need to clear the bit corresponding to
            // MM_F_NO_DUAL (0x002) in the `mapopt` field of the aligner.
            // this executes following https://github.com/lh3/minimap2/blob/618d33515e5853c4576d5a3d126fdcda28f0e8a4/main.c#L120
            self.mapopt.flag &= !0x002;
        }
        self
    }

    /// Sets the number of threads minimap2 will use for building the index
    pub fn with_index_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    pub fn check_opts(&self) -> Result<(), &'static str> {
        let result = unsafe { mm_check_opt(&self.idxopt, &self.mapopt) };

        if result == 0 {
            Ok(())
        } else {
            Err("Invalid options")
        }
    }

    /// Set index parameters for minimap2 using builder pattern
    /// Creates the index as well with the given number of threads (set at struct creation).
    /// You must set the number of threads before calling this function.
    ///
    /// Parameters:
    /// path: Location of pre-built index or FASTA/FASTQ file (may be gzipped or plaintext)
    /// Output: Option (None) or a filename
    ///
    /// Returns the aligner with the index set
    pub fn with_index<P>(mut self, path: P, output: Option<&str>) -> Result<Self, &'static str>
    where
        P: AsRef<Path>,
    {
        match self.set_index(path, output) {
            Ok(_) => Ok(self),
            Err(e) => Err(e),
        }
    }

    /// Set the index (in-place, without builder pattern)
    pub fn set_index<P>(&mut self, path: P, output: Option<&str>) -> Result<(), &'static str>
    where
        P: AsRef<Path>,
    {
        let path_str = match std::ffi::CString::new(path.as_ref().as_os_str().as_bytes()) {
            Ok(path) => path,
            Err(_) => {
                return Err("Invalid Path");
            }
        };

        // Confirm file exists
        if !path.as_ref().exists() {
            return Err("File does not exist");
        }

        // Confirm file is not empty
        if path.as_ref().metadata().unwrap().len() == 0 {
            return Err("File is empty");
        }

        let output = match output {
            Some(output) => match std::ffi::CString::new(output) {
                Ok(output) => output,
                Err(_) => return Err("Invalid Output"),
            },
            None => std::ffi::CString::new(Vec::new()).unwrap(),
        };

        let idx_reader = MaybeUninit::new(unsafe {
            mm_idx_reader_open(path_str.as_ptr(), &self.idxopt, output.as_ptr())
        });

        let idx;

        let idx_reader = unsafe { idx_reader.assume_init() };

        unsafe {
            // Just a test read? Just following: https://github.com/lh3/minimap2/blob/master/python/mappy.pyx#L147
            idx = MaybeUninit::new(mm_idx_reader_read(
                // self.idx_reader.as_mut().unwrap() as *mut mm_idx_reader_t,
                &mut *idx_reader as *mut mm_idx_reader_t,
                self.threads as i32,
            ));
            // Close the reader
            mm_idx_reader_close(idx_reader);
            // Set index opts
            mm_mapopt_update(&mut self.mapopt, *idx.as_ptr());
            // Idx index name
            mm_idx_index_name(idx.assume_init());
        }

        self.idx = Some(unsafe { idx.assume_init() });

        Ok(())
    }

    /// Aligns a given sequence (as bytes) to the index associated with this aligner
    ///
    /// Parameters:
    /// seq: Sequence to align
    /// cs: Whether to output CIGAR string
    /// MD: Whether to output MD tag
    /// max_frag_len: Maximum fragment length
    /// extra_flags: Extra flags to pass to minimap2 as `Vec<u64>`
    pub fn map(&self, seq: &[u8], query_name: Option<&[u8]>) -> Result<Vec<Mapping>, &'static str> {
        // Make sure index is set
        if !self.idx.is_some() {
            return Err("No index");
        }

        // Make sure sequence is not empty
        if seq.is_empty() {
            return Err("Sequence is empty");
        }

        let mut mm_reg: MaybeUninit<*mut mm_reg1_t> = MaybeUninit::uninit();

        // Number of results
        let mut n_regs: i32 = 0;
        // let mut map_opt = self.mapopt.clone();

        let qname = match query_name {
            None => std::ptr::null(),
            Some(qname) => qname.as_ptr() as *const ::std::os::raw::c_char,
        };

        let mappings = BUF.with(|buf| {
            let km = unsafe { mm_tbuf_get_km(buf.borrow_mut().get_buf()) };

            mm_reg = MaybeUninit::new(unsafe {
                mm_map(
                    self.idx.unwrap() as *const mm_idx_t,
                    seq.len() as i32,
                    seq.as_ptr() as *const ::std::os::raw::c_char,
                    &mut n_regs,
                    buf.borrow_mut().get_buf(),
                    &self.mapopt,
                    qname,
                )
            });
            let mut mappings = Vec::with_capacity(n_regs as usize);

            for i in 0..n_regs {
                unsafe {
                    let reg_ptr = (*mm_reg.as_ptr()).offset(i as isize);
                    let const_ptr = reg_ptr as *const mm_reg1_t;
                    let reg: mm_reg1_t = *reg_ptr;

                    let contig: *mut ::std::os::raw::c_char =
                        (*((*(self.idx.unwrap())).seq.offset(reg.rid as isize))).name;

                    let is_primary = reg.parent == reg.id;
                    let is_supplementary = reg.sam_pri() == 0;

                    mappings.push(Mapping {
                        target_name: Some(
                            std::ffi::CStr::from_ptr(contig)
                                .to_str()
                                .unwrap()
                                .to_string(),
                        ),
                        target_len: (*((*(self.idx.unwrap())).seq.offset(reg.rid as isize))).len
                            as i32,
                        target_start: reg.rs,
                        target_end: reg.re,
                        query_name: query_name.map(|q| String::from_utf8_lossy(q).to_string()),
                        query_len: NonZeroI32::new(seq.len() as i32),
                        query_start: reg.qs,
                        query_end: reg.qe,
                        strand: if reg.rev() == 0 {
                            Strand::Forward
                        } else {
                            Strand::Reverse
                        },
                        match_len: reg.mlen,
                        block_len: reg.blen,
                        mapq: reg.mapq(),
                        is_primary,
                        is_supplementary,
                    });
                    libc::free(reg.p as *mut c_void);
                }
            }
            mappings
        });
        // free some stuff here
        unsafe {
            let ptr: *mut mm_reg1_t = mm_reg.assume_init();
            let c_void_ptr: *mut c_void = ptr as *mut c_void;
            libc::free(c_void_ptr);
        }
        Ok(mappings)
    }
}
