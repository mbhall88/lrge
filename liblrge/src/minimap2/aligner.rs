use super::mapping::Mapping;
use super::thread_buf::BUF;
use super::{IdxOpt, MapOpt};
use minimap2_sys::*;
use std::ffi::c_void;
use std::mem::MaybeUninit;
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

#[derive(Clone)]
pub(crate) struct Aligner {
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
        if self.idx.is_none() {
            return Err("No index");
        }

        // Make sure sequence is not empty
        if seq.is_empty() {
            return Err("Sequence is empty");
        }

        let mut mm_reg: MaybeUninit<*mut mm_reg1_t> = MaybeUninit::uninit();

        // Number of results
        let mut n_regs: i32 = 0;

        let qname = match query_name {
            None => std::ptr::null(),
            Some(qname) => qname.as_ptr() as *const ::std::os::raw::c_char,
        };
        let query_name = query_name.map(|q| q.to_vec()).unwrap_or(b"*".to_vec());
        let query_len = seq.len() as i32;

        let mappings = BUF.with(|buf| {
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
                    let reg: mm_reg1_t = *reg_ptr;

                    let contig: *mut ::std::os::raw::c_char =
                        (*((*(self.idx.unwrap())).seq.offset(reg.rid as isize))).name;
                    let target_name = std::ffi::CStr::from_ptr(contig).to_bytes().to_vec();

                    let strand = if reg.rev() == 0 { '+' } else { '-' };

                    // tp:A:<CHAR> Type of aln: P/primary, S/secondary and I,i/inversion
                    let tp = match (reg.id == reg.parent, reg.inv() != 0) {
                        (true, true) => 'I',
                        (true, false) => 'P',
                        (false, true) => 'i',
                        (false, false) => 'S',
                    };
                    // cm:i:<INT> which is the number of minimizers on the chain
                    let cm = reg.cnt;
                    // s1:i:<INT> which is the number of residues in the matching chain (chaining score)
                    let s1 = reg.score;
                    // dv:f:<FLOAT> approximate per-base sequence divergence
                    let dv = reg.div;
                    // rl:i:<INT> Length of query regions harboring repetitive seeds
                    let rl = (*buf.borrow_mut().get_buf()).rep_len;

                    mappings.push(Mapping {
                        target_name,
                        target_len: (*((*(self.idx.unwrap())).seq.offset(reg.rid as isize))).len
                            as i32,
                        target_start: reg.rs,
                        target_end: reg.re,
                        query_name: query_name.clone(),
                        query_len,
                        query_start: reg.qs,
                        query_end: reg.qe,
                        strand,
                        match_len: reg.mlen,
                        block_len: reg.blen,
                        mapq: reg.mapq(),
                        tp,
                        cm,
                        s1,
                        dv,
                        rl,
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