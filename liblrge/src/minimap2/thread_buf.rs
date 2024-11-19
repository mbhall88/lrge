use minimap2_sys::{mm_tbuf_destroy, mm_tbuf_init, mm_tbuf_t};
use std::cell::RefCell;

// Thread local buffer (memory management) for minimap2
thread_local! {
    pub(crate) static BUF: RefCell<ThreadLocalBuffer> = RefCell::new(ThreadLocalBuffer::new());
}

/// ThreadLocalBuffer for minimap2 memory management
#[derive(Debug)]
pub(crate) struct ThreadLocalBuffer {
    buf: *mut mm_tbuf_t,
    max_uses: usize,
    uses: usize,
}

impl ThreadLocalBuffer {
    pub fn new() -> Self {
        let buf = unsafe { mm_tbuf_init() };
        Self {
            buf,
            max_uses: 15,
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
