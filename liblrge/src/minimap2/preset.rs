/// Preset options for minimap2
#[allow(dead_code)]
pub(crate) enum Preset {
    /// Align noisy long reads of ~10% error rate to a reference genome. This is the default mode.
    MapOnt,
    /// Align PacBio high-fidelity (HiFi) reads to a reference genome (-k19 -w19 -U50,500 -g10k -A1 -B4 -O6,26 -E2,1 -s200).
    MapHifi,
    /// Align older PacBio continuous long (CLR) reads to a reference genome (-Hk19).
    MapPb,
    /// accurate long reads (error rate <1%) against a reference genome
    LongReadHq,
    /// Long assembly to reference mapping (-k19 -w19 -U50,500 --rmq -r100k -g10k -A1 -B19 -O39,81 -E3,1 -s200 -z200 -N50). Typically, the alignment will not extend to regions with 5% or higher sequence divergence. Only use this preset if the average divergence is far below 5%.
    Asm5,
    /// Long assembly to reference mapping (-k19 -w19 -U50,500 --rmq -r100k -g10k -A1 -B9 -O16,41 -E2,1 -s200 -z200 -N50). Up to 10% sequence divergence.
    Asm10,
    /// Long assembly to reference mapping (-k19 -w10 -U50,500 --rmq -r100k -g10k -A1 -B4 -O6,26 -E2,1 -s200 -z200 -N50). Up to 20% sequence divergence.
    Asm20,
    /// Long-read spliced alignment (-k15 -w5 --splice -g2k -G200k -A1 -B2 -O2,32 -E1,0 -b0 -C9 -z200 -ub --junc-bonus=9 --cap-sw-mem=0 --splice-flank=yes). In the splice mode, 1) long deletions are taken as introns and represented as the ‘N’ CIGAR operator; 2) long insertions are disabled; 3) deletion and insertion gap costs are different during chaining; 4) the computation of the ‘ms’ tag ignores introns to demote hits to pseudogenes.
    Splice,
    /// Long-read splice alignment for PacBio CCS reads (-xsplice -C5 -O6,24 -B4).
    SpliceHq,
    /// Short single-end reads without splicing (-k21 -w11 --sr --frag=yes -A2 -B8 -O12,32 -E2,1 -b0 -r100 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g100 -2K50m --heap-sort=yes --secondary=no).
    ShortRead,
    /// PacBio CLR all-vs-all overlap mapping (-Hk19 -Xw5 -e0 -m100).
    AvaPb,
    /// Oxford Nanopore all-vs-all overlap mapping (-k15 -Xw5 -e0 -m100 -r2k).
    AvaOnt,
}

impl Preset {
    /// Get the preset name as a null-terminated byte literal. Intended for use with minimap2's `mm_set_opt` function.
    pub fn as_bytes(&self) -> &[u8] {
        match self {
            Preset::MapOnt => b"map-ont\0",
            Preset::MapHifi => b"map-hifi\0",
            Preset::MapPb => b"map-pb\0",
            Preset::LongReadHq => b"lr:hq\0",
            Preset::Asm5 => b"asm5\0",
            Preset::Asm10 => b"asm10\0",
            Preset::Asm20 => b"asm20\0",
            Preset::Splice => b"splice\0",
            Preset::SpliceHq => b"splice:hq\0",
            Preset::ShortRead => b"sr\0",
            Preset::AvaPb => b"ava-pb\0",
            Preset::AvaOnt => b"ava-ont\0",
        }
    }
}
