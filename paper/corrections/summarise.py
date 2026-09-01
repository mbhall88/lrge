#!/usr/bin/env python3
"""Summary comparisons for the PR #45 benchmark, read from the collected TSV."""
import csv, sys, collections

path = sys.argv[1]
rows = list(csv.DictReader(open(path), delimiter="\t"))
by = {(r["accession"], r["mode"]): r for r in rows}
accs = sorted({r["accession"] for r in rows})

hist = {}
for r in csv.DictReader(open("/scratch/user/uqmhal11/lrge/paper/corrections/rerun_estimates.tsv"), delimiter="\t"):
    hist[r["run"]] = r


def ratio(acc, mode):
    r = by.get((acc, mode))
    if not r or r["estimate_over_truth"] == "NA":
        return None
    return float(r["estimate_over_truth"])


def hist_ratio(acc, col):
    h = hist.get(acc)
    if not h or h[col] in ("NA", ""):
        return None
    return int(h[col]) / int(h["truth_size"])


print("== per accession: old default (rerun_default) -> auto, plus -F before/after ==")
print(f"{'accession':<13}{'class':<9}{'old':>7}{'auto':>8}{'oldF':>8}{'autoF':>8}{'skew':>9}{'kept%':>8}")
for acc in accs:
    cls = by.get((acc, "auto"), {}).get("class", "?")
    old, auto = hist_ratio(acc, "rerun_default"), ratio(acc, "auto")
    oldF, autoF = hist_ratio(acc, "rerun_F"), ratio(acc, "auto_F")
    r = by.get((acc, "auto"), {})
    kept = "NA"
    if r.get("retained_reads", "NA") not in ("NA", "") and r.get("total_reads", "NA") not in ("NA", ""):
        kept = f"{100 * int(r['retained_reads']) / int(r['total_reads']):.1f}"
    fmt = lambda v: f"{v:.3f}" if v is not None else "  NA "
    print(f"{acc:<13}{cls:<9}{fmt(old):>7}{fmt(auto):>8}{fmt(oldF):>8}{fmt(autoF):>8}"
          f"{r.get('skew_score','NA'):>9}{kept:>8}")

print("\n== band crossings (outliers only, old default vs auto) ==")
outl = [a for a in accs if by.get((a, "auto"), {}).get("class") == "outlier"]
for label, lo, hi in [(">=0.5x", 0.5, None), (">=0.8x", 0.8, None), ("0.8-1.2x", 0.8, 1.2)]:
    def count(f):
        n = 0
        for a in outl:
            v = f(a)
            if v is None:
                continue
            if v >= lo and (hi is None or v <= hi):
                n += 1
        return n
    print(f"{label:<10} old {count(lambda a: hist_ratio(a,'rerun_default'))}/{len(outl)}"
          f"   auto {count(lambda a: ratio(a,'auto'))}/{len(outl)}")

print("\n== mean absolute relative error, outliers ==")
for label, f in [("old default", lambda a: hist_ratio(a, "rerun_default")),
                 ("auto", lambda a: ratio(a, "auto")),
                 ("old -F", lambda a: hist_ratio(a, "rerun_F")),
                 ("auto -F", lambda a: ratio(a, "auto_F"))]:
    vals = [abs(f(a) - 1) for a in outl if f(a) is not None]
    print(f"{label:<12} n={len(vals):<3} mean |est/truth - 1| = {sum(vals)/len(vals):.3f}")

print("\n== benefit of -F, before vs after normalization (outliers) ==")
print(f"{'accession':<13}{'old F/default':>15}{'auto F/auto':>14}")
for a in outl:
    old, oldF, auto, autoF = hist_ratio(a, "rerun_default"), hist_ratio(a, "rerun_F"), ratio(a, "auto"), ratio(a, "auto_F")
    b1 = f"{oldF/old:.2f}x" if old and oldF else "NA"
    b2 = f"{autoF/auto:.2f}x" if auto and autoF else "NA"
    print(f"{a:<13}{b1:>15}{b2:>14}")

print("\n== controls: auto vs never ==")
for a in [x for x in accs if by.get((x, "auto"), {}).get("class") == "control"]:
    for mode in ("never", "auto"):
        r = by.get((a, mode))
        if r:
            print(f"{a:<13}{mode:<7} est={r['estimate_bp']:>10} ratio={r['estimate_over_truth']:>6} "
                  f"norm={r['normalization_engaged']:<4} noovl={r['queries_without_overlaps']:>4} "
                  f"iqr={r['interval_low_bp']}-{r['interval_high_bp']}")

print("\n== never vs historical rerun_default (exactness check) ==")
for a in accs:
    r = by.get((a, "never"))
    if not r:
        continue
    h = hist.get(a, {}).get("rerun_default", "NA")
    ok = "SAME" if r["estimate_bp"] == h else "DIFF"
    print(f"{a:<13} never={r['estimate_bp']:>10}  rerun_default={h:>10}  {ok}")

print("\n== runtime / peak memory ==")
agg = collections.defaultdict(list)
for r in rows:
    if r["elapsed"] == "NA" or r["max_rss_kb"] == "NA":
        continue
    parts = [float(p) for p in r["elapsed"].split(":")]
    secs = parts[0] * 60 + parts[1] if len(parts) == 2 else parts[0] * 3600 + parts[1] * 60 + parts[2]
    agg[r["mode"]].append((secs, int(r["max_rss_kb"]) / 1048576))
for mode, vals in sorted(agg.items()):
    s = [v[0] for v in vals]; m = [v[1] for v in vals]
    print(f"{mode:<8} n={len(vals):<3} wall median {sorted(s)[len(s)//2]:6.1f}s max {max(s):6.1f}s | "
          f"RSS median {sorted(m)[len(m)//2]:5.2f} GiB max {max(m):5.2f} GiB")

print("\n== per-accession auto vs never cost ==")
for a in accs:
    ra, rn = by.get((a, "auto")), by.get((a, "never"))
    if ra and rn and ra["max_rss_kb"] != "NA" and rn["max_rss_kb"] != "NA":
        print(f"{a:<13} never {rn['elapsed']:>10} {int(rn['max_rss_kb'])/1048576:5.2f} GiB   "
              f"auto {ra['elapsed']:>10} {int(ra['max_rss_kb'])/1048576:5.2f} GiB   "
              f"engaged={ra['normalization_engaged']}")
