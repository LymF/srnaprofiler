# Formulas

Every quantity `classify` reports: the formula, the null it is tested against, and what it
cannot tell you. Meant to be read next to the output table.

Two conventions run through all of it. A read mapping to *N* places contributes 1/*N* to each,
so every count is a weighted sum rather than a line count. And every feature comes from the
reads of its own locus; the library as a whole enters only in `reads_rpm`.

---

## Calling a locus

Depth is weighted coverage, `depth[p] = Σ 1/NH` over reads spanning *p*. A locus starts as a
run of positions above a threshold given in reads per million:

```
min_depth = max(1, round(mincov_RPM × lib_size / 1e6))
```

`lib_size` is `Σ 1/NH` over the whole library, not the number of alignment lines. The two can
differ by a factor of three.

Each island grows by `pad` on both sides, so neighbours tolerate a gap of up to 2 × `pad`
before touching. Touching is not merging. Two islands join only if the size class holding more
than half their read mass is the same, and |`frac_top`ₐ − `frac_top`ᵦ| ≤ 0.5. Neighbours that
fail split the gap between them.

Profiles are compared by mass rather than by modal length, because with a handful of reads per
raw island the mode is noise.

Intervals from `--mirna-loci` are carved out before merging and become loci of their own, with
the flanks continuing separately. A hairpin is about 90 bp and a locus is 600–1000 bp, so an
uncarved precursor gets whatever class the rest of the locus has.

---

## Per-locus features

```
reads         = Σ 1/NH                      both strands
reads_rpm     = reads / lib_size × 1e6
distinct_seqs = number of distinct sequences
frac_top      = Σ 1/NH(+) / reads
unique_frac   = Σ 1/NH over reads with NH=1 / reads
multimap_frac = alignment lines with NH>1 / all alignment lines
```

`distinct_seqs` is the *n* of every binomial test below. It is a sample size, not an
abundance, and the two diverge in repeat-derived loci.

### Size

Lengths are read in the read's own orientation, so minus-strand reads are reverse-complemented
first.

```
seq_lens[L]    = Σ 1/NH over reads of length L,  dicermin ≤ L ≤ dicermax
frac_in_window = (reads − short − long) / reads
dicer_call     = argmax_L seq_lens[L]   if frac_in_window ≥ 0.8, else "N"
L1             = argmax_L seq_lens[L]
band(L)        = seq_lens[L−1] + seq_lens[L] + seq_lens[L+1]
size_purity    = band(L1) / Σ seq_lens
L2             = argmax seq_lens[L] over L with |L − L1| ≥ 2
dicer_call_2   = L2, reported only if band(L2) / Σ seq_lens ≥ 0.2
```

Purity is measured over a band because Dicer products scatter by ±1 nt. On a single length, a
clean 21-nt locus looks impure and a second peak appears almost everywhere — that is the
spread of one population, not a second pathway.

### 5′ nucleotide

`pct_1u`, `pct_5a` and `pct_10a` are the weighted fraction of reads whose first base is U,
whose first base is A, or whose tenth base is A.

### Ping-pong

The 5′–5′ offset spectrum between strands:

```
counts[d] = Σ_p  plus5[p] × minus5[p + d − 1],   d = 1..dmax

pingpong_frac    = counts[10] / Σ_{d=1..20} counts[d],   only if Σ ≥ 10 pairs
overlap_dup_frac = counts[dicer_call − 2] / Σ_{d=1..dmax} counts[d]
```

The Dicer duplex offset is *L* − 2, from the 2-nt 3′ overhang, so it moves with read length.

The ten-pair guard matters when reading the column: below ten pairs `pingpong_frac` is written
as 0.0, meaning *not measurable*, not *absent*. In a near single-stranded locus there are no
opposing pairs to count. Filter on `pingpong_pairs10` before averaging, or you measure the
fraction of single-stranded loci instead.

### Phasing

Only reads whose length equals `dicer_call` enter the register — a 26-nt read in a 21-nt phase
is noise.

```
register(read) = (5′ − locus_start) mod L,  minus-strand 5′ shifted by −2
phase_number   = cycles occupied in the dominant register
phase_ratio    = mass in the dominant register / total mass
phase_score    = phase_ratio × phase_number × ln(register mass × 20e6 / lib_size)
```

Mass is normalized to 20 million. Without it the score tracks sequencing depth, and any
threshold calibrated on one library becomes a threshold on that library's depth.

---

## Evidence scores (`s_*`)

Each is `−log10(p)`, capped at 300.

| column | test | null |
|---|---|---|
| `s_1u`, `s_5a` | one-sided binomial, *k* = fraction × `distinct_seqs` | base composition of the locus DNA |
| `s_fractop` | two-sided binomial on `frac_top` | 0.5 |
| `s_fractop_bid` | one-sided binomial on the minority strand | 0.10 |
| `s_pingpong` | Poisson, P(X ≥ counts[10]) | mean rate of offsets 1–20 excluding 10 |
| `s_phase` | hypergeometric, scan-corrected | see below |

The 5′ base null is the locus DNA, not a flat 0.25. It answers: if a cut fell anywhere in this
locus, how often would the product start with this base?

```
null(base) = frac_top × f(base in DNA) + (1 − frac_top) × f(complement in DNA)
```

with *f* over every position of the locus, N excluded from the denominator. The `frac_top`
weighting is there because a minus-strand read starts with T where the genome has A. Real
values land near 0.27 in Drosophila and 0.33 in Arabidopsis, so 0.25 was optimistic in both.

The ping-pong null is the locus's own background. Offsets are not uniform, so the rate comes
from offsets 1–20 with the peak excluded.

Phasing is scanned, so it is corrected. An 11-cycle window slides one cycle at a time:

```
population = 2 × window_length − 1
successes  = 2 × cycles − 1
draws      = distinct occupied positions in the window
observed   = phase_number of the dominant register
p_window   = P(X ≥ observed)
p          = min(1, p_window × windows tested)        Bonferroni
```

The winning window is chosen by evidence — cycles in register, then score — never by p-value,
since the hypergeometric rewards an empty window. Reporting the best of many windows without
correcting inflates significance with locus length, and the effect is large: under a
permutation null, a long locus with no phase at all passed p ≤ 0.01 most of the time.

---

## Effect strengths (`e_*`)

These are what `score` uses: the same *k*, *n* and *p₀* as the `s_*`, read as effect size
rather than significance.

```
e_X = max(0, (L95(p̂, n) − p₀) / (1 − p₀))
```

L95 is the 95% Wilson lower bound. The result lies in [0, 1]: 0.61 means the effect is 61% of
the way from the null to the maximum, with 95% confidence. These are not TSV columns — they
feed `score`.

| strength | p̂ | n | p₀ |
|---|---|---|---|
| `e_1u`, `e_5a` | `pct_1u`, `pct_5a` | `distinct_seqs` | DNA composition |
| `e_fractop_bid` | minority strand | `distinct_seqs` | 0.10 |
| `e_pingpong` | `counts[10] / Σcounts` | `Σcounts` | background rate / `Σcounts` |
| `e_phase` | cycles in register / 11 | 11 | 1 / `phase_length` |
| `e_size` | class-window mass / windowed mass | `distinct_seqs` | window width / kept range |

`e_size` is the only class-dependent one. The rest mean the same thing for any class; size
means nothing except against a window.

A p-value mixes effect with sample size, which is why the score is built this way. The old
`−log10(p)` score ranked deep loci above strong ones.

---

## Classification

The class score uses the raw features, not the `s_*` columns.

```
_window(size, lo, hi, edge=3) = 1 inside the window, decaying 1 − d/edge outside
_mid(frac_top, σ=0.15)        = exp(−(frac_top − 0.5)² / 2σ²)
_sigmoid(x, x0, scale)        = 1 / (1 + exp(−(x − x0)/scale))

S(class) = w_size             × _window(dicer_call, lo, hi)
         + w_1u               × pct_1u
         + w_5a               × pct_5a
         + w_fractop_mid      × _mid(frac_top)
         + w_pingpong         × pingpong_frac × gate
         − w_pingpong_penalty × pingpong_frac × gate
         + w_phase            × phase_signal  × gate
```

`gate` is `_window(dicer_call, lo, hi)` of the class being scored, and every term is multiplied
by it. A signature cannot pull a locus into a class whose length range it does not occupy.
That is what lets one engine serve plants and arthropods: 5′ U is common to several small RNA
classes and means piRNA only together with the size.

`phase_signal` depends on the class's `phase_on`:

```
phase_on: score  (default) → _sigmoid(phase_score, phase_x0, phase_scale)
phase_on: cycles           → e_phase, the run in register with no abundance in it
```

Both modes exist because neither works everywhere. `phase_score` carries ln(mass), so against
a class that dominates the genome it selects depth rather than phase. But dropping abundance
loses real loci where the positives happen to be abundant. Hence one mode per class.

```
cls        = argmax_class S
confidence = softmax(S)[cls]
margin     = S(1st) − S(2nd)
review     = confidence < review_confidence OR margin < min_margin
             OR reads < min_reads OR dicer_call == "N"
score      = min(e_* over defining_evidence)
```

In `defining_evidence` a nested list is OR and takes the maximum; top-level entries are AND
and take the minimum. `[[1u, pingpong], size]` reads as 5′ U or ping-pong, and size. The
minimum is deliberate — one strong signal should not cover for a missing one.

`score` is 0 when `dicer_call` is "N", and NA for miRNA, whose class comes from annotation and
was never scored.

---

## Limits

`score` and the `s_*` columns answer different questions. The `s_*` ask whether a signature
would be surprising under the null, which grows with depth. `score` asks how strong the effect
is, which does not. Compare loci with `score`.

Score scales differ by class. Where a class is defined by one weak signature against a high
null, its scores sit low across the board — that is the score reporting weaker evidence, not a
calibration problem.

`s_fractop_bal` is obsolete. It tested equivalence to 50/50 when the signature is
bidirectionality, and a perfectly balanced locus could still fail it for lack of power. The
column stays for compatibility; `s_fractop_bid` replaced it.

`s_fractop_bid` is noisy. Single-stranded loci routinely clear the 0.10 null, and the
separation from double-stranded ones is under 2×.

`s_phase` is insensitive — it misses most annotated phased loci at p ≤ 0.01. The test itself is
sound and the scan correction is right; what does not work is ranking loci by that p-value.
That job belongs to `e_phase`.

`s_pingpong` assumes variance equals the mean. In loci with overdispersed offset backgrounds
it is anticonservative; a negative binomial would be the fix.

`e_phase` measures how far the phase runs, not how surprising it is, and it saturates on
occupancy rather than concentration. A locus carrying one read at every position of a window
reaches the maximum with no enrichment at all and a phasing p-value of 1.0. For `phasiRNA`,
where `score` *is* `e_phase`, that is worth knowing.

`e_pingpong` is optimistic: its *n* counts pairs, which scale with the square of coverage and
are not independent trials.

`phase_score` correlates with log(reads), so it is not comparable between loci. `pingpong_z` is
capped and unstable in narrow, deep loci — classify on `pingpong_frac`.

`overlap_dup_frac`, `size_purity`, `dicer_call_2` and `unique_frac` are reported and read by no
rule. Two of them discriminate well, and they still should not become class weights: they are
consequences of the class rather than evidence for it. piRNAs live in transposons, so their
`unique_frac` is near zero. Weighting that would teach the classifier that repetitive means
piRNA, which is a fact about where piRNAs live, not proof that a locus is one.

Scores from runs before the effect-strength change do not compare with current ones.
