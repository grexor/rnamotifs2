"""
Vectorised motif search — a fast path for rnamotifs2.search.

`coverage_batch()` reproduces, for a whole batch of sequences at once, exactly
what `rnamotifs2.search.coverage()` computes one sequence at a time:

    v1[i] = 1 if position i is covered by any occurrence of the motif(s)
    v2    = np.convolve(v1, [1]*(2*hw+1), "same")          (windowed count)
    v     = v1 * v2
    result = v[hw : len(seq)-hw]

`v17()` is a vectorised re-implementation of `rnamotifs2.search.v17` returning
the identical 6-tuple. Per-region sequence encoding and k-mer packing are cached
so a pool worker pays that cost once and then scores each motif with a single
array compare. A golden test (`tests/`) checks it matches the reference.
"""

import numpy as np
from collections import Counter

import rnamotifs2
import pybio

# nucleotide -> small int; anything else (N, "-", padding) -> 4, which can never
# equal a motif base code (0..3). Motif non-ACGT chars are coded 255.
_SEQ_CODE = np.full(256, 4, dtype=np.uint8)
for _b, _c in ((65, 0), (67, 1), (71, 2), (84, 3)):
    _SEQ_CODE[_b] = _c
    _SEQ_CODE[_b + 32] = _c

_BASE = 5  # radix for packing k-mers (codes 0..4)


def expand_motif(m):
    return pybio.sequence.expand(m)


def _encode_seqs(seqs, Lmax):
    n = len(seqs)
    S = np.full((n, Lmax), 4, dtype=np.uint8)
    lens = np.zeros(n, dtype=np.int64)
    for i, s in enumerate(seqs):
        if not s:
            continue
        a = np.frombuffer(s.encode("latin-1", "replace"), dtype=np.uint8)
        lens[i] = a.shape[0]
        S[i, :a.shape[0]] = _SEQ_CODE[a]
    return S, lens


def _pack(S, k):
    """(n, Lmax-k+1) int32: base-5 code of the k-mer starting at each position."""
    n, Lmax = S.shape
    W = Lmax - k + 1
    if W <= 0:
        return np.zeros((n, 0), dtype=np.int32)
    P = S[:, 0:W].astype(np.int32)
    for j in range(1, k):
        P = P * _BASE + S[:, j:j + W]
    return P


def _motif_code(em):
    c = 0
    for ch in em:
        b = _SEQ_CODE[ord(ch)]
        c = c * _BASE + (int(b) if b < 4 else 4)
    # a motif containing a non-ACGT char gets a code that real k-mers (bases
    # 0..3 only) can never produce, so it never matches — same as the
    # reference, where str.find of such a motif never hits.
    return c


def _starts_from_pack(packs, em):
    k = len(em)
    P = packs.get(k)
    if P is None or P.size == 0:
        n = packs["n"]
        return np.zeros((n, 0), dtype=bool)
    return P == _motif_code(em)


def _cover_from_starts(starts_list, Lmax):
    """OR the [p, p+k) spans of every start array into v1 (n, Lmax) bool."""
    n = None
    for st, k in starts_list:
        if st.shape[1]:
            n = st.shape[0]
            break
    if n is None:
        # figure n from any entry
        n = starts_list[0][0].shape[0] if starts_list else 0
    v1 = np.zeros((n, Lmax), dtype=bool)
    for st, k in starts_list:
        W = st.shape[1]
        if W:
            for d in range(k):
                v1[:, d:d + W] |= st
    return v1


def _window_sum(v1i, hw):
    n, Lmax = v1i.shape
    csum = np.zeros((n, Lmax + 1), dtype=np.int32)
    np.cumsum(v1i, axis=1, out=csum[:, 1:])
    idx = np.arange(Lmax)
    lo = np.maximum(idx - hw, 0)
    hi = np.minimum(idx + hw + 1, Lmax)
    return csum[:, hi] - csum[:, lo]


def _packs_for(S):
    return {"n": S.shape[0], 3: _pack(S, 3), 4: _pack(S, 4), 5: _pack(S, 5)}


def coverage_batch(seqs, hw, motif, strict=False, packs=None, S=None, lens=None):
    """Vectorised `rnamotifs2.search.coverage` for many sequences.

    Returns (V, lens); the caller slices row i as ``V[i, hw:lens[i]-hw]``.
    """
    motif_list = motif.split("_") if isinstance(motif, str) else list(motif)

    if S is None:
        Lmax = max((len(s) for s in seqs), default=0)
        S, lens = _encode_seqs(seqs, Lmax)
        packs = _packs_for(S)
    Lmax = S.shape[1]
    n = S.shape[0]
    if Lmax == 0:
        return np.zeros((n, 0), dtype=np.int32), lens

    starts_list = []
    row_ok = np.ones(n, dtype=bool)
    for m in motif_list:
        for em in expand_motif(m):
            k = len(em)
            if k not in packs:
                packs[k] = _pack(S, k)
            st = _starts_from_pack(packs, em)
            if strict:
                row_ok &= st.any(axis=1) if st.size else np.zeros(n, dtype=bool)
            starts_list.append((st, k))

    v1 = _cover_from_starts(starts_list, Lmax)
    if strict:
        v1[~row_ok, :] = False
    v1i = v1.astype(np.int32)
    return v1i * _window_sum(v1i, hw), lens


# ---------------------------------------------------------------------------
# vectorised v17
# ---------------------------------------------------------------------------

_R1 = slice(251 - 95, 251 - 55 + 1)
_R2A = slice(251 - 50, 251 - 20 + 1)
_R2B = slice(20, 50 + 1)
_R3 = slice(60, 100 + 1)

_region_cache = {}


def _region_data(region, hw):
    """All events of class {region[-1], c}, encoded. Independent of rfilter —
    the caller masks rows per call."""
    key = (region, id(rnamotifs2.data.data), hw, rnamotifs2.config.perms,
           rnamotifs2.perm.generation)
    hit = _region_cache.get(key)
    if hit is not None:
        return hit

    want = region[-1]
    rows, meta = [], []
    a2_raw, a3_raw, a2_seqs, a3_seqs = [], [], [], []
    nums_all = Counter()

    for (eid, chr, strand, skip_start, in_start, in_stop, skip_stop, event_class) in rnamotifs2.data.data:
        if event_class not in (want, "c"):
            continue
        ec = {"c": "c", "e": "t", "s": "t"}[event_class]
        nums_all[ec + ".all"] += 1
        (_, _, _, _), (a2_exon, a2_intron, _, _), (a3_exon, a3_intron, _, _), (_, _, _, _) = \
            rnamotifs2.sequence.coords(strand, skip_start, in_start, in_stop, skip_stop)
        s1, s2, s3, s4 = rnamotifs2.sequence.sequence[eid]
        rows.append((eid, ec))
        meta.append((a2_intron, a2_exon, a3_intron, a3_exon))
        a2_raw.append(s2)
        a3_raw.append(s3)
        a2_seqs.append("-" * (200 - a2_intron) + s2[hw:-hw] + "-" * (50 - a2_exon))
        a3_seqs.append("-" * (50 - a3_exon) + s3[hw:-hw] + "-" * (200 - a3_intron))

    need_a2 = region in ("r1s", "r1e", "r2s", "r2e")
    need_a3 = region in ("r2s", "r2e", "r3s", "r3e")
    Sa2 = Sa3 = Pa2 = Pa3 = None
    la2 = la3 = None
    if need_a2:
        Sa2, la2 = _encode_seqs(a2_raw, max((len(s) for s in a2_raw), default=0))
        Pa2 = _packs_for(Sa2)
    if need_a3:
        Sa3, la3 = _encode_seqs(a3_raw, max((len(s) for s in a3_raw), default=0))
        Pa3 = _packs_for(Sa3)

    is_t = np.array([ec == "t" for (_, ec) in rows], dtype=bool)
    meta_arr = np.array(meta, dtype=np.int64) if meta else np.zeros((0, 4), np.int64)

    data = dict(rows=rows, meta=meta_arr, is_t=is_t, nums_all=nums_all,
                a2_seqs=a2_seqs, a3_seqs=a3_seqs,
                Sa2=Sa2, Sa3=Sa3, Pa2=Pa2, Pa3=Pa3, la2=la2, la3=la3,
                need_a2=need_a2, need_a3=need_a3)

    # permutation label matrix, built once per region and reused for every
    # motif: column k tells you, for permutation p, whether the event at rank
    # k among a motif's "considered" events (see v17) is a random draw of the
    # target class or of "c" - see the long comment in v17 for why "rank"
    # rather than "row" is what rnamotifs2.perm.ec_perm actually encodes.
    perms = rnamotifs2.config.perms
    if perms > 0:
        letter = region[-1]
        n_region = len(rows)
        dc_arr = np.asarray(rnamotifs2.data.data_class)
        is_target = np.zeros((perms, n_region), dtype=bool)
        is_control = np.zeros((perms, n_region), dtype=bool)
        for p in range(perms):
            labels = dc_arr[rnamotifs2.perm.ec_perm[p][:n_region]]
            is_target[p] = labels == letter
            is_control[p] = labels == "c"
        data["perm_is_target"] = is_target
        data["perm_is_control"] = is_control

    _region_cache.clear()  # only ever one comparison/region live at a time
    _region_cache[key] = data
    return data


def _region_signal(region, rd, motif, hw):
    """(n, reg_len) int array: the region vector `r` for every kept event."""
    n = len(rd["rows"])
    Va2 = coverage_batch(None, hw, motif, S=rd["Sa2"], packs=rd["Pa2"],
                         lens=rd["la2"])[0] if rd["need_a2"] else None
    Va3 = coverage_batch(None, hw, motif, S=rd["Sa3"], packs=rd["Pa3"],
                         lens=rd["la3"])[0] if rd["need_a3"] else None

    reg_len = 62 if region in ("r2s", "r2e") else 41
    out = np.zeros((n, reg_len), dtype=np.int64)
    meta = rd["meta"]
    for i in range(n):
        a2_intron, a2_exon, a3_intron, a3_exon = meta[i]
        if rd["need_a2"]:
            cov2 = Va2[i, hw:rd["la2"][i] - hw]
            p2 = np.zeros(251, dtype=np.int64)
            p2[200 - a2_intron:200 - a2_intron + cov2.shape[0]] = cov2
        if rd["need_a3"]:
            cov3 = Va3[i, hw:rd["la3"][i] - hw]
            p3 = np.zeros(251, dtype=np.int64)
            p3[50 - a3_exon:50 - a3_exon + cov3.shape[0]] = cov3
        if region in ("r1s", "r1e"):
            out[i] = p2[_R1]
        elif region in ("r2s", "r2e"):
            out[i] = np.concatenate([p2[_R2A], p3[_R2B]])
        else:
            out[i] = p3[_R3]
    return out


def v17(comps, genome, region="r1s", motif="YCAY", hw=15, pth=4, rfilter={},
        step=0, base_motif="TTTT", base_h=None):
    rnamotifs2.data.read_config(comps)
    rd = _region_data(region, hw)
    rows, is_t = rd["rows"], rd["is_t"]
    n = len(rows)

    keep = np.array(
        [rfilter.get("%s.%s" % (ec, eid), None) is None for (eid, ec) in rows],
        dtype=bool)

    nums = Counter(rd["nums_all"])
    nums["t"] = int((keep & is_t).sum())
    nums["c"] = int((keep & ~is_t).sum())

    rmat = _region_signal(region, rd, motif, hw)
    rmax = rmat.max(axis=1) if n else np.zeros(0, dtype=np.int64)

    # ---- consider_exons --------------------------------------------------
    if step == 0:
        consider = keep.copy()
    else:
        consider = np.zeros(n, dtype=bool)
        basemotif_hmin = max(4, base_h // 2)
        a2s, a3s = rd["a2_seqs"], rd["a3_seqs"]
        if region in ("r1s", "r1e"):
            win = [s[251 - 95 - 15:251 - 55 + 1 + 15] for s in a2s]
        elif region in ("r2s", "r2e"):
            win = [a2s[i][251 - 50 - 15:251 - 20 + 1] + a3s[i][20:50 + 1 + 15] for i in range(n)]
        else:
            win = [s[60 - 15:100 + 1 + 15] for s in a3s]
        Wcov, _ = coverage_batch(win, hw, [base_motif, motif[0]], strict=True)
        for i in range(n):
            if not keep[i]:
                continue
            L = len(win[i])
            if L > 30:
                rr = Wcov[i, hw:L - hw]
                if rr.size and rr.max() >= basemotif_hmin:
                    consider[i] = True
                    nums["t1" if rows[i][1] == "t" else "c1"] += 1

    # ---- choose h ------------------------------------------------------------
    if step == 0:
        tmax = rmax[keep & is_t]
        distances = []
        for h in range(4, 32):
            perc = int((tmax >= h).sum()) / float(nums["t"]) * 100
            distances.append((h, perc, abs(perc - pth)))
        chosen_h = rnamotifs2.search.choose_h(
            distances, perc_from=rnamotifs2.data.perc_from,
            perc_to=rnamotifs2.data.perc_to)
        if chosen_h is None:
            print("no h found, ignoring motif %s" % motif)
            return None
    else:
        thr = max(2, int(0.04 * nums["t"]))
        thr50 = 0.5 * nums["t1"]
        cmax = rmax[consider & is_t]
        distances = []
        for h in range(4, 32):
            ep = int((cmax >= h).sum())
            if ep < thr50:
                distances.append((abs(ep - thr), h, ep, thr))
        distances.sort()
        if not distances:
            return None
        chosen_h = max(4, distances[0][1])

    print("%s.%s.%s: h=%s" % (comps, genome, motif, chosen_h))

    # ---- filter + sum + count -----------------------------------------------
    vectors_sum = {}
    rcounts = Counter()
    present = {}
    fmat = (rmat >= chosen_h).astype(np.int64)
    for i, (eid, ec) in enumerate(rows):
        if not consider[i]:
            continue
        f = fmat[i]
        cur = vectors_sum.get(ec)
        vectors_sum[ec] = f.copy() if cur is None else cur + f
        if f.sum() > 0:
            rcounts[ec] += 1
            present["%s.%s" % (ec, eid)] = 1
            if rmax[i] >= 14:
                rfilter["%s.%s" % (ec, eid)] = 1

    # ---- permutation counts (rnamotifs2.config.perms > 0 only) ---------------
    # rcounts["<letter>.p<p>"] / rcounts["c.p<p>"] are what compute.rtest reads
    # back to build p_emp. See the "perm_is_target"/"perm_is_control" comment
    # in _region_data for the rank-not-row indexing this replicates.
    if "perm_is_target" in rd:
        considered_sum_r = (fmat[consider].sum(axis=1) > 0).astype(np.int64)
        n_considered = considered_sum_r.shape[0]
        letter = region[-1]
        IT = rd["perm_is_target"][:, :n_considered].astype(np.int64)
        IC = rd["perm_is_control"][:, :n_considered].astype(np.int64)
        target_counts = IT.dot(considered_sum_r)
        control_counts = IC.dot(considered_sum_r)
        for p in range(target_counts.shape[0]):
            rcounts["%s.p%s" % (letter, p)] = int(target_counts[p])
            rcounts["c.p%s" % p] = int(control_counts[p])

    return vectors_sum, rcounts, chosen_h, rfilter, nums, present


# ---------------------------------------------------------------------------
# vectorised areas() (used by the drawing step)
# ---------------------------------------------------------------------------

_areas_cache = {}


def _areas_data(hw):
    key = (id(rnamotifs2.data.data), hw)
    hit = _areas_cache.get(key)
    if hit is not None:
        return hit

    rows = []
    raw = [[], [], [], []]        # r1..r4 genomic seqs
    off = []                      # (left pad for r1..r4)
    nums = Counter()
    for (eid, chr, strand, skip_start, in_start, in_stop, skip_stop, event_class) in rnamotifs2.data.data:
        nums[event_class] += 1
        (r1_exon, r1_intron, _, _), (r2_exon, r2_intron, _, _), \
            (r3_exon, r3_intron, _, _), (r4_exon, r4_intron, _, _) = \
            rnamotifs2.sequence.coords(strand, skip_start, in_start, in_stop, skip_stop)
        s = rnamotifs2.sequence.sequence[eid]
        for j in range(4):
            raw[j].append(s[j])
        rows.append((eid, event_class))
        off.append((50 - r1_exon, 200 - r2_intron, 50 - r3_exon, 200 - r4_intron))

    enc = []
    for j in range(4):
        S, lens = _encode_seqs(raw[j], max((len(x) for x in raw[j]), default=0))
        enc.append((S, lens, _packs_for(S)))
    data = dict(rows=rows, off=np.array(off, dtype=np.int64), nums=nums, enc=enc)
    _areas_cache.clear()
    _areas_cache[key] = data
    return data


def areas(comps, motif="YCAY", hw=15, h=None, pth=4):
    rnamotifs2.data.read_config(comps)
    ad = _areas_data(hw)
    rows, off, nums = ad["rows"], ad["off"], ad["nums"]
    n = len(rows)

    # (n, 4, 251) padded coverage
    P = np.zeros((n, 4, 251), dtype=np.int64)
    for j in range(4):
        S, lens, packs = ad["enc"][j]
        V = coverage_batch(None, hw, motif, S=S, packs=packs, lens=lens)[0]
        for i in range(n):
            cov = V[i, hw:lens[i] - hw]
            lo = off[i, j]
            P[i, j, lo:lo + cov.shape[0]] = cov

    is_ec = np.array([ec for (_, ec) in rows])
    pmax = P.max(axis=2)  # (n, 4)

    if h is None:
        print("%s.%s: looking for h closest to threshold" % (comps, motif))
        notc = is_ec != "c"
        m = pmax[notc].max(axis=1)
        denom = float(nums["e"] + nums["s"])
        distances = []
        for hh in range(4, 32):
            perc = int((m >= hh).sum()) / denom * 100
            distances.append((hh, perc, abs(perc - pth)))
        h = rnamotifs2.search.choose_h(
            distances, perc_from=rnamotifs2.data.perc_from,
            perc_to=rnamotifs2.data.perc_to)

    print("%s.%s: h=%s" % (comps, motif, h))

    F = (P >= h).astype(np.int64)                 # (n, 4, 251)
    stats = Counter()
    v2, v3 = F[:, 1, :], F[:, 2, :]
    r1p = v2[:, 251 - 95:251 - 55 + 1].sum(axis=1) > 0
    r2p = (v2[:, 251 - 50:251 - 20 + 1].sum(axis=1) + v3[:, 20:50 + 1].sum(axis=1)) > 0
    r3p = v3[:, 60:100 + 1].sum(axis=1) > 0
    for i, (eid, ec) in enumerate(rows):
        stats[ec] += 1
        stats["r1%s" % ec] += int(r1p[i])
        stats["r2%s" % ec] += int(r2p[i])
        stats["r3%s" % ec] += int(r3p[i])

    vectors_sum = {}
    for ec in set(is_ec.tolist()):
        sel = is_ec == ec
        s = F[sel].sum(axis=0)  # (4, 251)
        vectors_sum[ec] = (s[0], s[1], s[2], s[3])
    return vectors_sum, h, stats
