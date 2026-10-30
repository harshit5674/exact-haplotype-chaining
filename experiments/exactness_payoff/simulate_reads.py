"""Simulate HiFi-like reads that are mosaics of 1-3 haplotypes of a panel (graph.gfa.gz walks).

Model (seeded, deterministic):
- Read length L ~ Uniform{20000..40000}. k haplotypes (k = 1 controls, 2 or 3 mosaics).
- Breakpoints: k-1 cut points in the error-free mosaic, uniform in [4000, L-4000] with spacing >= 5000.
- Haplotype order: h_1 uniform; h_{s+1} uniform among the other haplotypes (no conditioning on alleles).
- A switch from h to h' at hap position q of h: if the vertex containing q also lies on h', the read continues
  on h' at the same vertex and offset; otherwise it moves forward along h to the next vertex shared with h'
  and switches at its first base (the cut moves right). Cuts that would move more than 8 kb are resampled.
  So every read is a valid walk of the graph, and the true haplotype of every read base is known.
- Errors on the mosaic, per base: substitution 0.10%, insertion 0.05%, deletion 0.05%; indel length 1
  (prob 0.6) or 2 (prob 0.4); total about 0.2% events. Bases uniform over ACGT (substitutions differ).
- Strand: '+' or '-' with probability 1/2; '-' reads are written reverse-complemented.
- Read name: sim<idx>_k<k>_<strand>; the truth file has everything else.
- Informative mode (--min-walk-switches m): mosaics (k >= 2) are kept only if their own vertex walk,
  restricted to vertices of at least 17 bp, cannot be explained with fewer than m haplotype switches
  (intersect-or-reset scan over the panel's vertex memberships). Controls (k = 1) are unaffected.

Outputs: <prefix>.fq.gz and <prefix>.truth.json (per read: haplotype segments in mosaic coordinates,
breakpoints in mosaic and final forward-read coordinates, and the error events, which reconstruct the
read-to-mosaic coordinate map). Coordinates are 0-based; "forward read" means before reverse complement.

  python3.13 experiments/exactness_payoff/simulate_reads.py --panel hprc60 --prefix .../sim/hprc60 --seed 20261006
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import json
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'experiments/hprc_mhc_panel'))
import convert  # noqa: E402  (pure loader: load_graph; no side effects on import)

PANELS = {name: ROOT / 'experiments/hprc_mhc_panel/data' / name for name in ('hprc60', 'hprc463')}
RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
SUB, INS, DEL = 0.0010, 0.0005, 0.0005
MAX_SHIFT = 8000


class Haplotypes:
    def __init__(self, panel_dir: Path):
        graph, seqs, ids = convert.load_graph(panel_dir / 'graph.gfa.gz')
        self.graph = graph
        self.seq = [seqs[name] for name in graph['names']]  # by topological vertex id
        self.paths = graph['paths']
        self.names = graph['haplotypes']
        self.H = len(self.paths)
        self.starts, self.length, self.index = [], [], []
        for p in self.paths:
            st, pos = [], 0
            for v in p:
                st.append(pos)
                pos += len(self.seq[v])
            self.starts.append(st)
            self.length.append(pos)
            self.index.append({v: k for k, v in enumerate(p)})

    def locate(self, h: int, q: int) -> tuple[int, int, int]:
        """Hap position q -> (walk index, vertex, offset in vertex)."""
        k = bisect.bisect_right(self.starts[h], q) - 1
        v = self.paths[h][k]
        return k, v, q - self.starts[h][k]

    def substring(self, h: int, a: int, b: int) -> str:
        out, q = [], a
        while q < b:
            k, v, off = self.locate(h, q)
            take = min(b - q, len(self.seq[v]) - off)
            out.append(self.seq[v][off:off + take])
            q += take
        return ''.join(out)


def plan_read(hs: Haplotypes, rng: random.Random, k: int):
    """Return (segments, cuts) of an error-free mosaic, or None to resample."""
    L = rng.randint(20000, 40000)
    if k == 1:
        cuts = []
    else:
        for _ in range(1000):
            cuts = sorted(rng.randint(4000, L - 4000) for _ in range(k - 1))
            if all(b - a >= 5000 for a, b in zip(cuts, cuts[1:])):
                break
        else:
            return None
    haps = [rng.randrange(hs.H)]
    for _ in range(k - 1):
        haps.append(rng.choice([h for h in range(hs.H) if h != haps[-1]]))
    h = haps[0]
    start = rng.randint(0, hs.length[h] - L - 20000)
    segs, pos, rpos = [], start, 0
    for s in range(k):
        h = haps[s]
        if s == k - 1:
            end = pos + (L - rpos)
            if end > hs.length[h] or L - rpos < 2000:
                return None
            segs.append(dict(read_from=rpos, read_to=rpos + end - pos, hap=h, hap_from=pos))
            break
        q = pos + (cuts[s] - rpos)
        if q >= hs.length[h]:
            return None
        widx, v, off = hs.locate(h, q)
        h2 = haps[s + 1]
        if v not in hs.index[h2]:
            # move forward along h to the next vertex shared with h2
            nxt = next((w for w in range(widx + 1, len(hs.paths[h])) if hs.paths[h][w] in hs.index[h2]), None)
            if nxt is None:
                return None
            q2 = hs.starts[h][nxt]
            if q2 - q > MAX_SHIFT:
                return None
            q, v, off = q2, hs.paths[h][nxt], 0
        seg_len = q - pos
        if seg_len < 2000:
            return None
        segs.append(dict(read_from=rpos, read_to=rpos + seg_len, hap=h, hap_from=pos))
        rpos += seg_len
        pos = hs.starts[h2][hs.index[h2][v]] + off
    return segs


def walk_min_switches(hs: Haplotypes, segs, masks, minlen: int = 17) -> int:
    """Min label changes needed to explain the mosaic's own vertex walk with single panel haplotypes,
    counting only vertices of length >= minlen (shorter ones cannot hold a 17-mer anchor)."""
    seq_masks = []
    for s in segs:
        h, a = s['hap'], s['hap_from']
        b = a + s['read_to'] - s['read_from']
        k0 = hs.locate(h, a)[0]
        k1 = hs.locate(h, b - 1)[0]
        for k in range(k0, k1 + 1):
            v = hs.paths[h][k]
            if len(hs.seq[v]) >= minlen:
                seq_masks.append(masks[v])
    if not seq_masks:
        return 0
    cur, resets = seq_masks[0], 0
    for m in seq_masks[1:]:
        if cur & m:
            cur &= m
        else:
            cur, resets = m, resets + 1
    return resets


def apply_errors(seq: str, rng: random.Random):
    """Return (read, events); events are (mosaic_pos, kind, bases) in increasing position."""
    out, events = [], []
    i = 0
    n = len(seq)
    while i < n:
        u = rng.random()
        if u < SUB:
            b = rng.choice([x for x in 'ACGT' if x != seq[i].upper()])
            out.append(b)
            events.append((i, 'S', b))
            i += 1
        elif u < SUB + INS:
            ln = 1 if rng.random() < 0.6 else 2
            ins = ''.join(rng.choice('ACGT') for _ in range(ln))
            out.append(seq[i])
            out.append(ins)  # inserted after base i
            events.append((i, 'I', ins))
            i += 1
        elif u < SUB + INS + DEL:
            ln = 1 if rng.random() < 0.6 else 2
            events.append((i, 'D', seq[i:i + ln]))
            i += ln
        else:
            out.append(seq[i])
            i += 1
    return ''.join(out), events


def read_to_mosaic(events, mosaic_len: int) -> list[int]:
    """Forward-read index -> mosaic index (-1 for inserted bases), rebuilt from the error events."""
    ev = {e[0]: e for e in events}
    out, i = [], 0
    while i < mosaic_len:
        e = ev.get(i)
        if e is None or e[1] == 'S':
            out.append(i)
            i += 1
        elif e[1] == 'I':
            out.append(i)
            out.extend([-1] * len(e[2]))
            i += 1
        else:  # deletion
            i += len(e[2])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--panel', choices=sorted(PANELS), default='hprc60')
    ap.add_argument('--prefix', required=True)
    ap.add_argument('--seed', type=int, default=20261006)
    ap.add_argument('--counts', default='80,110,110', help='reads with k = 1, 2, 3 haplotypes')
    ap.add_argument('--haps', default='', help='optional comma-separated haplotype indices to draw from')
    ap.add_argument('--min-walk-switches', type=int, default=0,
                    help='informative mode: reject mosaics (k >= 2) whose vertex walk needs fewer switches')
    args = ap.parse_args()
    prefix = Path(args.prefix)
    fq, truth_path = Path(str(prefix) + '.fq.gz'), Path(str(prefix) + '.truth.json')
    for p in (fq, truth_path):
        if p.exists():
            raise SystemExit(f'{p} exists')
    prefix.parent.mkdir(parents=True, exist_ok=True)
    hs = Haplotypes(PANELS[args.panel])
    vmasks = [0] * len(hs.seq)
    for h, path in enumerate(hs.paths):
        for v in path:
            vmasks[v] |= 1 << h
    allowed = [int(x) for x in args.haps.split(',')] if args.haps else None
    rng = random.Random(args.seed)
    counts = [int(x) for x in args.counts.split(',')]
    reads = []
    idx = tried = 0
    with gzip.open(fq, 'wt') as f:
        for k, cnt in zip((1, 2, 3), counts):
            made = 0
            while made < cnt:
                segs = plan_read(hs, rng, k)
                if segs is None:
                    continue
                if allowed is not None and any(s['hap'] not in allowed for s in segs):
                    continue
                wms = walk_min_switches(hs, segs, vmasks)
                tried += 1
                if k >= 2 and wms < args.min_walk_switches:
                    continue
                mosaic = ''.join(hs.substring(s['hap'], s['hap_from'], s['hap_from'] + s['read_to'] - s['read_from'])
                                 for s in segs).upper()
                read, events = apply_errors(mosaic, rng)
                rmap = read_to_mosaic(events, len(mosaic))
                assert len(rmap) == len(read)
                inv = {}
                for r, m in enumerate(rmap):
                    if m >= 0 and m not in inv:
                        inv[m] = r
                cuts = [s['read_from'] for s in segs[1:]]

                def fwd(m):  # first forward-read index at or after mosaic index m
                    while m not in inv:
                        m += 1
                    return inv[m]
                strand = '+' if rng.random() < 0.5 else '-'
                name = f'sim{idx:04d}_k{k}_{strand}'
                out_seq = read if strand == '+' else read.translate(RC)[::-1]
                f.write(f'@{name}\n{out_seq}\n+\n{"I" * len(out_seq)}\n')
                reads.append(dict(name=name, k=k, strand=strand, read_length=len(read), mosaic_length=len(mosaic),
                                  segments=[dict(s, hap_name=hs.names[s['hap']]) for s in segs],
                                  breakpoints_mosaic=cuts, breakpoints_read=[fwd(m) for m in cuts],
                                  walk_min_switches=wms,
                                  error_events=[list(e) for e in events]))
                idx += 1
                made += 1
    doc = dict(panel=args.panel, seed=args.seed, counts=counts, min_walk_switches=args.min_walk_switches,
               planned_mosaics_tried=tried, model=dict(sub=SUB, ins=INS, dele=DEL,
               indel_lengths={'1': 0.6, '2': 0.4}, length=[20000, 40000], min_cut=4000, min_spacing=5000,
               max_shift=MAX_SHIFT), haplotypes=hs.names, reads=reads)
    truth_path.write_text(json.dumps(doc) + '\n')
    print(json.dumps(dict(reads=len(reads), fastq=str(fq), truth=str(truth_path))))


if __name__ == '__main__':
    main()
