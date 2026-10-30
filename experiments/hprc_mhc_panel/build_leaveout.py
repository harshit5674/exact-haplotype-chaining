"""Leave-out panel: all 463 usable haplotypes except CHM13 (the reads' own haplotype).

The reads are CHM13 PacBio reads, and every panel built by build_panel.py contains CHM13,
as Minichain's 60-haplotype MHC graph does. This panel removes it, so no haplotype in the
graph is the reads' source and chains must follow other haplotypes or switch. It reuses
build_panel.py's region, threading and checks unchanged and writes data/hprc462_noCHM13/.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_panel as bp  # noqa: E402


def main():
    out = HERE / 'panels_leaveout.json'
    if out.exists():
        raise SystemExit('refusing to overwrite')
    for n, src in bp.SOURCES.items():
        assert bp.digest(bp.DATA / n, 'md5') == src['md5'], n
    region = bp.locate_region()
    samples, bubbles = bp.read_calls(region)
    walks, status = bp.thread(samples, bubbles)
    del walks['CHM13#0']
    names = {s for w in walks.values() for s in w}
    seqs, links = bp.load_segments(names)
    panel = bp.write_panel('hprc462_noCHM13', walks, seqs, links)
    out.write_text(json.dumps(dict(script='experiments/hprc_mhc_panel/build_leaveout.py', region=region,
                                   removed='CHM13#0', panel=panel), indent=1) + '\n')
    print({k: v for k, v in panel.items() if k != 'haplotypes'})


if __name__ == '__main__':
    main()
