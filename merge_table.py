#!/usr/bin/env python3
import re
import sys
from pathlib import Path

import numpy as np

ETAMIN, ETAMAX, NETA = -12.0, 12.0, 121
NT = 301
HALF = (NETA - 1) // 2          # index of y = 0; only y >= 0 is computed

FILES = ["gluon_energy_density.dat"] + [f"quark_construction_n{k}_{q}.dat"
                                        for k in (0, 1) for q in "uds"]
ROWS = {f: (NT * NT if f.startswith("gluon") else NT) for f in FILES}


def check_symmetry(path):
    # gluon(+y,T1,T2) should equal gluon(-y,T2,T1): the integrand maps into
    # itself under that exchange, and T1 and T2 run over the same grid.
    g = np.loadtxt(path)
    key = {(round(y, 4), round(t1, 4), round(t2, 4)): v for y, t1, t2, v in g}
    worst = 0.0
    for (y, t1, t2), v in key.items():
        m = key.get((round(-y, 4), t2, t1))
        if m is not None and max(abs(v), abs(m)) > 0:
            worst = max(worst, abs(v - m) / max(abs(v), abs(m)))
    return worst


def mirror(name, text):
    # The negative rapidities are not computed. For the gluon table the mirror
    # is e(-y,T1,T2) = e(+y,T2,T1), so the (T1,T2) plane is transposed. For the
    # quark tables the row already carries f(+y) and f(-y) side by side, so the
    # mirror only swaps the two pairs of columns.
    rows = [l.split() for l in text.splitlines() if l.strip()]
    y = "%g" % -float(rows[0][0])
    if not name.startswith("gluon"):
        return "".join("\t".join([y, r[1], r[4], r[5], r[2], r[3]]) + "\n"
                       for r in rows)
    val = {(r[1], r[2]): r[3] for r in rows}
    return "".join("\t".join([y, r[1], r[2], val[(r[2], r[1])]]) + "\n"
                   for r in rows)


out = Path("tabs/IPSat/Set0")
out.mkdir(parents=True, exist_ok=True)
srcs = [Path(f"table_build/eta_{i:03d}/tabs/IPSat/Set0") for i in range(HALF, NETA)]

for name in FILES:
    texts = []
    for j, s in enumerate(srcs):
        text = (s / name).read_text()
        n = sum(1 for line in text.splitlines() if line.strip())
        if n != ROWS[name]:
            sys.exit(f"eta_{HALF+j:03d}/{name}: {n} rows, expected {ROWS[name]}")
        # make_cuhre_1C writes a zero when Cuhre returns NaN and says so in the
        # job log, but a NaN reaching the file some other way would read back as
        # a table the shape checks accept.
        low = text.lower()
        if "nan" in low or "inf" in low:
            sys.exit(f"eta_{HALF+j:03d}/{name}: non-finite value")
        texts.append(text)
    with (out / name).open("w") as fout:
        for text in reversed(texts[1:]):        # y = -12 .. -0.2, mirrored
            fout.write(mirror(name, text))
        for text in texts:                      # y = 0 .. 12, as computed
            fout.write(text)

cfg = (srcs[0] / "config.yaml").read_text()
cfg = re.sub(r"    NETA: .*", f"    NETA: {NETA}", cfg)
cfg = re.sub(r"    ETA_RANGE: .*", f"    ETA_RANGE: [{ETAMIN:g},{ETAMAX:g}]", cfg)
(out / "config.yaml").write_text(cfg)

print(f"merged {len(srcs)} computed slices into {NETA} rapidities in {out}")
print(f"symmetry: worst rel diff = {check_symmetry(out / 'gluon_energy_density.dat'):.3e}")
