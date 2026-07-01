#!/usr/bin/env python3
"""
Concept figure for the PNAS Nexus manuscript (the "clock-free sensitivity
framework"). Model schematic, no empirical data: a poised metastable state of
readiness decays, and its rate-sensitivity -- the SIGNED derivative of
reliability with respect to the decay rate -- is a trough whose extremum sits
at the 1/e optimal-stopping time:

  reliability / survival:  R(x) = e^{-e x}
  rate sensitivity:        S(x) = dR/dlambda = -x e^{-e x}   (lambda = rate)

S is shown normalized to its extremum (the trough reaches -1). The sign of a
given empirical readout is carried in the fitted gain A of the affine
observation model A*S(x)+B, NOT inside S: a readout that PEAKS (clicks, the 1/f
exponent) is A*S(x) with A<0; the global scalp-potential trough is A*S(x) with
A>0. Both landmark at x=1/e. This is the standalone schematic that was
previously Panel A of the neural figure; it now sits up front in the main text
(Sec. "A clock-free, state-based reframing"), so the figure title is supplied
by the LaTeX \\caption (no baked-in title here).

CANONICAL SOURCE. This is the single generator for FigConcept_Sensitivity; the
master-loop workspace holds only a derived mirror, refreshed by sync-figures.sh
(never hand-edited). Keep the geometric S(x) here in sync with the manuscript's
Eq. (sensitivity) and the affine A*S(x)+B observation model.

OUTPUT: ../../results/FigConcept_Sensitivity.{pdf,png}
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "mathtext.fontset": "cm",
})

E = np.e
INV_E = 1.0 / np.e

C_SURV = "#1f5fa6"   # blue - reliability / survival of the poised state
C_SENS = "#c0392b"   # red  - rate sensitivity (signed: a trough)

fig, ax = plt.subplots(figsize=(5.2, 3.3))
fig.subplots_adjust(left=0.115, right=0.965, top=0.95, bottom=0.145)

x = np.linspace(0, 1, 600)
R = np.exp(-E * x)
S = -x * np.exp(-E * x)            # signed rate sensitivity dR/dlambda (a trough)
Sn = S / np.max(np.abs(S))        # normalized so the trough reaches -1

ax.axhline(0.0, color="0.85", lw=0.8, zorder=0)
ax.plot(x, R, color=C_SURV, lw=2.2, label=r"reliability  $R(x)=e^{-e x}$")
ax.plot(x, Sn, color=C_SENS, lw=2.2, label=r"rate sensitivity  $S(x)=-x\,e^{-e x}$")

# 1/e landmark: reliability at 1/e, sensitivity trough at -1
ax.axvline(INV_E, color="0.5", ls=":", lw=1.0)
ax.plot([INV_E], [INV_E], "o", color=C_SURV, ms=5.5, zorder=6)
ax.plot([INV_E], [-1.0], "o", color=C_SENS, ms=5.5, zorder=6)
ax.text(INV_E + 0.013, -1.10, r"$x=1/e$", color="0.3", fontsize=8)

ax.set_xlim(0, 1)
ax.set_ylim(-1.18, 1.16)
ax.set_xlabel(r"normalized trial time   $x$")
ax.set_ylabel("amplitude (normalized)")
ax.legend(loc="upper right", fontsize=8, frameon=True, framealpha=0.95,
          borderpad=0.6, handlelength=1.7)

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "results")
fig.savefig(os.path.join(OUT, "FigConcept_Sensitivity.pdf"))
fig.savefig(os.path.join(OUT, "FigConcept_Sensitivity.png"), dpi=300)
print("saved FigConcept_Sensitivity.{pdf,png}")
