#!/usr/bin/env python3
"""
Concept figure for the PNAS Nexus manuscript (the "clock-free sensitivity
framework"). Model schematic, no empirical data: a prepared metastable state
of readiness decays, and its sensitivity to the decay rate peaks at the 1/e
optimal-stopping point.

  reliability / survival:  R(x)   = e^{-e x}
  rate sensitivity:        S(x)   = x e^{-e x}   (shown normalized to its max)

Both landmark at x = 1/e. This is the standalone version of the schematic that
was previously Panel A of the neural figure; it now sits up front in the main
text (Sec. "A clock-free, state-based reframing"), so the figure title is
supplied by the LaTeX \\caption (no baked-in title here).

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
    "axes.spines.top": False,
    "axes.spines.right": False,
})

E = np.e
INV_E = 1.0 / np.e

C_SURV = "#1f5fa6"   # blue - reliability / survival of the metastable state
C_SENS = "#c0392b"   # red  - rate sensitivity

fig, ax = plt.subplots(figsize=(5.2, 3.3))
fig.subplots_adjust(left=0.11, right=0.965, top=0.95, bottom=0.16)

x = np.linspace(0, 1, 600)
R = np.exp(-E * x)
S = x * np.exp(-E * x)
Sn = S / S.max()

ax.plot(x, R, color=C_SURV, lw=2.3, label=r"reliability  $R(x)=e^{-e x}$")
ax.plot(x, Sn, color=C_SENS, lw=2.3, label=r"rate sensitivity  $S(x)=x\,e^{-e x}$")

# 1/e landmark
ax.axvline(INV_E, color="0.5", ls=":", lw=1.0)
ax.plot([INV_E], [INV_E], "o", color=C_SURV, ms=5.5, zorder=6)
ax.plot([INV_E], [1.0], "o", color=C_SENS, ms=5.5, zorder=6)
ax.text(INV_E + 0.015, 0.04, r"$x = 1/e$", color="0.3", fontsize=8)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1.16)
ax.set_xlabel(r"normalized trial time   $x$")
ax.set_ylabel("amplitude (normalized)")
ax.legend(loc="upper right", fontsize=8, frameon=True, framealpha=0.95,
          borderpad=0.6, handlelength=1.7)

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "results")
fig.savefig(os.path.join(OUT, "FigConcept_Sensitivity.pdf"))
fig.savefig(os.path.join(OUT, "FigConcept_Sensitivity.png"), dpi=300)
print("saved FigConcept_Sensitivity.{pdf,png}")
