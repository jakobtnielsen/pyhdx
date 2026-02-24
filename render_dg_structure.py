"""
Render SecB crystal structure (1QYN) colored by residue-level ΔG from PyHDX.
Shows the full homotetramer (chains A/B/C/D).
Color ramp matches the scatter plot: blue (low ΔG) → white → red (high ΔG).
Residues with no ΔG data are shown in grey.
Output: outputs/pyhdx_dG_structure.png
"""

import pandas as pd
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from pathlib import Path

# ── Load ΔG data ──────────────────────────────────────────────────────────────
dg_kj = pd.read_csv("outputs/pyhdx_dG.csv", index_col="r_number")["dG_kJ_mol"]
dg_min = dg_kj.dropna().min()
dg_max = dg_kj.dropna().max()
print(f"ΔG range: {dg_min:.1f} – {dg_max:.1f} kJ/mol  ({len(dg_kj.dropna())} residues)")

# ── Launch PyMOL in headless mode ─────────────────────────────────────────────
import pymol
from pymol import cmd

pymol.finish_launching(["pymol", "-cq"])   # -c = no GUI, -q = quiet

# Load the full tetramer (chains A/B/C/D)
PDB = "test_data/HDX_D9096080/data/SecB_structure.pdb"
cmd.load(PDB, "secb")

# Cartoon representation, white background
cmd.hide("everything", "secb")
cmd.show("cartoon", "secb")
cmd.bg_color("white")

# Set all b-factors to -1 (sentinel = no data → will be grey)
cmd.alter("secb", "b = -1.0")

# Assign ΔG values to all four chains identically (homotetramer)
for resi, dg in dg_kj.dropna().items():
    cmd.alter(f"secb and resi {int(resi)}", f"b = {dg:.4f}")

cmd.rebuild()

# "hot fire" (inferno) color ramp — matches the scatter plot exactly.
# PyMOL has no native inferno, so we generate colors from matplotlib and assign per residue.
cmap = cm.get_cmap("inferno")
norm = mcolors.Normalize(vmin=dg_min, vmax=dg_max)

cmd.color("grey70", "secb")   # default for uncovered residues

for resi, dg in dg_kj.dropna().items():
    r, g, b_ch, _ = cmap(norm(dg))
    color_name = f"dg_{int(resi)}"
    cmd.set_color(color_name, [r, g, b_ch])
    cmd.color(color_name, f"secb and resi {int(resi)}")

# View: orient to show the flat face of the tetramer β-sheet
cmd.orient("secb")
cmd.zoom("secb", 3)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("ray_shadows", 0)
cmd.set("ambient", 0.45)
cmd.set("ray_opaque_background", 1)

# Render and save
OUT = "outputs/pyhdx_dG_structure.png"
cmd.ray(2400, 1600)
cmd.png(OUT, dpi=150)
print(f"Saved: {OUT}  ({Path(OUT).stat().st_size / 1024:.1f} kB)")

cmd.quit()
