"""
Render SecB crystal structure (1QYN) colored by residue-level ΔG from PyHDX.
Color ramp: blue (low ΔG, fast-exchanging) → white → red (high ΔG, well-protected).
Residues with no ΔG data are shown in grey.
Output: outputs/pyhdx_dG_structure.png
"""

import pandas as pd
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

# Load the bundled PDB (1QYN, tetramer A/B/C/D)
PDB = "/tmp/hdxms_data/HDX_D9096080/data/SecB_structure.pdb"
cmd.load(PDB, "secb")

# Cartoon representation, white background
cmd.hide("everything", "secb")
cmd.show("cartoon", "secb")
cmd.bg_color("white")

# Only show chain A (one monomer is sufficient for interpretation)
cmd.remove("secb and not chain A")

# Set all b-factors to -1 (sentinel = no data → will be grey)
cmd.alter("secb", "b = -1.0")

# Assign ΔG values as b-factor; resi selector needs integer string
for resi, dg in dg_kj.dropna().items():
    cmd.alter(f"secb and resi {int(resi)}", f"b = {dg:.4f}")

cmd.rebuild()

# Color by ΔG: grey for no-data, blue→white→red for data range
cmd.color("grey70", "secb")
cmd.spectrum("b", "blue_white_red", "secb and b > -0.5",
             minimum=dg_min, maximum=dg_max)

# Nice view
cmd.orient("secb")
cmd.zoom("secb", 3)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("ray_shadows", 0)
cmd.set("ambient", 0.4)
cmd.set("ray_opaque_background", 1)

# Render and save
OUT = "outputs/pyhdx_dG_structure.png"
cmd.ray(2400, 1600)
cmd.png(OUT, dpi=150)
print(f"Saved: {OUT}  ({Path(OUT).stat().st_size / 1024:.1f} kB)")

cmd.quit()
