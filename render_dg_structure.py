"""
Render two PyMOL structure figures in a single session:
  1. outputs/pyhdx_dG_structure.png  — inferno colormap (ΔG, tetramer)
  2. outputs/pyhdx_ddG_structure.png — RdBu_r diverging (ΔΔG dimer−tetramer)
Residues with no data are grey. Full homotetramer shown in both.
"""

import pandas as pd
import matplotlib.cm as mcm
import matplotlib.colors as mcolors
from pathlib import Path

# ── Load data ─────────────────────────────────────────────────────────────────
dg_kj = pd.read_csv("outputs/pyhdx_dG.csv", index_col="r_number")["dG_kJ_mol"].dropna()
ddg    = pd.read_csv("outputs/pyhdx_ddG.csv", index_col="r_number")["ddG_kJ_mol"].dropna()

dg_min, dg_max = dg_kj.min(), dg_kj.max()
abs_max = max(abs(ddg.min()), abs(ddg.max()))

print(f"ΔG  range: {dg_min:.1f} – {dg_max:.1f} kJ/mol")
print(f"ΔΔG range: {ddg.min():.1f} – {ddg.max():.1f} kJ/mol")

# ── Colormaps ─────────────────────────────────────────────────────────────────
cmap_dg  = mcm.get_cmap("inferno")
norm_dg  = mcolors.Normalize(vmin=dg_min, vmax=dg_max)

cmap_ddg = mcolors.LinearSegmentedColormap.from_list(
    "RdBu5", ["#2166ac", "#d1e5f0", "#f7f7f7", "#fddbc7", "#b2182b"]
)
norm_ddg = mcolors.TwoSlopeNorm(vmin=-abs_max, vcenter=0, vmax=abs_max)

# ── Single PyMOL session ──────────────────────────────────────────────────────
import pymol
from pymol import cmd

pymol.finish_launching(["pymol", "-cq"])

PDB = "test_data/HDX_D9096080/data/SecB_structure.pdb"

def load_and_color(obj_name, data, cmap, norm, out_path):
    """Load PDB, colour by data series, ray-trace and save."""
    cmd.load(PDB, obj_name)
    cmd.hide("everything", obj_name)
    cmd.show("cartoon", obj_name)
    cmd.bg_color("white")
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.45)
    cmd.set("ray_opaque_background", 1)

    # Default grey for uncovered residues
    cmd.color("grey70", obj_name)

    # Per-residue colours from matplotlib colormap
    for resi, val in data.items():
        r, g, b, _ = cmap(norm(val))
        cname = f"{obj_name}_r{int(resi)}"
        cmd.set_color(cname, [r, g, b])
        cmd.color(cname, f"{obj_name} and resi {int(resi)}")

    cmd.orient(obj_name)
    cmd.zoom(obj_name, 3)
    cmd.ray(1200, 800)
    cmd.png(out_path, dpi=150)
    print(f"Saved: {out_path}  ({Path(out_path).stat().st_size / 1024:.1f} kB)")
    cmd.delete(obj_name)

# Figure 1 — ΔG (inferno)
load_and_color("secb_dg", dg_kj, cmap_dg, norm_dg, "outputs/pyhdx_dG_structure.png")

# Figure 2 — ΔΔG (diverging)
load_and_color("secb_ddg", ddg, cmap_ddg, norm_ddg, "outputs/pyhdx_ddG_structure.png")

cmd.quit()
