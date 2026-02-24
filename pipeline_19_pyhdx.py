"""
Pipeline #19 — PyHDX HDX-MS residue-level ΔG fitting
Dataset: SecB WT apo (HDX_D9096080) from HDXMS-database
Outputs: outputs/pyhdx_dG.csv, outputs/pyhdx_dG_plot.png
"""

# ── Standard and third-party imports ──────────────────────────────────────────
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no GUI
import matplotlib.pyplot as plt
import pandas as pd
import polars as pl

import pyhdx
from hdxms_datasets import RemoteDataBase, load_peptides
from pyhdx import HDXMeasurement
from pyhdx.fitting import fit_rates_half_time_interpolate, fit_gibbs_global
from pyhdx.process import apply_control, correct_d_uptake

# ── Suppress deprecation noise from pyhdx internals ──────────────────────────
warnings.filterwarnings("ignore", category=DeprecationWarning)

# ── Output directory ──────────────────────────────────────────────────────────
OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

# ===========================================================================
# Step 1 — Confirm installation
# ===========================================================================
print(f"\n[Step 1] PyHDX version: {pyhdx.__version__}")


# ===========================================================================
# Step 2 — Download / load the SecB dataset from the HDXMS remote database
# ===========================================================================
print("\n[Step 2] Loading SecB dataset (HDX_D9096080) …")

DB_DIR = Path("test_data")

db = RemoteDataBase(DB_DIR)
print(f"  Using local dataset from {DB_DIR.resolve()}")

dataset = db.load_dataset("HDX_D9096080")
print(f"  Dataset: {dataset.description[:60]}")
print(f"  States: {[s.name for s in dataset.states]}")

# Use 'Tetramer' state (SecB WT apo, partially-deuterated experiment)
state = dataset.states[0]
print(f"  Selected state: {state.name!r}")
print(f"  Peptide sets: {len(state.peptides)}")
print(f"    [0] {state.peptides[0].deuteration_type.name}, "
      f"exposures={state.peptides[0].filters.get('Exposure')}")
print(f"    [1] {state.peptides[1].deuteration_type.name}, "
      f"exposures={state.peptides[1].filters.get('Exposure')}")


# ===========================================================================
# Step 3 — Load and pre-process peptide data
# ===========================================================================
print("\n[Step 3] Pre-processing peptide data …")

base_dir = Path(".")  # data_file paths in dataset.json are already relative to cwd

# Load narwhals DataFrames and convert to pandas
# pyhdx uses 'stop' as an exclusive (one-past-end) residue number,
# while hdxms-datasets returns 'end' as the last inclusive residue.
def to_pandas(nw_df) -> pd.DataFrame:
    """Convert narwhals → polars → pandas and add exclusive 'stop' column."""
    df = nw_df.to_native().to_pandas()
    df["stop"] = df["end"] + 1  # pyhdx requires exclusive stop
    return df

apo_pd = to_pandas(load_peptides(state.peptides[0], base_dir=base_dir))
fd_pd  = to_pandas(load_peptides(state.peptides[1], base_dir=base_dir))

print(f"  Partially-deuterated peptides: {apo_pd.shape[0]} rows, "
      f"{apo_pd['exposure'].nunique()} timepoints")
print(f"  FD control peptides:           {fd_pd.shape[0]} rows")
print(f"  Timepoints (s): {sorted(apo_pd['exposure'].unique())}")

# Back-exchange normalization (computes 'rfu' and 'rfu_sd')
peptides_corr = apply_control(apo_pd, fd_pd)

# Drop first residue (N-terminal back-exchange correction), compute 'uptake_corrected'
peptides_final = correct_d_uptake(peptides_corr)

print(f"  Columns after processing: {list(peptides_final.columns)}")


# ===========================================================================
# Step 4 — Create HDXMeasurement object
# ===========================================================================
print("\n[Step 4] Creating HDXMeasurement object …")

SECB_SEQ = state.protein_state.sequence

hdxm = HDXMeasurement(
    peptides_final,
    temperature=303.15,   # 30 °C in Kelvin
    pH=8.0,
    sequence=SECB_SEQ,
)
print(hdxm)


# ===========================================================================
# Step 5 — Fit exchange rates (half-time interpolation)
# ===========================================================================
print("\n[Step 5] Fitting exchange rates via half-time interpolation …")

rate_result = fit_rates_half_time_interpolate(hdxm)
rate_series  = rate_result.output["rate"]   # pd.Series, index = r_number
print(f"  Rate result shape: {rate_result.output.shape}")
print(f"  Rate range: {rate_series.min():.4g} – {rate_series.max():.4g} s⁻¹")


# ===========================================================================
# Step 6 — Convert rates to ΔG initial guesses
# ===========================================================================
print("\n[Step 6] Computing ΔG initial guesses …")

gibbs_guess = hdxm.guess_deltaG(rate_series)   # pd.Series in J/mol
print(f"  Guess shape: {gibbs_guess.shape}, NaN count: {gibbs_guess.isna().sum()}")


# ===========================================================================
# Step 7 — Global PyTorch regularized ΔG fit
# ===========================================================================
print("\n[Step 7] Running global PyTorch ΔG fit …")
print("  (epochs=100000, r1=0.1, lr=1e5, stop_loss=1e-5, patience=50)")

gibbs_result = fit_gibbs_global(
    hdxm,
    gibbs_guess,
    epochs=100_000,
    r1=0.1,
    stop_loss=1e-5,
    patience=50,
    lr=1e5,
)

print(f"  Fit complete. Epochs run: {gibbs_result.metadata.get('epochs_run')}")
print(f"  Total loss:   {gibbs_result.metadata.get('total_loss'):.6g}")
print(f"  MSE loss:     {gibbs_result.metadata.get('mse_loss'):.6g}")


# ===========================================================================
# Step 8 — Save ΔG results to CSV
# ===========================================================================
print("\n[Step 8] Saving ΔG results …")

# dG DataFrame: index = r_number, columns = [state_name]
dg_df = gibbs_result.dG.copy()
dg_df.index.name = "r_number"
dg_df.columns.name = "state"

# Convert J/mol → kJ/mol
dg_kj = dg_df / 1000.0
dg_kj.columns = [f"dG_kJ_mol"]

csv_path = OUTPUT_DIR / "pyhdx_dG.csv"
dg_kj.to_csv(csv_path)
print(f"  Saved: {csv_path}  ({len(dg_kj)} rows)")
print(dg_kj.head(5))


# ===========================================================================
# Step 9 — Save residue-level ΔG plot
# ===========================================================================
print("\n[Step 9] Saving ΔG plot …")

dg_values = dg_kj["dG_kJ_mol"].dropna()

# Use the same blue→white→red colormap as the PyMOL structure figure
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# "hot fire" colormap: black/purple (low ΔG) → orange → yellow (high ΔG)
cmap = cm.get_cmap("inferno")
norm = mcolors.Normalize(vmin=dg_values.min(), vmax=dg_values.max())
colors = cmap(norm(dg_values.values))

fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(dg_values.index, dg_values.values, color="lightgray", linewidth=0.8, zorder=1)
ax.scatter(dg_values.index, dg_values.values, s=25, c=colors, zorder=2)
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
plt.colorbar(sm, ax=ax, label="ΔG (kJ/mol)", fraction=0.03, pad=0.02)
ax.set_xlabel("Residue number")
ax.set_ylabel("ΔG (kJ/mol)")
ax.set_title("PyHDX — Residue-level ΔG (SecB WT apo)")
ax.axhline(0, color="gray", linewidth=0.6, linestyle="--")
plt.tight_layout()

png_path = OUTPUT_DIR / "pyhdx_dG_plot.png"
fig.savefig(png_path, dpi=150)
plt.close()
print(f"  Saved: {png_path}  ({png_path.stat().st_size / 1024:.1f} kB)")


# ===========================================================================
# Step 10 — Fit dimer state ΔG
# (Dimer has no FD control; reuse tetramer FD control — standard practice)
# ===========================================================================
print("\n[Step 10] Fitting dimer state ΔG …")

state_dimer = dataset.states[1]   # SecB his dimer apo (Y109A/T115A/S119A)
dim_pl = load_peptides(state_dimer.peptides[0], base_dir=base_dir).to_native()
dim_pd = dim_pl.to_pandas()
dim_pd["stop"] = dim_pd["end"] + 1

# Reuse tetramer FD control for back-exchange normalisation
dim_corr  = apply_control(dim_pd, fd_pd)
dim_final = correct_d_uptake(dim_corr)

DIMER_SEQ = state_dimer.protein_state.sequence
hdxm_dim = HDXMeasurement(dim_final, temperature=303.15, pH=8.0, sequence=DIMER_SEQ)
print(hdxm_dim)

rate_dim        = fit_rates_half_time_interpolate(hdxm_dim)
gibbs_guess_dim = hdxm_dim.guess_deltaG(rate_dim.output["rate"])
gibbs_dim       = fit_gibbs_global(hdxm_dim, gibbs_guess_dim,
                                   epochs=100_000, r1=0.1,
                                   stop_loss=1e-5, patience=50, lr=1e5)

dg_dim_kj = gibbs_dim.dG.copy() / 1000.0
dg_dim_kj.columns = ["dG_kJ_mol"]
dg_dim_kj.index.name = "r_number"

csv_dim_path = OUTPUT_DIR / "pyhdx_dG_dimer.csv"
dg_dim_kj.to_csv(csv_dim_path)
print(f"  Saved: {csv_dim_path}  ({len(dg_dim_kj)} rows)")


# ===========================================================================
# Step 11 — Compute ΔΔG = ΔG_dimer − ΔG_tetramer and save CSV
# ===========================================================================
print("\n[Step 11] Computing ΔΔG (dimer − tetramer) …")

# Align on common residues
common = dg_kj.index.intersection(dg_dim_kj.index)
ddg = (dg_dim_kj.loc[common, "dG_kJ_mol"] - dg_kj.loc[common, "dG_kJ_mol"]).dropna()
ddg.name = "ddG_kJ_mol"
ddg.index.name = "r_number"

csv_ddg_path = OUTPUT_DIR / "pyhdx_ddG.csv"
ddg.to_csv(csv_ddg_path, header=True)
print(f"  ΔΔG range: {ddg.min():.1f} – {ddg.max():.1f} kJ/mol  ({len(ddg)} residues)")
print(f"  Saved: {csv_ddg_path}")


# ===========================================================================
# Step 12 — ΔΔG scatter plot (diverging colormap, 0-centred)
# ===========================================================================
print("\n[Step 12] Saving ΔΔG plot …")

# Symmetric colour scale centred at zero
abs_max = max(abs(ddg.min()), abs(ddg.max()))
norm_ddg = mcolors.TwoSlopeNorm(vmin=-abs_max, vcenter=0, vmax=abs_max)
cmap_ddg = cm.get_cmap("RdBu_r")   # red = dimer more protected, blue = tetramer more protected
colors_ddg = cmap_ddg(norm_ddg(ddg.values))

fig, ax = plt.subplots(figsize=(12, 4))
ax.axhline(0, color="gray", linewidth=0.8, linestyle="--", zorder=1)
ax.plot(ddg.index, ddg.values, color="lightgray", linewidth=0.8, zorder=2)
ax.scatter(ddg.index, ddg.values, s=25, c=colors_ddg, zorder=3)
sm2 = cm.ScalarMappable(cmap=cmap_ddg, norm=norm_ddg)
sm2.set_array([])
plt.colorbar(sm2, ax=ax, label="ΔΔG (kJ/mol)", fraction=0.03, pad=0.02)
ax.set_xlabel("Residue number")
ax.set_ylabel("ΔΔG (kJ/mol)")
ax.set_title("PyHDX — ΔΔG dimer − tetramer (SecB)  |  red = dimer more protected")
plt.tight_layout()

png_ddg_path = OUTPUT_DIR / "pyhdx_ddG_plot.png"
fig.savefig(png_ddg_path, dpi=150)
plt.close()
print(f"  Saved: {png_ddg_path}  ({png_ddg_path.stat().st_size / 1024:.1f} kB)")


# ===========================================================================
# Verification summary
# ===========================================================================
print("\n=== Verification ===")
print(f"  pyhdx installed:            {pyhdx.__version__}")
print(f"  Tetramer peptides/residues: {hdxm.Np} / {hdxm.Nr}")
print(f"  Dimer    peptides/residues: {hdxm_dim.Np} / {hdxm_dim.Nr}")
print(f"  Tetramer CSV rows:          {len(dg_kj)}")
print(f"  Dimer    CSV rows:          {len(dg_dim_kj)}")
print(f"  ΔΔG CSV rows:               {len(ddg)}")
print(f"  outputs/pyhdx_dG.csv:           {csv_path.exists()}")
print(f"  outputs/pyhdx_dG_dimer.csv:     {csv_dim_path.exists()}")
print(f"  outputs/pyhdx_ddG.csv:          {csv_ddg_path.exists()}")
print(f"  outputs/pyhdx_dG_plot.png:      {png_path.exists()}")
print(f"  outputs/pyhdx_ddG_plot.png:     {png_ddg_path.exists()}")
print("\nDone.")
