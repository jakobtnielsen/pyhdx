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

DB_DIR = Path("/tmp/hdxms_data")
DB_DIR.mkdir(exist_ok=True)

db = RemoteDataBase(DB_DIR)

# Fetch from remote only if not already cached
if "HDX_D9096080" not in db.local_datasets:
    ok, msg = db.fetch_dataset("HDX_D9096080")
    print(f"  Fetch status: {ok}  {msg}")
else:
    print("  Using cached dataset.")

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

base_dir = db.database_dir / "HDX_D9096080"

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

fig, ax = plt.subplots(figsize=(12, 4))
ax.scatter(dg_values.index, dg_values.values, s=12, color="steelblue", alpha=0.85)
ax.plot(dg_values.index, dg_values.values)
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
# Verification summary
# ===========================================================================
print("\n=== Verification ===")
print(f"  pyhdx installed:         {pyhdx.__version__}")
print(f"  HDXMeasurement peptides: {hdxm.Np}")
print(f"  HDXMeasurement residues: {hdxm.Nr}")
print(f"  CSV rows:                {len(dg_kj)}")
print(f"  CSV non-NaN ΔG values:   {dg_values.shape[0]}")
print(f"  CSV path exists:         {csv_path.exists()}")
print(f"  PNG path exists:         {png_path.exists()}")
print(f"  PNG size (bytes):        {png_path.stat().st_size}")
print("\nDone.")
