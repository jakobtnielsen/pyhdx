# Data Guide — Pipeline #19 PyHDX

## Biological context

**Protein:** SecB — *Escherichia coli* molecular chaperone (UniProt [P0AG86](https://www.uniprot.org/uniprot/P0AG86))

SecB is a cytoplasmic chaperone that prevents premature folding of pre-secretory proteins and
delivers them to the SecA ATPase for translocation through the SecYEG translocon. In solution it
forms a **homotetramer** (dimer of dimers) with a characteristic β-sheet core and a long
hydrophobic groove that binds unfolded client peptides.

Two oligomeric states are measured in this dataset:
- **Tetramer** — wild-type SecB (SecB WT apo), 4 identical subunits
- **Dimer** — engineered mutant (Y109A / T115A / S119A) that locks the protein in a dimeric state

---

## Source data

| Property | Value |
|---|---|
| Dataset ID | `HDX_D9096080` |
| Repository | [HDXMS-database](https://github.com/Jhsmit/HDXMS-database) |
| Publication | Krishnamurthy (2022), *Anal. Chem.* [DOI: 10.1021/acs.analchem.1c02155](https://doi.org/10.1021/acs.analchem.1c02155) |
| Author | Srinath Krishnamurthy (KU Leuven, ORCID: 0000-0001-5492-4450) |
| License | CC0 (public domain) |

### How to download

The dataset is fetched automatically by `pipeline_19_pyhdx.py` via the `hdxms-datasets` package.
To download it manually in Python:

```python
from hdxms_datasets import RemoteDataBase
from pathlib import Path

db = RemoteDataBase(Path("/tmp/hdxms_data"))
ok, msg = db.fetch_dataset("HDX_D9096080")
# Files land in /tmp/hdxms_data/HDX_D9096080/
```

The pipeline caches to `/tmp/hdxms_data/`. Change `DB_DIR` at the top of
`pipeline_19_pyhdx.py` to persist across reboots.

---

## Raw data files

All files live under `HDX_D9096080/` after download.

### `data/ecSecB_apo.csv` — Tetramer experiment (66 kB)

DynamX v3 state-data CSV. Contains **both** the partial-deuteration experiment and the
full-deuteration control, distinguished by the `State` column.

| Column | Description |
|---|---|
| `Protein` | Protein label (always `Accession` here) |
| `Start` | First residue number of peptide (1-indexed, inclusive) |
| `End` | Last residue number of peptide (inclusive) |
| `Sequence` | Amino acid sequence of the peptide |
| `MaxUptake` | Maximum exchangeable amide hydrogens |
| `MHP` | Monoisotopic mass of the undeuterated peptide |
| `State` | `SecB WT apo` (experiment) or `Full deuteration control` |
| `Exposure` | Labelling time in **minutes** (0, 0.167, 0.5, 1, 10, 100) |
| `Center` | Centroid m/z of the isotope envelope |
| `Center SD` | Standard deviation of centroid m/z across replicates |
| `Uptake` | Absolute deuterium uptake in Da |
| `Uptake SD` | Standard deviation of uptake |
| `RT` | Chromatographic retention time (minutes) |
| `RT SD` | Standard deviation of RT |

**States present:**
- `Full deuteration control` — single time point (0.167 min), used for back-exchange normalisation
- `SecB WT apo` — five time points (0.167, 0.5, 1, 10, 100 min)

**Coverage:** 63 unique peptides, 88% sequence coverage (residues 9–155 of 155 total)

### `data/ecSecB_dimer.csv` — Dimer experiment (59 kB)

Same format as above. State label: `SecB his dimer apo`. Five time points.
Mutant sequence carries Y109A / T115A / S119A substitutions vs. WT.

### `data/SecB_structure.pdb` — Crystal structure (379 kB)

PDB entry **1QYN** — crystal structure of *E. coli* SecB at 2.35 Å resolution.
Four chains (A, B, C, D), each covering approximately residues 9–144.
This pipeline uses **chain A** only for structural visualisation.

### `dataset.json` — Machine-readable metadata

Describes states, peptide filters, sequences, mutations, experimental conditions,
and maps data files to protein states. Consumed by `hdxms-datasets.load_dataset()`.

---

## Experimental conditions

| Parameter | Value |
|---|---|
| Labelling buffer pH | 8.0 |
| Temperature | 303.15 K (30 °C) |
| Deuterium percentage | 90% D₂O |
| Labelling times | 0.167, 0.5, 1, 10, 100 min |
| Quench / back-exchange correction | Full-deuteration control (0.167 min) |

> **Note:** Exposure times are stored in **minutes** in the raw CSV files.
> `hdxms-datasets.load_peptides()` converts them to **seconds** automatically
> (multiply by 60), so downstream code works in seconds.

---

## Processed / output data

### `outputs/pyhdx_dG.csv`

Residue-level Gibbs free energy of H/D exchange protection, produced by PyHDX global fitting.

| Column | Description |
|---|---|
| `r_number` | Residue number (1-indexed, matches PDB numbering) |
| `dG_kJ_mol` | ΔG of exchange protection in **kJ/mol** |

- 146 rows (covered residues 10–155)
- ΔG range: 10.1 – 42.9 kJ/mol
- Higher ΔG = more protected = slower exchange = more buried / stable
- Residues not covered by any peptide are absent (no row)

### `outputs/pyhdx_dG_plot.png`

Line + scatter plot of ΔG (kJ/mol) vs. residue number for SecB WT apo.

### `outputs/pyhdx_dG_structure.png`

PyMOL cartoon rendering of SecB chain A (PDB 1QYN) coloured by ΔG:
**blue** (low ΔG, fast exchange) → **white** → **red** (high ΔG, well protected).
Grey residues have no peptide coverage.

---

## Key column naming caveat (for developers)

The raw DynamX CSV uses `End` (inclusive last residue).
`hdxms-datasets.load_peptides()` returns this as `end` (still inclusive).
PyHDX internally uses `stop` as an **exclusive** end (one-past-end, Python convention).
When bridging the two packages you must add:

```python
df["stop"] = df["end"] + 1
```

See `pipeline_19_pyhdx.py` (`to_pandas()` helper) for the full conversion.
