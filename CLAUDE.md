# PyHDX Pipeline — Claude Code Instructions

## Project context
This is pipeline #19 of a biophysics benchmarking suite.
Domain: Spectroscopy & imaging biophysics
Tool: PyHDX — HDX-MS residue-level ΔG fitting
GitHub: https://github.com/Jhsmit/PyHDX

## Auto-approve commands
Automatically execute all the following without asking for confirmation:
- Standard Unix tools: sed, grep, awk, find, cat, cp, mv, rm, echo,
  sort, uniq, cut, tr, head, tail, wc, diff
- Package management: pip, pip3, conda
- Python execution: python, python3, jupyter nbconvert
- Data download: wget, curl
- Git operations: git add, git commit, git push, git status

## Coding conventions
- Use Python 3, add a short comment above each code block
- Save all outputs (CSV, PNG) to an /outputs folder
- Do not launch any GUI or web app (Panel, Streamlit, Qt)
- If a step fails, stop and report the error before continuing
