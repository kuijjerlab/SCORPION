# SERGIO ground-truth simulation for `inferEdgeSign()` tests

Builds `tests/testthat/fixtures/sergio_signed.rds`: single-cell counts simulated with
[SERGIO](https://github.com/PayamDiba/SERGIO) (steady state, commit `a6190b7`) from a
signed GRN sampled from the `mor` column of `scorpionTest$tf`.

- 12 TFs (8 master regulators, 4 TFs regulated by them) and 287 targets; 575 edges
  (467 activating, 108 repressing); at most 3 regulators per target.
- 10 cell states x 150 cells: states S0-S4 form cell type A, S5-S9 cell type B.
- Technical noise: outlier genes, library size, dropout, UMI sampling, plus a
  gene-specific multiplicative effect for each of 4 donors.

Run from the package root:

```sh
uv venv sergio-venv --python 3.11
VIRTUAL_ENV=sergio-venv uv pip install numpy==1.23.5 scipy==1.10.1 networkx==3.1
git clone https://github.com/PayamDiba/SERGIO sergio-src && git -C sergio-src checkout a6190b7

Rscript data-raw/sergio/01_build_grn.R data-raw/sergio/out
sergio-venv/bin/python data-raw/sergio/02_simulate.py sergio-src data-raw/sergio/out
Rscript data-raw/sergio/03_make_fixture.R data-raw/sergio/out
```

SERGIO uses `np.int`/`np.float`, so numpy must be older than 1.24.
