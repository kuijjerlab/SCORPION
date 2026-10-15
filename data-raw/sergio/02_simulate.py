"""Simulate single-cell counts with SERGIO (steady state) from the GRN built by 01_build_grn.R.

Usage: python 02_simulate.py <sergioRepoDir> <outDir>
Requires SERGIO (https://github.com/PayamDiba/SERGIO, commit a6190b7) with
numpy==1.23.5, scipy==1.10.1, networkx==3.1 (SERGIO uses np.int / np.float).
Writes: counts.csv (genes x cells, UMI counts) and cells.csv (cell metadata).
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, sys.argv[1])
from SERGIO.sergio import sergio  # noqa: E402

outDir = sys.argv[2]
np.random.seed(20261014)

nBins = 10  # bins 0-4: cell type A, bins 5-9: cell type B
nCells = 150  # cells per bin
nDonors = 4

with open(os.path.join(outDir, "genes.csv")) as f:
    genes = [row["gene"] for row in csv.DictReader(f)]
nGenes = len(genes)

sim = sergio(number_genes=nGenes, number_bins=nBins, number_sc=nCells,
             noise_params=1, decays=0.8, sampling_state=15, noise_type="dpd")
sim.build_graph(input_file_taregts=os.path.join(outDir, "targets.csv"),
                input_file_regs=os.path.join(outDir, "regs.csv"), shared_coop_state=2)
sim.simulate()
expr = sim.getExpressions()  # bins x genes x cells

# Technical noise (SERGIO steady-state recipe): outlier genes, library size, dropout
expr = sim.outlier_effect(expr, outlier_prob=0.01, mean=0.8, scale=1)
_, expr = sim.lib_size_effect(expr, mean=4.6, scale=0.4)
expr = np.multiply(sim.dropout_indicator(expr, shape=8, percentile=45), expr)
expr = np.concatenate(expr, axis=1)  # genes x (bins * cells)

# Donor effects: random donor per cell and a gene-specific multiplicative donor shift
donor = np.random.randint(nDonors, size=expr.shape[1])
donorFactor = np.random.lognormal(mean=0, sigma=0.25, size=(nGenes, nDonors))
expr = expr * donorFactor[:, donor]
counts = sim.convert_to_UMIcounts(expr)

cellBin = np.repeat(np.arange(nBins), nCells)
cellIds = ["cell%05d" % i for i in range(expr.shape[1])]
with open(os.path.join(outDir, "counts.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gene"] + cellIds)
    for g, row in zip(genes, counts):
        w.writerow([g] + row.astype(int).tolist())
with open(os.path.join(outDir, "cells.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["cell_id", "state", "cell_type", "donor"])
    for c, b, d in zip(cellIds, cellBin, donor):
        w.writerow([c, "S%d" % b, "A" if b < 5 else "B", "D%d" % (d + 1)])

print("Simulated %d genes x %d cells; zero fraction %.2f" % (nGenes, counts.shape[1], np.mean(counts == 0)))
