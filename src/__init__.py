"""
ms_percolation_framework
========================
Research framework for analysing Multiple Sclerosis progression via an
abstract percolation approach applied to the directed state-transition
graph of a continuous-time multi-state Markov model (CT-MSM).

Scientific motivation
---------------------
Kannan et al. (2017) [Math. Biosci. 289:1-8] formalised the RRMS→SPMS
transition as a noise-driven threshold crossing in a two-compartment
dynamical system.  Mirkov et al. (2026) [arXiv:2602.08576] reviewed the
percolation / phase-transition hypothesis for MS more broadly, showing
that cognitive semantic networks in MS fragment faster (lower percolation
integrals) and that structural connectomes degrade in a hub-centric,
non-linear fashion consistent with threshold-like behaviour.

This framework bridges those ideas by:
  1. Fitting a CT-MSM to longitudinal EDSS data → infinitesimal generator Q.
  2. Representing Q as a directed weighted graph G(Q).
  3. Sweeping an edge-weight threshold θ, computing directed percolation
     metrics (GSCC size, in/out-component sizes, susceptibility) at each θ.
  4. Identifying the critical threshold θ_c and derived summary statistics
     (proximity δ(t), percolation integrals).
  5. Augmenting the MSM with percolation covariates for trajectory modelling.

Package layout
--------------
src/
  utils.py          – logging, seeding, config I/O
  data.py           – load / preprocess / synthetic generation
  msm_fit.py        – CT-MSM fitting, Q matrix estimation
  graph_utils.py    – build DiGraph from Q, thresholding
  percolation.py    – threshold sweep, GSCC, criticality, integrals
  augmentation.py   – percolation-augmented MSM, trajectory models
  visualization.py  – publication-quality figures
"""

__version__ = "0.1.0"
__author__ = "MS Percolation Framework Contributors"
