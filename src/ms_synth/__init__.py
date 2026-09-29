"""
ms-synth
========
A known-truth generator of synthetic longitudinal multiple sclerosis (MS)
cohorts, expressed as a continuous-time multi-state Markov model (CT-MSM).

Every simulated cohort carries its complete data-generating process: the
generator matrix, patient-level scaling factors and latent class labels are
known exactly, so methods fitted to the observed (irregular, noisy, censored)
data can be scored against the truth.

Modules
-------
synthetic_data  - generator: Q construction, patient scaling, Gillespie
                  simulation, visit/observation process, validation helpers
msm_fit         - CT-MSM estimation (crude count/time and panel likelihood)
data            - schema validation and transition extraction
utils           - logging, seeding, configuration I/O
graph_utils, percolation, augmentation, visualization
                - optional downstream analysis of G(Q) (percolation
                  application; requires the ``percolation`` extra)
"""

__version__ = "1.0.0"
__author__ = "Nikola Mirkov and contributors"

__all__ = ["__version__"]
