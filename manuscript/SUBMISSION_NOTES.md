# Submission notes — ms-synth (GigaScience Technical Note)

## Target and fit
**GigaScience, Technical Note.** Chosen because its requirements match this work:
open-source software with reproducible examples (submissions that reviewers cannot
reproduce are rejected without review), a structured abstract (Background /
Findings / Conclusions, ≤ 250 words), references numbered in order of citation,
software cited with version and persistent identifier, and data deposition in GigaDB.
Fallback: *PLOS Computational Biology* (Software) — the sections map directly
(Background → Introduction; Findings → Design and Implementation + Results).

| Requirement | Status |
|---|---|
| Structured abstract ≤ 250 words, no citations, URL at end | Done (230 words) |
| 3–10 keywords | Done (8) |
| Numbered references in order of citation, with DOIs | Done (50 refs, 47 DOIs; 3 have none by nature) |
| "Availability of supporting source code and requirements" block | Done — RRID, bio.tools ID and Zenodo DOI are placeholders |
| Data availability, abbreviations, declarations | Done — some items need author input (below) |
| Defined code version, test data, expected outputs | Done: tag `v1.0.0`, `reference/`, `verify_reproduction.py` |
| Figures at end of TeX file for submission | `make submission` (endfloat) |
| OUP LaTeX template | **Author action**: port text into the OUP template at submission |

## Author actions before submission (in order)
1. Review the commits in the bundle, then push commits **and** the `v1.0.0` tag to
   `github.com/nikola-m/ms-synth`.
2. Confirm the GitHub Actions CI run passes (configured, not yet run on GitHub —
   macOS/Windows jobs untested). Build the Dockerfile once (untested here).
3. Enable the Zenodo–GitHub integration, create the GitHub release `v1.0.0`, and copy
   the minted DOI into: `CITATION.cff`, `manuscript/references.bib` (`mssynth`), and
   the Availability section of the manuscript.
4. Register the software at SciCrunch (RRID) and bio.tools; fill in both IDs.
5. Complete: corresponding-author e-mail and all author e-mails (title page); CRediT
   roles per author; funding confirmation; **generative-AI disclosure** in
   Acknowledgements (GigaScience requires tool, model version, dates, user).
6. Run `python manuscript/check_references.py` (needs internet) — resolves every DOI
   via Crossref and flags any title mismatch.
7. Port into the OUP template; submit `ms_synth_manuscript_submission.pdf` as the
   "Reference PDF", with `references.bib`, figures and Additional file 1.
8. On acceptance: deposit the code snapshot and reference outputs in GigaDB.

## Substantive corrections relative to the previous report
- **Q-recovery interpretation.** The previous claim that the 0.31 relative error
  reflected finite-sample variance was wrong. With known truth it decomposes into a
  heterogeneity floor (oracle error 0.25) plus a −18% panel-observation bias of the
  crude estimator. This is now a headline result (estimator benchmark, Fig. 5).
- **Panel MLE claim.** The report stated N_ij/T_i is the MLE for panel data; it is the
  MLE only for continuously observed paths. Corrected.
- **Provenance.** Removed the abstract's residual claim of consistency with Danish
  registry/MSBase data: the generator matrix is hand-calibrated to aggregate targets.
- **SPMS-phase relapse benchmark.** The "0.23–0.41" range had no source; Big MS Data
  reports mean 0.23 (SD 0.34). The synthetic 0.281 is slightly above; Fig. 3D corrected.
- **Misattributed citations.** Vasanthaprasad 2022 (SPMS *prevalence*) was cited for
  20-year conversion rates — removed. The Danish registry citation ("Mult Scler
  28:706–719") was invalid — replaced by Magyari et al. 2021. Two grey-literature
  items (a 2019 newsletter, an MS Coalition PDF) replaced by peer-reviewed sources.
- **Milestone definitions.** Stated explicitly that synthetic milestones are
  first-reached and unconfirmed, and that 13.6% of patients have EDSS ≥ 3 at the first
  visit — this partly explains the shorter EDSS-3 median (8.4 vs 10 years).
- **Unverifiable "published ranges"** in the descriptive table (e.g. "70–85%") removed;
  only values traceable to a named cohort remain.

## Reference verification log (50 entries)
- **Web-verified this session** (publisher, PubMed/PMC, DOAJ or indexed record):
  magyari2021, butzkueven2006, vukusic2020, trojano2019, hillert2015, guo2026,
  kappos2020pira, ribbons2015, tilling2016, brown2019, mandel2007, lange2015,
  lublin2022, stadler2022 (no DOI; USENIX URL), signori2023.
- **Verified from the project PDFs:** mirkov2026, kannan2017.
- **Verified in the earlier session** (verified-data compilation / web lookup):
  scalfari2010, kappos2018, weinshenker1989, tremlett2006, koch2010, palace2014,
  leray2010, confavreux2006, polman2006, samjoo2023, mcginley2021.
- **Cited with DOI on the GigaScience author page:** smith2016.
- **Standard references, DOIs from canonical records, not individually re-checked
  this session:** walton2020, lublin2014, kurtzke1983, nowok2016, walonoski2018,
  burton2006, morris2019, kalbfleisch1985, jackson2011, gillespie1977, newman2010,
  tutuncu2013, harris2020, virtanen2020, mckinney2010, davidsonpilon2019, hunter2007,
  sandve2013, wilkinson2016, gdpr2016 (URL), mssynth (own software; DOI pending).
  → covered by `check_references.py` (step 6 above).
- Co-author lists are truncated ("et al.") to authors confirmed in verified records.
