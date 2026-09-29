#!/usr/bin/env python3
"""Resolve every DOI in references.bib via the Crossref API and compare titles.

Standard library only. Run from this directory:  python check_references.py
Prints OK / MISMATCH / NOT FOUND per entry; exit status 1 if any problem.
(Written for the submission check; requires internet access.)
"""
import difflib, json, re, sys, time, urllib.error, urllib.parse, urllib.request

bib = open("references.bib", encoding="utf-8").read()
entries = re.findall(r"@\w+\{([^,]+),(.*?)\}\s*(?=@|\Z)", bib, re.S)
bad = 0
for key, body in entries:
    doi = re.search(r"doi=\{([^}]*)\}", body)
    title = re.search(r"title=\{(.*?)\},(?:journal|booktitle|howpublished|publisher|year)", body, re.S)
    if not doi:
        print(f"-- {key}: no DOI (URL reference)"); continue
    t_bib = re.sub(r"[{}\\]", "", title.group(1)).lower() if title else ""
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi.group(1))
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "ms-synth-refcheck"}), timeout=20) as r:
            msg = json.load(r)["message"]
        t_cr = " ".join(msg.get("title", [""])).lower()
        sim = difflib.SequenceMatcher(None, t_bib, t_cr).ratio()
        status = "OK" if sim > 0.85 else "MISMATCH"
        bad += status != "OK"
        print(f"{status:8s} {key}: similarity {sim:.2f}" + ("" if status == "OK" else f"\n         bib: {t_bib}\n         crossref: {t_cr}"))
    except urllib.error.HTTPError as e:
        # arXiv DOIs (10.48550) are registered with DataCite, not Crossref
        if "10.48550" in doi.group(1):
            print(f"-- {key}: DataCite DOI (check https://doi.org/{doi.group(1)})")
        else:
            bad += 1; print(f"NOT FOUND {key}: {doi.group(1)} (HTTP {e.code})")
    time.sleep(0.2)
sys.exit(1 if bad else 0)
