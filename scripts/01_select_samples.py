"""Apply the QC exclusions and write the list of samples per MST subset to output/mst/samples_used.tsv."""
from __future__ import annotations

import argparse
from pathlib import Path

from common import OUTPUT, SUBSETS, alignment_path, load_samples
from kansasii_paths import ASSEMBLIES, SCREENING_MAP_RESULTS

ap = argparse.ArgumentParser()
ap.add_argument("--screening-map-results", type=Path, default=SCREENING_MAP_RESULTS)
ap.add_argument("--assemblies", type=Path, default=ASSEMBLIES, help="<sample>/5_typing/kansasii_snippy/ per sample")
args = ap.parse_args()

OUTPUT.mkdir(parents=True, exist_ok=True)
s = load_samples(args.screening_map_results)
missing = [n for n in s[s.keep]["sample"] if not alignment_path(n, args.assemblies).exists()]
assert not missing, f"no snippy alignment for {missing}"
s["subsets"] = [",".join(k for k, sp in SUBSETS.items() if r.keep and (sp is None or r.species in sp))
                for r in s.itertuples()]
s.to_csv(OUTPUT / "samples_used.tsv", sep="\t", index=False)
print("excluded:\n", s[~s.keep][["sample", "species", "reason"]].to_string(index=False))
for k in SUBSETS:
    print(f"{k}: {s.subsets.str.split(',').apply(lambda x, k=k: k in x).sum()} samples")
