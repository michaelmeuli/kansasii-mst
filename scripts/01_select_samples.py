"""Apply the QC exclusions and write the list of samples per MST subset to output/mst/samples_used.tsv."""
from __future__ import annotations

from common import OUT, SUBSETS, alignment_path, load_samples

OUT.mkdir(parents=True, exist_ok=True)
s = load_samples()
missing = [n for n in s[s.keep]["sample"] if not alignment_path(n).exists()]
assert not missing, f"no snippy alignment for {missing}"
s["subsets"] = [",".join(k for k, sp in SUBSETS.items() if r.keep and (sp is None or r.species in sp))
                for r in s.itertuples()]
s.to_csv(OUT / "samples_used.tsv", sep="\t", index=False)
print("excluded:\n", s[~s.keep][["sample", "species", "reason"]].to_string(index=False))
for k in SUBSETS:
    print(f"{k}: {s.subsets.str.split(',').apply(lambda x, k=k: k in x).sum()} samples")
