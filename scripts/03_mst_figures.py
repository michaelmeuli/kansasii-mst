"""Minimum spanning tree figures from the core-SNP distance matrices (one per subset).
Genomes with identical core SNPs (distance 0) are collapsed, and MST edges of at most N SNPs are contracted (single
linkage, N per subset in common.MERGE_SNPS or --merge N), so each node is a group of genomes whose area scales with the
genome count. The remaining edges are original MST edges; edge labels are SNP differences. --merge 0 = identical only."""
from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.lines import Line2D
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree

from common import MERGE_SNPS, OUT, SPECIES_COLORS, SUBSETS

TITLES = {"complex": "M. kansasii complex", "kansasii": "M. kansasii (sensu stricto)",
          "pseudokansasii": "M. pseudokansasii", "persicum": "M. persicum"}
ap = argparse.ArgumentParser()
ap.add_argument("--merge", type=int, default=None, help="merge nodes joined by MST edges of <= N SNPs (default: per subset)")
ap.add_argument("--subset", choices=list(SUBSETS), default=None, help="only this subset (default: all)")
args = ap.parse_args()
samples = pd.read_csv(OUT / "samples_used.tsv", sep="\t").set_index("sample")

for name in [args.subset] if args.subset else SUBSETS:
    D = pd.read_csv(OUT / f"core_snp_distance_{name}.csv", index_col=0)
    # collapse identical haplotypes (rows of the distance matrix equal <=> distance 0)
    hap: list[list[str]] = []
    for n in D.index:
        for g in hap:
            if D.loc[n, g[0]] == 0:
                g.append(n)
                break
        else:
            hap.append([n])
    reps = [g[0] for g in hap]
    H = D.loc[reps, reps].to_numpy(dtype=float)
    if len(reps) == 1:
        print(f"{name}: single haplotype, skipped")
        continue
    # zero entries would be dropped by the sparse MST; distinct haplotypes differ by >= 1 so this is safe
    T = minimum_spanning_tree(csr_matrix(H)).tocoo()
    merge = MERGE_SNPS[name] if args.merge is None else args.merge
    short = T.data <= merge
    _, comp = connected_components(csr_matrix((np.ones(int(short.sum())), (T.row[short], T.col[short])), shape=H.shape),
                                   directed=False)
    groups: list[list[str]] = [[] for _ in range(int(comp.max()) + 1)]
    for h, c in zip(hap, comp.tolist()):
        groups[c].extend(h)
    k: int = len(groups)
    G = nx.Graph()
    mixed = 0
    for i, g in enumerate(groups):
        sp = samples.loc[g, "species"].value_counts()
        mixed += len(sp) > 1
        G.add_node(i, n=len(g), species=sp.index[0])
    if mixed:
        print(f"WARNING {name}: {mixed} merged nodes mix species (coloured by majority)")
    for i, j, w in zip(T.row[~short].tolist(), T.col[~short].tolist(), T.data[~short].tolist()):
        G.add_edge(int(comp[i]), int(comp[j]), weight=w)
    if k == 1:
        print(f"{name}: everything merged at <= {merge} SNPs, skipped")
        continue
    assert G.number_of_edges() == k - 1 and nx.is_connected(G)
    # layout: tree path lengths, compressed (sqrt) so long and short edges are both visible
    Gs = nx.Graph((u, v, {"weight": np.sqrt(d["weight"])}) for u, v, d in G.edges(data=True))
    dist = dict(nx.all_pairs_dijkstra_path_length(Gs))
    pos = nx.kamada_kawai_layout(Gs, dist=dist)

    big = k > 100
    fig, ax = plt.subplots(figsize=(24, 20) if big else (11, 9) if k <= 30 else (16, 13))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#888888", width=1.0)
    species_nodes = sorted({d["species"] for _, d in G.nodes(data=True)})
    nodes = nx.draw_networkx_nodes(G, pos, ax=ax, node_color=[SPECIES_COLORS.get(G.nodes[v]["species"], "#999999") for v in G],
                           node_size=[min(60 + 90 * G.nodes[v]["n"], 4000) for v in G], edgecolors="black", linewidths=0.6)
    def node_label(g: list[str]) -> str:
        return "\n".join(g[:3]) + (f"\n+{len(g) - 3} more" if len(g) > 3 else "") + (f" (n={len(g)})" if len(g) > 1 else "")
    pd.DataFrame({"node": range(k), "n": [len(g) for g in groups],
                  "species": [G.nodes[i]["species"] for i in range(k)],
                  "max_internal_snps": [int(D.loc[g, g].to_numpy().max()) for g in groups],
                  "samples": [",".join(g) for g in groups]}).to_csv(OUT / f"mst_nodes_{name}.tsv", sep="\t", index=False)
    thr = 0 if k <= 60 else float(np.percentile([d["weight"] for *_, d in G.edges(data=True)], 75))
    edge_txt = nx.draw_networkx_edge_labels(G, pos, ax=ax, font_size=5, rotate=False, bbox={"alpha": 0, "pad": 0},
                                 edge_labels={(u, v): f"{int(d['weight'])}" for u, v, d in G.edges(data=True)
                                              if d["weight"] > thr})
    # node labels start just outside their node (pointing away from their tree neighbours), then are pushed apart
    ax.axis("off")
    fig.tight_layout()
    fig.canvas.draw()  # transforms are only final after a draw
    texts = []
    for i, g in enumerate(groups):
        r_px = np.sqrt(min(60 + 90 * G.nodes[i]["n"], 4000) / np.pi) * fig.dpi / 72
        cpx = ax.transData.transform(pos[i])
        # direction = middle of the widest angular gap between the node's tree neighbours (leaf: away from its parent)
        ang = np.sort([np.arctan2(*(ax.transData.transform(pos[j]) - cpx)[::-1]) for j in G[i]])
        gaps = np.diff(np.append(ang, ang[0] + 2 * np.pi))
        a = ang[int(np.argmax(gaps))] + gaps.max() / 2
        vec = np.array([np.cos(a), np.sin(a)])
        x, y = ax.transData.inverted().transform(cpx + vec * (r_px + 10 + 12 * max(G.degree(i) - 2, 0)))
        texts.append(ax.text(x, y, node_label(g), fontsize=5 if big else 6 if k > 30 else 7, ha="center", va="center",
                             bbox={"facecolor": "white", "alpha": 0.7, "pad": 0.3, "edgecolor": "none"}))
    adjust_text(texts, x=[pos[i][0] for i in range(k)], y=[pos[i][1] for i in range(k)], objects=list(edge_txt.values()),
                ax=ax, expand=(1.1, 1.3), force_text=(0.3, 0.6), force_static=(0.5, 0.8), iter_lim=1000,
                arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.4})
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc=SPECIES_COLORS.get(sp, "#999999"), mec="black", label=sp)
                       for sp in species_nodes], title="species", loc="best", frameon=False)
    nvar = len(pd.read_csv(OUT / f"core_snp_sites_{name}.csv"))
    ax.set_title(f"{TITLES[name]}: MST of core-genome SNPs\n{len(D)} genomes, {k} nodes (merged at <= {merge} SNPs), {nvar} variable core sites "
                 f"(edge labels = SNP differences)", fontsize=11)
    ax.axis("off")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"mst_{name}.{ext}", dpi=200)
    plt.close(fig)
    print(f"{name}: wrote mst_{name}.png/pdf ({len(D)} genomes, {len(hap)} haplotypes, {k} nodes at <= {merge} SNPs)")
