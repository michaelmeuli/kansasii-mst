"""Minimum spanning tree figures from the core-SNP distance matrices (one per subset).
Genomes with identical core SNPs (distance 0) are collapsed into one node whose area scales with the genome count.
Edge labels are SNP differences."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree

from common import OUT, SPECIES_COLORS, SUBSETS

TITLES = {"complex": "M. kansasii complex", "kansasii": "M. kansasii (sensu stricto)",
          "pseudokansasii": "M. pseudokansasii", "persicum": "M. persicum"}
samples = pd.read_csv(OUT / "samples_used.tsv", sep="\t").set_index("sample")

for name in SUBSETS:
    D = pd.read_csv(OUT / f"core_snp_distance_{name}.csv", index_col=0)
    # collapse identical haplotypes (rows of the distance matrix equal <=> distance 0)
    groups: list[list[str]] = []
    for n in D.index:
        for g in groups:
            if D.loc[n, g[0]] == 0:
                g.append(n)
                break
        else:
            groups.append([n])
    reps = [g[0] for g in groups]
    H = D.loc[reps, reps].to_numpy(dtype=float)
    k: int = len(reps)
    if k == 1:
        print(f"{name}: single haplotype, skipped")
        continue
    # zero entries would be dropped by the sparse MST; distinct haplotypes differ by >= 1 so this is safe
    T = minimum_spanning_tree(csr_matrix(H)).tocoo()
    G = nx.Graph()
    for i, g in enumerate(groups):
        sp = samples.loc[g, "species"].value_counts()
        G.add_node(i, n=len(g), species=sp.index[0])
    for i, j, w in zip(T.row.tolist(), T.col.tolist(), T.data.tolist()):
        G.add_edge(i, j, weight=w)
    assert G.number_of_edges() == k - 1 and nx.is_connected(G)
    # layout: tree path lengths, compressed (sqrt) so long and short edges are both visible
    Gs = nx.Graph((u, v, {"weight": np.sqrt(d["weight"])}) for u, v, d in G.edges(data=True))
    dist = dict(nx.all_pairs_dijkstra_path_length(Gs))
    pos = nx.kamada_kawai_layout(Gs, dist=dist)

    fig, ax = plt.subplots(figsize=(11, 9))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#888888", width=1.0)
    species_nodes = sorted({d["species"] for _, d in G.nodes(data=True)})
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=[SPECIES_COLORS.get(G.nodes[v]["species"], "#999999") for v in G],
                           node_size=[60 + 70 * G.nodes[v]["n"] for v in G], edgecolors="black", linewidths=0.6)
    thr = 0 if k <= 60 else float(np.percentile([d["weight"] for *_, d in G.edges(data=True)], 75))
    nx.draw_networkx_edge_labels(G, pos, ax=ax, font_size=6, rotate=False, bbox={"alpha": 0, "pad": 0},
                                 edge_labels={(u, v): f"{int(d['weight'])}" for u, v, d in G.edges(data=True)
                                              if d["weight"] > thr})
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc=SPECIES_COLORS.get(sp, "#999999"), mec="black", label=sp)
                       for sp in species_nodes], title="species", loc="best", frameon=False)
    nvar = len(pd.read_csv(OUT / f"core_snp_sites_{name}.csv"))
    ax.set_title(f"{TITLES[name]}: MST of core-genome SNPs\n{len(D)} genomes, {k} haplotypes, {nvar} variable core sites "
                 f"(edge labels = SNP differences)", fontsize=11)
    ax.axis("off")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"mst_{name}.{ext}", dpi=200)
    plt.close(fig)
    print(f"{name}: wrote mst_{name}.png/pdf ({len(D)} genomes, {k} haplotypes)")
