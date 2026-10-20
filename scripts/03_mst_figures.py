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
from matplotlib.patches import Rectangle
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree

from common import MERGE_SNPS, OUT, SPECIES_COLORS, SUBSETS

TITLES = {"complex": "M. kansasii complex", "kansasii": "M. kansasii (sensu stricto)",
          "pseudokansasii": "M. pseudokansasii", "persicum": "M. persicum"}
ap = argparse.ArgumentParser()
ap.add_argument("--merge", type=int, default=None, help="merge nodes joined by MST edges of <= N SNPs (default: per subset)")
ap.add_argument("--subset", choices=list(SUBSETS), default=None, help="only this subset (default: all)")
args = ap.parse_args()


def node_area(n: int) -> float:
    """Marker area (pt^2) of a node: grows with its genome count only (not edges or degree), as sqrt so big nodes stay moderate."""
    return float(150.0 * np.sqrt(n))


def sample_number(sample: str) -> int:
    return int(sample.rsplit("-", 1)[1])  # Mkan329-001 -> 1


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

    big = len(D) > 150  # the complex figure only: bigger canvas for its crowded centre
    fig, ax = plt.subplots(figsize=(24, 20) if big else (11, 9) if k <= 30 else (16, 13))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#888888", width=1.0)
    species_nodes = sorted({d["species"] for _, d in G.nodes(data=True)})
    nodes = nx.draw_networkx_nodes(G, pos, ax=ax, node_color=[SPECIES_COLORS.get(G.nodes[v]["species"], "#999999") for v in G],
                           node_size=[node_area(G.nodes[v]["n"]) for v in G], edgecolors="black", linewidths=0.6)
    def node_label(g: list[str]) -> str:
        nums = [f"M{x:03d}" for x in sorted(sample_number(n) for n in g)]
        shown = nums[:8]
        lines = [", ".join(shown[i:i + 3]) for i in range(0, len(shown), 3)]
        if len(nums) > 8:
            lines.append(f"+{len(nums) - 8} more")
        return "\n".join(lines + ([f"(n={len(g)})"] if len(g) > 1 else []))
    pd.DataFrame({"node": range(k), "n": [len(g) for g in groups],
                  "species": [G.nodes[i]["species"] for i in range(k)],
                  "max_internal_snps": [int(D.loc[g, g].to_numpy().max()) for g in groups],
                  "numbers": [",".join(str(x) for x in sorted(sample_number(n) for n in g)) for g in groups],
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
        r_px = np.sqrt(node_area(G.nodes[i]["n"])) / 2 * fig.dpi / 72  # marker diameter = sqrt(area) pt
        cpx = ax.transData.transform(pos[i])
        # direction = middle of the widest angular gap between the node's tree neighbours (leaf: away from its parent)
        ang = np.sort([np.arctan2(*(ax.transData.transform(pos[j]) - cpx)[::-1]) for j in G[i]])
        gaps = np.diff(np.append(ang, ang[0] + 2 * np.pi))
        a = ang[int(np.argmax(gaps))] + gaps.max() / 2
        vec = np.array([np.cos(a), np.sin(a)])
        x, y = ax.transData.inverted().transform(cpx + vec * (r_px + 10 + 12 * max(G.degree(i) - 2, 0)))
        texts.append(ax.text(x, y, node_label(g), fontsize=6 if k > 30 else 7, ha="center", va="center",
                             # lists of big nodes get a framed box in the node's colour, joined to the node by a leader line
                             bbox={"facecolor": "white", "alpha": 0.85 if len(g) >= 4 else 0.7, "pad": 1.5 if len(g) >= 4 else 0.3,
                                   "edgecolor": SPECIES_COLORS.get(G.nodes[i]["species"], "#999999") if len(g) >= 4 else "none"}))
    # invisible squares the size of each node so labels keep off the discs (scatter extents are NaN in adjustText)
    node_boxes = []
    for v in range(k):
        half = np.sqrt(node_area(G.nodes[v]["n"])) / 2 * fig.dpi / 72
        cx, cy = ax.transData.transform(pos[v])
        (x0, y0), (x1, y1) = ax.transData.inverted().transform([(cx - half, cy - half), (cx + half, cy + half)])
        node_boxes.append(ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec="none")))
    adjust_text(texts, x=[pos[i][0] for i in range(k)], y=[pos[i][1] for i in range(k)], objects=[*node_boxes, *edge_txt.values()],
                ax=ax, expand=(1.1, 1.3), force_text=(0.3, 0.6), force_static=(0.5, 0.8), iter_lim=1000)
    # adjustText only softly avoids the node discs; push any label still touching a disc radially outward until it is clear
    fig.canvas.draw()
    centres = [np.asarray(ax.transData.transform(pos[v]), dtype=float) for v in range(k)]
    radii = [float(np.sqrt(node_area(G.nodes[v]["n"])) / 2 * fig.dpi / 72) for v in range(k)]
    inv = ax.transData.inverted()
    for _ in range(60):
        moved = False
        for t in texts:
            bb = t.get_window_extent()
            for cen, rad in zip(centres, radii):
                dx = cen[0] - min(max(cen[0], bb.x0), bb.x1)
                dy = cen[1] - min(max(cen[1], bb.y0), bb.y1)
                d = float(np.hypot(dx, dy))
                if d < rad + 3:
                    tc = np.array([(bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2])
                    vec = tc - cen
                    norm = float(np.hypot(*vec))
                    vec = vec / norm if norm > 1e-6 else np.array([0.0, -1.0])
                    t.set_position(tuple(inv.transform(tc + vec * max(rad + 3 - d, 4.0))))
                    moved = True
                    break
        if not moved:
            break
    # leader lines: every label that is not touching its node is joined to the node edge by a line to the nearest point of the label
    for t, cen, rad in zip(texts, centres, radii):
        bb = t.get_window_extent()
        nearest = np.array([min(max(cen[0], bb.x0), bb.x1), min(max(cen[1], bb.y0), bb.y1)])
        gap = float(np.hypot(*(nearest - cen)))
        if gap > rad + 18:  # only labels that are clearly detached from their node get a leader line
            start = cen + (nearest - cen) / gap * rad
            (x0, y0), (x1, y1) = inv.transform([tuple(start), tuple(nearest)])
            ax.add_line(Line2D([x0, x1], [y0, y1], color="#444444", lw=0.9, ls=(0, (1.5, 1.5)), zorder=1.5))  # dotted, so it is not mistaken for a tree edge
    size_leg = ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec="black", markersize=float(np.sqrt(node_area(m))),
                                         label=str(m)) for m in (1, 5, 10, 20) if m <= max(len(g) for g in groups) * 2],
                         title="genomes", loc="upper left", bbox_to_anchor=(1.01, 0.82), frameon=False, labelspacing=1.2)
    ax.add_artist(size_leg)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc=SPECIES_COLORS.get(sp, "#999999"), mec="black", label=sp)
                       for sp in species_nodes], title="species", loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    nvar = len(pd.read_csv(OUT / f"core_snp_sites_{name}.csv"))
    ax.set_title(f"{TITLES[name]}: MST of core-genome SNPs\n{len(D)} genomes, {k} nodes (merged at <= {merge} SNPs), {nvar} variable core sites\n"
                 f"PE/PPE, IS, phage and repeats masked; node labels = sample ID (MNNN); edge labels = SNP differences", fontsize=11)
    ax.axis("off")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"mst_{name}.{ext}", dpi=200)
    plt.close(fig)
    print(f"{name}: wrote mst_{name}.png/pdf ({len(D)} genomes, {len(hap)} haplotypes, {k} nodes at <= {merge} SNPs)")
