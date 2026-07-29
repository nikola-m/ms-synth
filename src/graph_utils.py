"""
graph_utils.py
==============
Build a directed weighted graph from the CT-MSM generator matrix Q and
provide edge-weight thresholding routines.

Scientific motivation
---------------------
Following Mirkov et al. (2026), MS progression can be abstracted as
progressive edge removal in a directed connectivity graph.  Here the
graph is not a structural connectome but the *state-transition graph*
of the CT-MSM: nodes are clinical macrostates and edge weight λ_{ij}
quantifies the hazard of transitioning from state i to state j.

Progressive thresholding of this graph (removing weak edges as θ
increases) abstracts the loss of dynamical reversibility observed as
MS transitions from RRMS to SPMS (Kannan et al. 2017).

Key functions
-------------
build_graph        : Q  →  nx.DiGraph with 'weight' edge attribute
threshold_graph    : (G, θ) → thresholded DiGraph
normalise_weights  : map weights to [0, 1] for 'normalized_prob' mode
get_weight_range   : return (w_min, w_max) for the off-diagonal entries
"""

from __future__ import annotations

from typing import Any

import networkx as nx
import numpy as np

from .utils import get_logger

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------


def build_graph(
    Q: np.ndarray,
    state_labels: dict[int, str] | list[str] | None = None,
    include_self_loops: bool = False,
) -> nx.DiGraph:
    """Construct a directed weighted graph from the generator matrix Q.

    Nodes correspond to disease macrostates; directed edge i → j carries
    the transition intensity λ_{ij} = Q[i, j] (i ≠ j) as its weight.

    Parameters
    ----------
    Q : np.ndarray, shape (n, n)
        Infinitesimal generator matrix.  Off-diagonal entries must be ≥ 0.
    state_labels : dict[int, str] or list[str], optional
        Human-readable labels for each node (used for visualisation).
    include_self_loops : bool
        Whether to add self-loops with weight |Q[i, i]|.  Default False.

    Returns
    -------
    G : nx.DiGraph
        Nodes labelled 0 … n-1 with optional ``label`` attribute.
        Edges carry ``weight`` = λ_{ij}.
    """
    n = Q.shape[0]
    G = nx.DiGraph()

    # Add nodes
    for i in range(n):
        attrs: dict[str, Any] = {"state_idx": i}
        if state_labels is not None:
            if isinstance(state_labels, dict):
                attrs["label"] = state_labels.get(i, str(i))
            else:
                attrs["label"] = state_labels[i] if i < len(state_labels) else str(i)
        G.add_node(i, **attrs)

    # Add directed edges
    for i in range(n):
        for j in range(n):
            w = float(Q[i, j])
            if i == j:
                if include_self_loops:
                    G.add_edge(i, j, weight=abs(w))
                continue
            if w > 0:
                G.add_edge(i, j, weight=w)

    logger.debug(
        "Built DiGraph: %d nodes, %d edges  (max_weight=%.4f)",
        G.number_of_nodes(),
        G.number_of_edges(),
        max((d["weight"] for _, _, d in G.edges(data=True)), default=0.0),
    )
    return G


# ---------------------------------------------------------------------------
# Weight normalisation
# ---------------------------------------------------------------------------


def normalise_weights(G: nx.DiGraph, method: str = "max") -> nx.DiGraph:
    """Return a copy of G with weights normalised to [0, 1].

    Parameters
    ----------
    G : nx.DiGraph
    method : str
        "max"   – divide by max weight across all edges.
        "row"   – divide each edge by the sum of out-weights of its source
                  node (row-normalisation, converts to transition probs).

    Returns
    -------
    G_norm : nx.DiGraph (copy)
    """
    G_norm = G.copy()
    weights = [d["weight"] for _, _, d in G_norm.edges(data=True)]
    if not weights:
        return G_norm

    if method == "max":
        w_max = max(weights)
        if w_max > 0:
            for u, v in G_norm.edges():
                G_norm[u][v]["weight"] /= w_max

    elif method == "row":
        for node in G_norm.nodes():
            out_edges = list(G_norm.out_edges(node, data=True))
            total = sum(d["weight"] for _, _, d in out_edges)
            if total > 0:
                for u, v, d in out_edges:
                    G_norm[u][v]["weight"] = d["weight"] / total
    else:
        raise ValueError(f"Unknown normalisation method: '{method}'")

    return G_norm


# ---------------------------------------------------------------------------
# Thresholding
# ---------------------------------------------------------------------------


def threshold_graph(
    G: nx.DiGraph,
    theta: float,
    method: str = "absolute",
) -> nx.DiGraph:
    """Remove edges whose weight falls below threshold θ.

    Parameters
    ----------
    G : nx.DiGraph
        Original directed graph (weights are NOT modified in place).
    theta : float
        Threshold value.
    method : str
        "absolute"           – keep edges with weight ≥ θ.
        "fraction_strongest" – keep the top fraction θ of edges by weight
                               (θ ∈ (0, 1]; θ=1 → keep all).
        "normalized_prob"    – G must already have row-normalised weights;
                               keep edges with weight ≥ θ.

    Returns
    -------
    G_thresh : nx.DiGraph
        Subgraph containing only nodes and surviving edges.  All original
        nodes are preserved (isolated nodes may result).
    """
    if method == "absolute":
        edges_to_keep = [
            (u, v) for u, v, d in G.edges(data=True) if d["weight"] >= theta
        ]
    elif method == "fraction_strongest":
        all_edges = sorted(
            G.edges(data=True), key=lambda e: e[2]["weight"], reverse=True
        )
        n_keep = max(1, int(np.ceil(len(all_edges) * theta)))
        edges_to_keep = [(u, v) for u, v, _ in all_edges[:n_keep]]
    elif method == "normalized_prob":
        edges_to_keep = [
            (u, v) for u, v, d in G.edges(data=True) if d["weight"] >= theta
        ]
    else:
        raise ValueError(f"Unknown threshold method: '{method}'")

    # Build subgraph preserving all nodes
    G_thresh = nx.DiGraph()
    G_thresh.add_nodes_from(G.nodes(data=True))
    for u, v in edges_to_keep:
        G_thresh.add_edge(u, v, **G[u][v])

    logger.debug(
        "Threshold θ=%.4f (%s): %d/%d edges kept",
        theta,
        method,
        G_thresh.number_of_edges(),
        G.number_of_edges(),
    )
    return G_thresh


def get_weight_range(G: nx.DiGraph) -> tuple[float, float]:
    """Return (min, max) of edge weights in G.

    Parameters
    ----------
    G : nx.DiGraph

    Returns
    -------
    (w_min, w_max) : tuple[float, float]
        Returns (0.0, 0.0) if G has no edges.
    """
    weights = [d["weight"] for _, _, d in G.edges(data=True)]
    if not weights:
        return 0.0, 0.0
    return float(min(weights)), float(max(weights))


def build_theta_grid(
    G: nx.DiGraph,
    n_theta: int = 100,
    theta_min: float = 0.0,
    theta_max: float | None = None,
) -> np.ndarray:
    """Build a linearly-spaced threshold grid spanning the weight range of G.

    Parameters
    ----------
    G : nx.DiGraph
    n_theta : int
    theta_min : float
    theta_max : float, optional
        If None, uses max edge weight in G.

    Returns
    -------
    theta_grid : np.ndarray, shape (n_theta,)
    """
    w_min, w_max = get_weight_range(G)
    if theta_max is None:
        theta_max = w_max
    theta_max = max(theta_max, theta_min + 1e-9)
    grid = np.linspace(theta_min, theta_max, n_theta)
    return grid


# ---------------------------------------------------------------------------
# Directed component helpers
# ---------------------------------------------------------------------------


def get_gscc(G: nx.DiGraph) -> frozenset[int]:
    """Return the node set of the Giant Strongly Connected Component (GSCC).

    Parameters
    ----------
    G : nx.DiGraph

    Returns
    -------
    frozenset of node indices belonging to the largest SCC.
    """
    sccs = list(nx.strongly_connected_components(G))
    if not sccs:
        return frozenset()
    return frozenset(max(sccs, key=len))


def get_in_component(G: nx.DiGraph, source_nodes: list[int]) -> set[int]:
    """Return all nodes that can reach at least one node in ``source_nodes``.

    This implements the "giant in-component" relative to RRMS-like starting
    states: nodes from which the progressive sink states are reachable.

    Parameters
    ----------
    G : nx.DiGraph
    source_nodes : list[int]
        Representative starting states (e.g., mild RRMS macrostates).

    Returns
    -------
    set of int
    """
    in_comp: set[int] = set()
    for s in source_nodes:
        if s in G:
            in_comp.update(nx.ancestors(G, s))
            in_comp.add(s)
    return in_comp


def get_out_component(G: nx.DiGraph, source_nodes: list[int]) -> set[int]:
    """Return all nodes reachable from at least one node in ``source_nodes``.

    Parameters
    ----------
    G : nx.DiGraph
    source_nodes : list[int]

    Returns
    -------
    set of int
    """
    out_comp: set[int] = set()
    for s in source_nodes:
        if s in G:
            out_comp.update(nx.descendants(G, s))
            out_comp.add(s)
    return out_comp


def condensation_dag(G: nx.DiGraph) -> nx.DiGraph:
    """Return the condensation DAG of G (one super-node per SCC).

    Parameters
    ----------
    G : nx.DiGraph

    Returns
    -------
    nx.DiGraph
        Each node carries a ``members`` attribute listing the original
        node indices in that SCC.
    """
    cond = nx.condensation(G)
    # Add member list attribute
    scc_list = list(nx.strongly_connected_components(G))
    mapping = {}
    for scc in scc_list:
        for node in scc:
            mapping[node] = frozenset(scc)
    return cond
