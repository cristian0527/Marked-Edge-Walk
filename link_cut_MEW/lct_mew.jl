# lct_mew.jl — Marked Edge Walk using Link-Cut Trees
#
# Drop-in algorithmic replacement for beano2.3.jl.
# Key changes vs. the SimpleGraph version (see profiling notes):
#
#   • Eliminates deepcopy(tree) each iteration (~45 % of runtime in profile2.jlprof).
#     Proposals are applied in-place to the LCT; undone with inverse link/cut on reject.
#   • find_cycle_edges: O(log n) LCT path query instead of O(n) cycle_basis().
#   • Node-to-district lookup: O(1) array instead of O(n) connected_components().
#   • calculate_score caches the MH weight across accepted steps, so only k Cholesky
#     calls are needed per iteration (vs. 2k in the original).
#   • all_parts_edgelists: O(|E(g)|) single pass instead of O(n_district^2) combinations.

include("splaytrees.jl")
include("linkcuttrees.jl")

using Graphs, StatsBase, LinearAlgebra, SparseArrays, DataFrames

using ProgressBars

const Edge64 = Graphs.SimpleGraphs.SimpleEdge{Int64}

# ── canonical edge ─────────────────────────────────────────────────────────

simple_edge(u::Int, v::Int) = u < v ? Edge64(u, v) : Edge64(v, u)
simple_edge(e::Edge64)      = src(e) < dst(e) ? e : Edge64(dst(e), src(e))

# ── state ──────────────────────────────────────────────────────────────────

"""
Full MCMC state for the LCT-based Marked Edge Walk.

  g            — dual graph (immutable throughout the chain)
  lct          — link-cut tree representing the current spanning tree T
  tree_adj     — adjacency-list copy of T; kept in sync with lct for O(1)
                 degree queries and BFS without touching LCT internals
  nontree_edges — E(g) \\ E(T); sampled to pick edge_plus each step
  marked_edges — the k-1 tree edges whose removal partitions T into k districts
  node_to_dist — node index → district label (1..k); updated on accept
  k, n         — number of districts / nodes
"""
mutable struct LCTState
    g             :: SimpleGraph{Int64}
    lct           :: LinkCutTree{Int64}
    tree_adj      :: Vector{Set{Int64}}
    nontree_edges :: Set{Edge64}
    marked_edges  :: Vector{Edge64}
    node_to_dist  :: Vector{Int64}
    k             :: Int64
    n             :: Int64
end

# ── LCT link / cut wrappers ───────────────────────────────────────────────

function link_edge!(s::LCTState, e::Edge64)
    u, v = src(e), dst(e)
    evert!(s.lct.nodes[u])
    link!(s.lct.nodes[u], s.lct.nodes[v])
    push!(s.tree_adj[u], v)
    push!(s.tree_adj[v], u)
    delete!(s.nontree_edges, e)
end

function cut_edge!(s::LCTState, e::Edge64)
    u, v = src(e), dst(e)
    evert!(s.lct.nodes[u])
    cut!(s.lct.nodes[v])
    delete!(s.tree_adj[u], v)
    delete!(s.tree_adj[v], u)
    push!(s.nontree_edges, e)
end

# ── cycle finding via LCT path query ─────────────────────────────────────

"""
Return the edges on the unique tree path from u to v.
After conceptually adding edge_plus=(u,v), these edges plus edge_plus form
the cycle.  Uses O(log n) LCT expose operations instead of O(n) cycle_basis.
"""
function find_cycle_edges(lct::LinkCutTree, u::Int, v::Int)
    evert!(lct.nodes[u])          # re-root at u
    expose!(lct.nodes[v])         # preferred path = u → v in one splay tree
    path  = traverseSubtree(lct.nodes[v], "in-order")   # shallowest → deepest
    verts = [nd.vertex for nd in path]
    return [simple_edge(verts[i], verts[i+1]) for i in 1:length(verts)-1]
end

# ── partition from tree adjacency ────────────────────────────────────────

"""
BFS over tree_adj, treating marked edges as district boundaries.
Returns a node_to_dist vector (labels 1..k).
Replaces connected_components(tree) + rem_edge!/add_edge! dance.
"""
function compute_node_to_dist(
    tree_adj   :: Vector{Set{Int64}},
    marked_set :: Set{Edge64},
    n          :: Int64
)
    node_to_dist = zeros(Int64, n)
    dist = 0
    for start in 1:n
        node_to_dist[start] != 0 && continue
        dist += 1
        queue = [start]
        while !isempty(queue)
            v = pop!(queue)
            node_to_dist[v] != 0 && continue
            node_to_dist[v] = dist
            for w in tree_adj[v]
                node_to_dist[w] == 0 &&
                    simple_edge(v, w) ∉ marked_set &&
                    push!(queue, w)
            end
        end
    end
    return node_to_dist
end

# ── spanning-tree count (Cholesky) ───────────────────────────────────────

"""
Build per-district edge lists in a single O(|E(g)|) pass.
Much faster than the old combinations(part, 2) approach for large districts.
"""
function all_parts_edgelists(
    g            :: SimpleGraph{Int64},
    node_to_dist :: Vector{Int64},
    k            :: Int64
)
    lists = [Edge64[] for _ in 1:k]
    for e in edges(g)
        u, v = src(e), dst(e)
        d = node_to_dist[u]
        d == node_to_dist[v] && push!(lists[d], simple_edge(u, v))
    end
    return lists
end

"""
Log number of spanning trees of the subgraph defined by edgelist.
Uses your log_τ_sparse approach: stays fully sparse (no Matrix() conversion),
builds COO → SparseMatrixCSC directly, then sparse Cholesky + logdet.
Vertices are remapped to 1..m so the matrix is always district-sized.
"""
function log_τ_sparse(edgelist::Vector{Edge64})
    isempty(edgelist) && return 0.0

    unique_nodes = Set{Int64}()
    for e in edgelist
        push!(unique_nodes, src(e), dst(e))
    end
    m = length(unique_nodes)
    m <= 1 && return 0.0

    node_map = Dict(node => i for (i, node) in enumerate(unique_nodes))

    rows = Int64[]; cols = Int64[]; vals = Int64[]
    for e in edgelist
        i, j = node_map[src(e)], node_map[dst(e)]
        push!(rows, i, j); push!(cols, j, i); push!(vals, 1, 1)
    end

    A  = sparse(rows, cols, vals, m, m)
    L  = spdiagm(0 => vec(sum(A, dims=2))) - A
    Λ  = cholesky(L[2:end, 2:end])
    return logdet(Λ)
end

"""
Log spanning tree count of the quotient graph.
The quotient graph has k nodes; districts i and j are connected by
border_count_{ij} parallel edges (equivalently, a weighted edge of that weight).
For k=2 this reduces to log(border_count), matching the original formula.
For k>2 this is the correct Matrix-Tree computation — sum(log(border_counts))
would be wrong there.
"""
function log_τ_quotient(borders::Dict{Tuple{Int64,Int64},Int64}, k::Int64)
    k <= 1 && return 0.0
    isempty(borders) && return 0.0
    L = zeros(Float64, k, k)
    for ((i, j), w) in borders
        L[i, i] += w;  L[j, j] += w
        L[i, j] -= w;  L[j, i] -= w
    end
    k == 2 && return log(L[1, 1])   # fast path: 1×1 reduced Laplacian
    Λ = cholesky(Symmetric(L[2:end, 2:end]))
    return logdet(Λ)
end

"""
Count border edges between every adjacent district pair.
Returns Dict{(di, dj), count} with di < dj.
"""
function compute_borders(g::SimpleGraph{Int64}, node_to_dist::Vector{Int64})
    b = Dict{Tuple{Int64, Int64}, Int64}()
    for e in edges(g)
        di, dj = node_to_dist[src(e)], node_to_dist[dst(e)]
        di == dj && continue
        key = di < dj ? (di, dj) : (dj, di)
        b[key] = get(b, key, 0) + 1
    end
    return b
end

"""
Combined MH score: Σ log τ(district_i) + log τ(quotient graph).
The quotient graph term is the correct k-district generalisation of the border
length correction — for k=2 it equals log(border_count), for k>2 it is the
full Matrix-Tree determinant on the k-node contracted graph.
Cached across iterations; only recomputed for the proposed state each step.
"""
function calculate_score(
    g            :: SimpleGraph{Int64},
    node_to_dist :: Vector{Int64},
    k            :: Int64
)
    lists    = all_parts_edgelists(g, node_to_dist, k)
    τ_trees  = sum(log_τ_sparse(el) for el in lists)
    borders  = compute_borders(g, node_to_dist)
    τ_quotient = log_τ_quotient(borders, k)
    return τ_trees + τ_quotient
end

# ── transition probability ────────────────────────────────────────────────

"""
Metropolis-Hastings forward/reverse proposal ratio.
Uses degree arrays (O(1) lookup) instead of neighbors(tree, v) (O(degree)).
Signature matches the original transition_probability in beano2.3.jl.
"""
function transition_probability(
    cycle_edges      :: Vector{Edge64},
    edge_plus        :: Edge64,
    old_marked       :: Edge64,
    new_marked       :: Edge64,
    old_marked_edges :: Vector{Edge64},
    new_marked_edges :: Vector{Edge64},
    degree_old       :: Vector{Int64},
    degree_new       :: Vector{Int64}
)
    new_marked == edge_plus && return 0.0

    w, x = src(old_marked), dst(old_marked)
    u, v = src(new_marked), dst(new_marked)

    if Set([u, v]) == Set([w, x])
        d_u  = degree_old[u];  d_u_p = degree_new[u]
        d_v  = degree_old[v];  d_v_p = degree_new[v]
        pm   = (d_u + d_v) / (d_u_p + d_v_p) * d_u_p / d_u * d_v_p / d_v
    else
        shared = intersect([u, v], [w, x])[1]
        pm     = degree_new[shared] / degree_old[shared]
    end

    cycle = Set(cycle_edges)
    l     = length(setdiff(cycle, Set(old_marked_edges)))
    l_p   = length(setdiff(cycle, Set(new_marked_edges)))
    l_p == 0 && return 0.0

    return (l / l_p) * pm
end

# ── helpers ───────────────────────────────────────────────────────────────

function tally(df::DataFrame, col::String, node_to_dist::Vector{Int64}, k::Int64)
    vals   = df[:, col]
    totals = zeros(k)
    for i in eachindex(node_to_dist)
        totals[node_to_dist[i]] += vals[i]
    end
    return totals
end

function within_percent_of_ideal(vals::Vector{Float64}, ideal::Float64, epsilon::Float64)
    maximum(abs.((vals .- ideal) ./ ideal)) < epsilon
end





# ENERGY FUNCTION ───────────────────────────────────────────────────────────────


function county_splits(g, ntd, county_ids)
    splits = 0
    for cty in Set(county_ids)
        cty_nodes = findall(county_ids .== cty)
        cty_districts = Set(ntd[i] for i in cty_nodes)
        if length(cty_districts) > 1
            splits += 1
        end
    end
    return splits
end



# K is the number of districts
function make_county_splits_energy(
    beta :: Float64,
    county_ids :: Vector{Any}
)
    return function (g, new_ntd, old_ntd, county_ids, beta)
        cnty_splits_new = county_splits(g, new_ntd, county_ids)
        cnty_splits_old = county_splits(g, old_ntd, county_ids)
        return -beta * (cnty_splits_new - cnty_splits_old)
    end
end

function make_combined_energy(
    beta_county_splits :: Number,
    county_ids :: Vector{Any},
    beta_cuts :: Number,
    cty_target :: Number,
    target_cuts :: Number
    # perim_dit :: Dict{SimpleEdge, Float64}
)
    return function combined_energy(g, new_ntd, old_ntd)
        cnty_splits_new = county_splits(g, new_ntd, county_ids)
        cnty_splits_old = county_splits(g, old_ntd, county_ids)

        cuts_new = length(BeanoInit.cut_edges(new_ntd, g))
        cuts_old = length(BeanoInit.cut_edges(old_ntd, g))

        return -beta_county_splits * ((cnty_splits_new - cty_target)^2 - (cnty_splits_old - cty_target)^2) -
        beta_cuts * ((cuts_new - target_cuts)^2 - (cuts_old - target_cuts)^2)
    end
end




# ------------------ NEW SETUP --------------------
# | REPUBLICAN VOTESHARE | 
# | ATTICUS CODE | 

"""
Per-district score: continuous tent function peaked at threshold.
  share ≤ threshold :  y = share / threshold          (rises linearly to 1.0)
  share >  threshold :  y = 1 - slope_down*(share-threshold)  (falls slowly)
Every point of share below threshold is rewarded equally.
Every point above threshold is penalized at rate slope_down.
Keep slope_down small relative to 1/threshold so the packing penalty is
much weaker than the approaching-majority reward.
Default threshold=0.501, slope_down=0.5:
  gain 45% → 47%  = +0.040
  loss 70% → 75%  = -0.025  (~1.6× smaller)
"""
function _district_score(share::Float64;
                         threshold::Float64=0.501,
                         slope_down::Float64=0.5)
    share <= threshold ? share / threshold : 1.0 - slope_down * (share - threshold)
end
"""
Only scores the `n` least-Republican districts, each on a standard
_district_score tent peaking at `threshold` (default 55% R). The other
k-n districts are ignored entirely.
Useful when you only care about packing/protecting a fixed number of
safe-R seats and don't want to penalize or reward the rest.
"""
function rep_voteshare_score_topn(df::DataFrame, node_to_dist::Vector{Int64}, k::Int64, n::Int64;
                                  d_col::String="G24PREDHAR", r_col::String="G24PRERTRU",
                                  threshold::Float64=0.55, slope_down::Float64=0.5)
    d = tally(df, d_col, node_to_dist, k)
    r = tally(df, r_col, node_to_dist, k)
    shares = Float64[]
    for i in 1:k
        tot = d[i] + r[i]
        tot == 0 && continue
        push!(shares, r[i] / tot)
    end
    sort!(shares; rev=false)
    score = 0.0
    for i in 1:min(n, length(shares))
        score += _district_score(shares[i]; threshold, slope_down)
    end
    return score
end


function make_party_combined_energy(
    beta_county_splits :: Number,
    county_ids :: Vector{Any},
    beta_cuts :: Number,
    cty_target :: Number,
    target_cuts :: Number;
    df                   = nothing,
    k                    :: Int64   = 0,
    beta_voteshare :: Number = 0.0,
    n_top :: Number = 3,
    voteshare_threshold :: Float64 = 0.45,
    voteshare_slope_down :: Float64 = 0.5,
    d_col                :: String  = "G24PREDHAR",
    r_col                :: String  = "G24PRERTRU"

    # perim_dit :: Dict{SimpleEdge, Float64}
)
    if beta_voteshare != 0.0
        @assert !isnothing(df) && k > 0 "make_combined_energy: df and k are required when beta_voteshare != 0"
    end
    return function triple_combined_energy(g, new_ntd, old_ntd)
        cty_splits_new = county_splits(g, new_ntd, county_ids)
        cty_splits_old = county_splits(g, old_ntd, county_ids)

        cuts_new = length(BeanoInit.cut_edges(new_ntd, g))
        cuts_old = length(BeanoInit.cut_edges(old_ntd, g))

        energy = -beta_county_splits * ((cty_splits_new - cty_target)^2 - (cty_splits_old - cty_target)^2) -
                 beta_cuts          * ((cuts_new - target_cuts)^2 - (cuts_old - target_cuts)^2)
        if beta_voteshare != 0.0
            score_new = rep_voteshare_score_topn(df, new_ntd, k, n_top;
                                                  d_col=d_col, r_col=r_col,
                                                  threshold=voteshare_threshold, slope_down=voteshare_slope_down)
            score_old = rep_voteshare_score_topn(df, old_ntd, k, n_top;
                                                  d_col=d_col, r_col=r_col,
                                                  threshold=voteshare_threshold, slope_down=voteshare_slope_down)
            energy += beta_voteshare * (score_new - score_old)
        end
        return energy
    end
end




### GAUSSIAN REP VOTE SHARES

function rep_voteshare_score_vector_gaussian(df::DataFrame, node_to_dist::Vector{Int64}, k::Int64, targets::Vector{Float64}, uses::Vector{Symbol}, slope_down::Number=0.5;
                                  d_col::String="G24PREDHAR", r_col::String="G24PRERTRU")
    d = tally(df, d_col, node_to_dist, k)
    r = tally(df, r_col, node_to_dist, k)
    shares = Float64[]

    for i in 1:k
        tot = d[i] + r[i]
        tot == 0 && continue
        push!(shares, r[i] / tot)
    end

    sort!(shares; rev=false)

    score = 0.0

    for i in 1:length(shares)
        if uses[i] == :do_nothing
            continue
        elseif uses[i] == :equal
            score += shares[i] <= targets[i] ? shares[i] / targets[i] : 1.0 - slope_down * (shares[i] - targets[i])^2
        elseif uses[i] == :less
            continue
            # score += (1-shares[i]) <= (1-targets[i]) ? (1-shares[i]) / (1-targets[i]) : 1.0
        elseif uses[i] ==:greater
            continue
            # score += shares[i] <= targets[i] ? shares[i] / targets[i] : 1.0
        end
    end
    return score
end

function make_combined_gaussian_party_energy(
    county_ids          :: Vector{Any},
    beta_county_splits :: Number,
    beta_cuts           :: Number,
    cty_target          :: Number,
    target_cuts         :: Number;
    df                   = nothing,
    k                    :: Int64   = 0,
    targets             :: Vector{Float64},
    uses                :: Vector{Symbol},
    beta_voteshare      :: Number  = 0.0,
    voteshare_slope_down :: Number = 0.5,
    d_col                :: String  = "G24PREDHAR",
    r_col                :: String  = "G24PRERTRU"
)
    if beta_voteshare != 0.0
        @assert !isnothing(df) && k > 0 "make_combined_energy: df and k are required when beta_voteshare != 0"
    end
    return function combined_energy(g, new_ntd, old_ntd)
        cty_splits_new = county_splits(g, new_ntd, county_ids)
        cty_splits_old = county_splits(g, old_ntd, county_ids)
        cuts_new = count(e -> new_ntd[src(e)] != new_ntd[dst(e)], edges(g))
        cuts_old = count(e -> old_ntd[src(e)] != old_ntd[dst(e)], edges(g))
        energy = -beta_county_splits * ((cty_splits_new - cty_target)^2 - (cty_splits_old - cty_target)^2) -
                 beta_cuts          * ((cuts_new - target_cuts)^2 - (cuts_old - target_cuts)^2)
        if beta_voteshare != 0.0
            score_new = rep_voteshare_score_vector_gaussian(df, new_ntd, k, targets, uses, voteshare_slope_down;
                                                  d_col=d_col, r_col=r_col)
            score_old = rep_voteshare_score_vector_gaussian(df, old_ntd, k, targets, uses, voteshare_slope_down;
                                                  d_col=d_col, r_col=r_col)
            energy -= beta_voteshare * (score_new - score_old)
        end
        return energy
    end
end

# example of atticus'

"""
county_splits / cuts terms are Gaussian-style: log-density is
-beta*(deviation from a fixed target)^2, so they penalize moving away from
cty_target/target_cuts in either direction.
The optional voteshare term (beta_voteshare != 0) is different in kind: it's
exponential-tilted rather than Gaussian. rep_voteshare_score_topn already
peaks at voteshare_threshold via its own tent shape (see _district_score), so
the stationary density is just reweighted by exp(beta_voteshare * score) —
log-density contribution beta_voteshare*(score_new - score_old), monotonic in
score rather than penalizing distance from a target. (A Gaussian penalty
directly on each top-n district's raw share against its own target voteshare
is the natural next step if this doesn't converge well enough on its own.)
"""
# function make_combined_energy(
#    county_ids          :: Vector{Any},
#    beta_county_splits :: Number,
#    beta_cuts           :: Number,
#    cty_target          :: Number,
#    target_cuts         :: Number;
#    df                   = nothing,
#    k                    :: Int64   = 0,
#    beta_voteshare      :: Number  = 0.0,
#    n_top                :: Int64   = 3,
#    voteshare_threshold :: Float64 = 0.45,
#    voteshare_slope_down :: Float64 = 0.5,
#    d_col                :: String  = "G24PREDHAR",
#    r_col                :: String  = "G24PRERTRU",
#)
#    if beta_voteshare != 0.0
#        @assert !isnothing(df) && k > 0 "make_combined_energy: df and k are required when beta_voteshare != 0"
#    end
#    return function combined_energy(g, new_ntd, old_ntd)
#        cty_splits_new = county_splits(g, new_ntd, county_ids)
#        cty_splits_old = county_splits(g, old_ntd, county_ids)
#        cuts_new = count(e -> new_ntd[src(e)] != new_ntd[dst(e)], edges(g))
#        cuts_old = count(e -> old_ntd[src(e)] != old_ntd[dst(e)], edges(g))
#        energy = -beta_county_splits * ((cty_splits_new - cty_target)^2 - (cty_splits_old - cty_target)^2) -
#                 beta_cuts          * ((cuts_new - target_cuts)^2 - (cuts_old - target_cuts)^2)
#        if beta_voteshare != 0.0
#            score_new = rep_voteshare_score_topn(df, new_ntd, k, n_top;
#                                                  d_col=d_col, r_col=r_col,
#                                                  threshold=voteshare_threshold, slope_down=voteshare_slope_down)
#            score_old = rep_voteshare_score_topn(df, old_ntd, k, n_top;
#                                                  d_col=d_col, r_col=r_col,
#                                                  threshold=voteshare_threshold, slope_down=voteshare_slope_down)
#            energy += beta_voteshare * (score_new - score_old)
#        end
#        return energy
#    end
#end


# ── energy function ───────────────────────────────────────────────────────

"""
mcmc_step! / run_chain! take a single `energy_fn` argument with signature

    energy_fn(g, new_ntd, old_ntd) -> Float64   # log(E_new / E_old)

so switching which energy function drives the chain means changing one
line at the call site (which `make_*_energy` you call), not the signatures
of mcmc_step!/run_chain! or every place that calls them. Build one with
make_cuts_energy(...) or make_polsby_popper_energy(...) below, or write
your own closure with the same (g, new_ntd, old_ntd) -> Float64 shape.
"""

"""
Cut-count energy  E ∝ exp(-beta * (n_cuts - target_cuts)^2).
Returns a closure usable as `energy_fn`.
"""
function make_cuts_energy(beta::Float64, target_cuts::Int64)
    return function (g::SimpleGraph{Int64}, new_ntd::Vector{Int64}, old_ntd::Vector{Int64})
        cuts_new = count(e -> new_ntd[src(e)] != new_ntd[dst(e)], edges(g))
        cuts_old = count(e -> old_ntd[src(e)] != old_ntd[dst(e)], edges(g))
        -beta * ((cuts_new - target_cuts)^2 - (cuts_old - target_cuts)^2)
    end
end

"""
Mean Polsby-Popper compactness (4πA/P², averaged over districts) for a
partition given as node_to_dist. Mirrors polsby_popper(areas, boundary_lengths,
perim_dict, ptition, g) in beano2.2_WI.jl, but works directly off node_to_dist
(no connected_components / cut_edges dance) so it can be called from
mcmc_step! without rebuilding a partition Dict.

  areas             — per-node area (n-vector), e.g. node["area"] from the
                       dual-graph JSON
  boundary_lengths  — per-node length of *state*-boundary perimeter, 0 for
                       interior nodes; see analysis.jl for how this is built
                       from node["boundary_node"] / node["boundary_perim"]
  perim_dict        — Dict keyed by simple_edge(u, v) -> shared perimeter
                       length between u and v, e.g. link["shared_perim"]

Each cut edge's shared perimeter borders exactly the two districts it
separates, so it's added to both district_boundaries entries, not split
across all k districts.
"""
function mean_polsby_popper(
    g                :: SimpleGraph{Int64},
    node_to_dist     :: Vector{Int64},
    k                :: Int64,
    areas            :: Vector{Float64},
    boundary_lengths :: Vector{Float64},
    perim_dict       :: Dict
)
    district_areas      = zeros(Float64, k)
    district_boundaries = zeros(Float64, k)

    for v in eachindex(node_to_dist)
        d = node_to_dist[v]
        district_areas[d]      += areas[v]
        district_boundaries[d] += boundary_lengths[v]
    end

    for e in edges(g)
        u, v   = src(e), dst(e)
        di, dj = node_to_dist[u], node_to_dist[v]
        di == dj && continue
        shared = perim_dict[simple_edge(u, v)]
        district_boundaries[di] += shared
        district_boundaries[dj] += shared
    end

    pps = (4 * pi * district_areas[d] / district_boundaries[d]^2 for d in 1:k)
    return sum(pps) / k
end

"""
Mean Polsby-Popper compactness energy, in the same target-seeking form as
make_cuts_energy: E ∝ exp(-beta * (mean_pp - target_pp)^2).

target_pp defaults to 1/3, matching polsby_popper_energy_function in
beano2.2_WI.jl. areas / boundary_lengths / perim_dict follow the same
convention as polsby_popper() there (see analysis.jl for how they're loaded
from the dual-graph JSON). Returns a closure usable as `energy_fn`.
"""
function make_polsby_popper_energy(
    beta             :: Float64,
    areas            :: Vector{Float64},
    boundary_lengths :: Vector{Float64},
    perim_dict       :: Dict,
    k                :: Int64;
    target_pp        :: Float64 = 1/3
)
    return function (g::SimpleGraph{Int64}, new_ntd::Vector{Int64}, old_ntd::Vector{Int64})
        pp_new = mean_polsby_popper(g, new_ntd, k, areas, boundary_lengths, perim_dict)
        pp_old = mean_polsby_popper(g, old_ntd, k, areas, boundary_lengths, perim_dict)
        -beta * ((pp_new - target_pp)^2 - (pp_old - target_pp)^2)
    end
end

# ── initialisation ────────────────────────────────────────────────────────

"""
Wilson's algorithm for a uniform random spanning tree.
Returns a Set{Edge64} of tree edges.
"""
function wilsons(g::SimpleGraph{Int64})
    n       = nv(g)
    visited = falses(n)
    spanning = Set{Edge64}()
    start    = rand(1:n)
    visited[start] = true

    while count(visited) < n
        current = rand(findall(.!visited))
        path    = [current]
        while !visited[current]
            next = rand(all_neighbors(g, current))
            idx  = findfirst(==(next), path)
            isnothing(idx) ? push!(path, next) : resize!(path, idx)
            current = next
        end
        for i in 1:length(path)-1
            visited[path[i]] = true
            push!(spanning, simple_edge(path[i], path[i+1]))
        end
    end
    return spanning
end

"""
Construct an LCTState from a dual graph, an initial spanning tree edge set,
and k-1 marked edges.  Links tree edges into the LCT in BFS order to avoid
same-component link errors.
"""
function build_lct_state(
    g             :: SimpleGraph{Int64},
    tree_edge_set :: Set{Edge64},
    marked_edges  :: Vector{Edge64},
    k             :: Int64
)
    n = nv(g)

    # adjacency list for the spanning tree
    tree_adj = [Set{Int64}() for _ in 1:n]
    for e in tree_edge_set
        push!(tree_adj[src(e)], dst(e))
        push!(tree_adj[dst(e)], src(e))
    end

    # build LCT in BFS order so link! never joins an already-connected pair
    lct     = LinkCutTree{Int64}(n)
    visited = falses(n)
    visited[1] = true
    queue   = [1]
    while !isempty(queue)
        v = popfirst!(queue)
        for w in tree_adj[v]
            visited[w] && continue
            visited[w] = true
            evert!(lct.nodes[w])
            link!(lct.nodes[w], lct.nodes[v])
            push!(queue, w)
        end
    end

    # non-tree edges for sampling edge_plus
    nontree = Set{Edge64}()
    for e in edges(g)
        se = simple_edge(src(e), dst(e))
        se ∉ tree_edge_set && push!(nontree, se)
    end

    node_to_dist = compute_node_to_dist(tree_adj, Set(marked_edges), n)

    return LCTState(g, lct, tree_adj, nontree, copy(marked_edges), node_to_dist, k, n)
end

# ── single MCMC step ──────────────────────────────────────────────────────

"""
One Metropolis-Hastings step.  Returns (updated_score, accepted::Bool).

No deepcopy of the tree: proposals are applied in-place to the LCT and
tree_adj, then reversed with inverse link/cut operations on rejection.
old_node_to_dist (a Vector{Int64} copy, ~n × 8 bytes) is the only
allocation that scales with graph size, replacing the O(n × degree) copy
of the full SimpleGraph adjacency list in the original code.

energy_fn is a (g, new_ntd, old_ntd) -> Float64 closure built by
make_cuts_energy / make_polsby_popper_energy / a custom closure — see the
"energy function" section above.
"""
function mcmc_step!(
    s                :: LCTState,
    df               :: DataFrame,
    epsilon          :: Float64,
    energy_fn,
    score_cur        :: Float64,
    pop_ideal        :: Float64;
    max_tries        :: Int     = 200
)
    n, k = s.n, s.k

    # snapshot degrees and partition BEFORE any proposal
    degree_old    = [length(s.tree_adj[v]) for v in 1:n]
    old_node_to_dist = copy(s.node_to_dist)   # ~n×8 bytes; much cheaper than copy(tree)

    local cycle_edges, edge_plus, edge_minus, old_marked, new_marked, new_ntd, m_idx

    found = false
    for tries in 1:max_tries

        # ── cycle basis step ──────────────────────────────────────────────
        edge_plus   = rand(s.nontree_edges)
        u, v        = src(edge_plus), dst(edge_plus)
        cycle_edges = find_cycle_edges(s.lct, u, v)
        marked_set  = Set(s.marked_edges)
        possible    = setdiff(Set(cycle_edges), marked_set)
        if isempty(possible)
            continue
        end
        edge_minus = rand(possible)
        cut_edge!(s, edge_minus)
        link_edge!(s, edge_plus)

        # ── marked edge step ──────────────────────────────────────────────
        old_marked = rand(s.marked_edges)
        v1, v2     = src(old_marked), dst(old_marked)
        chosen_v   = rand((v1, v2))
        new_marked = simple_edge(chosen_v, rand(s.tree_adj[chosen_v]))
        m_idx      = findfirst(==(old_marked), s.marked_edges)
        s.marked_edges[m_idx] = new_marked

        # ── population balance check ──────────────────────────────────────
        new_ntd = compute_node_to_dist(s.tree_adj, Set(s.marked_edges), n)
        pops    = tally(df, "population", new_ntd, k)
        if within_percent_of_ideal(pops, pop_ideal, epsilon)
            found = true
            break
        end

        # ── undo this attempt ─────────────────────────────────────────────
        s.marked_edges[m_idx] = old_marked
        cut_edge!(s, edge_plus)
        link_edge!(s, edge_minus)

        if tries == max_tries
            println("max tries reached")
        end
    end

    !found && return score_cur, false

    # ── compute acceptance probability ────────────────────────────────────
    degree_new = [length(s.tree_adj[v]) for v in 1:n]

    # reconstruct old_marked_edges vector (before the marked step)
    old_m_vec        = copy(s.marked_edges)
    old_m_vec[m_idx] = old_marked

    # original cycle_basis_step includes edge_plus in cycle_edges (it's a vertex
    # cycle, so edge_plus closes the loop); replicate that here so the l/l_p
    # ratio in transition_probability matches exactly
    a1 = transition_probability(
        [cycle_edges; edge_plus], edge_plus, old_marked, new_marked,
        old_m_vec, s.marked_edges, degree_old, degree_new
    )

    score_new = calculate_score(s.g, new_ntd, k)
    a2        = score_cur - score_new
    a3_log    = energy_fn(s.g, new_ntd, old_node_to_dist)

    # println(a3_log)

    # ── accept / reject ───────────────────────────────────────────────────
    if log(rand(Float64)) < log(a1) + a2 + a3_log
        s.node_to_dist = new_ntd
        return score_new, true
    else
        s.marked_edges[m_idx] = old_marked
        cut_edge!(s, edge_plus)
        link_edge!(s, edge_minus)
        return score_cur, false
    end
end

# ── main chain loop ───────────────────────────────────────────────────────

"""
Run the LCT-based Marked Edge Walk for num_iterations steps.
Returns a Vector of node_to_dist snapshots (one per step).

Arguments:
  s                — LCTState built by build_lct_state()
  df               — DataFrame with at least a "population" column
  num_iterations   — number of MCMC steps
  energy_fn        — (g, new_ntd, old_ntd) -> Float64 closure, e.g. from
                      make_cuts_energy(...) or make_polsby_popper_energy(...)
  epsilon          — population balance tolerance (fraction of ideal)
"""
function run_chain!(
    s                :: LCTState,
    df               :: DataFrame,
    num_iterations   :: Int,
    energy_fn;
    epsilon          :: Float64 = 0.05)
    pop_ideal  = sum(df[:, "population"]) / s.k
    score_cur  = calculate_score(s.g, s.node_to_dist, s.k)

    partitions = Vector{Vector{Int64}}(undef, num_iterations)
    accepts    = 0

    
    for i in ProgressBar(1:num_iterations, printing_delay=0.001) # progress bar displayed in terminal
        score_cur, accepted = mcmc_step!(
            s, df, epsilon, energy_fn, score_cur, pop_ideal
        )
        accepted && (accepts += 1)
        partitions[i] = copy(s.node_to_dist)
    end

    println("Acceptance rate: $(round(accepts / num_iterations, digits = 4))")
    return partitions
end

# ── convenience: build state from a SimpleGraph tree ─────────────────────

"""
Convert a Graphs.SimpleGraph spanning tree + marked edges (the format used
in beano2.3.jl) into an LCTState.  Useful for warm-starting from an existing
chain or testing against the original code.
"""
function from_simplegraph(
    g            :: SimpleGraph{Int64},
    tree         :: SimpleGraph{Int64},
    marked_edges :: Vector{Edge64},
    k            :: Int64
)
    tree_edge_set = Set(simple_edge(src(e), dst(e)) for e in edges(tree))
    return build_lct_state(g, tree_edge_set, marked_edges, k)
end
