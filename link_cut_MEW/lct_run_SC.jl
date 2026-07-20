# lct_run_NH.jl — LCT-based NH runner
# include("../choniverse.jl")
include("../energy_functions.jl")


#
# Drop-in replacement for Marked_edges/bludef_uniform_to_doob.jl.
# Uses LCTState internally for fast MCMC steps but returns the SAME
# serialization format: (initial_partition, flips, tree, marked_edges)
# so existing .jls files and deserialization code are fully compatible.
#

include("lct_mew.jl")

# Wrap the original init utilities to avoid name collisions with lct_mew.jl
module BeanoInit
    using ProgressBars
    include(joinpath(@__DIR__, "..", "Marked_edges", "beano2.2_WI.jl"))
end

using JSON, DataFrames, Serialization, Graphs, ProgressBars

# ── helpers ───────────────────────────────────────────────────────────────────

function lct_state_to_simplegraph(s::LCTState)
    t = SimpleGraph(s.n)
    for v in 1:s.n
        for w in s.tree_adj[v]
            w > v && add_edge!(t, v, w)
        end
    end
    return t
end

function ntd_to_partition_dict(ntd::Vector{Int64})
    Dict{Int, Int}(i => ntd[i] for i in eachindex(ntd))
end

function compute_flip_ntd(old_ntd::Vector{Int64}, new_ntd::Vector{Int64})
    flips = Dict{Int, Int}()
    for i in eachindex(old_ntd)
        old_ntd[i] != new_ntd[i] && (flips[i] = new_ntd[i])
    end
    return flips
end

# ── main ──────────────────────────────────────────────────────────────────────

const K          = 7
const EPSILON    = 0.02
const N_ITERS    = 10_000

function main(; initialization=nothing)

    data  = JSON.parsefile("SC/SC_dual_graph_stripped.json")
    nodes = data["nodes"]
    links = data["links"]

    is_boundary = [node["boundary_node"] for node in nodes]
    boundary_lengths = zeros(length(nodes))
    for i in 1:length(nodes)
        if is_boundary[i]
            boundary_lengths[i] = nodes[i]["boundary_perim"]
        end
    end
    areas = [node["area"] for node in nodes]

    perim_dict = Dict()
    g = Graphs.SimpleGraph(length(nodes))
    for link in links
        u, v = link["source"], link["target"]
        add_edge!(g, simple_edge(u + 1, v + 1))
        perim_dict[simple_edge(u + 1, v + 1)] = link["shared_perim"]
    end

    df_dict = Dict()
    for node in nodes
        for (key, value) in node
            key == "id" && continue
            if !haskey(df_dict, key)
                df_dict[key] = []
            end
            push!(df_dict[key], value)
        end
    end
    max_len = maximum(length(v) for v in values(df_dict))
    for (key, vals) in df_dict
        while length(vals) < max_len
            push!(vals, missing)
        end
    end
    df = DataFrame(df_dict)
    insertcols!(df, 1, :id => 1:nrow(df))

    if isnothing(initialization)
        districts = nothing
        while isnothing(districts)
            districts = BeanoInit.find_k_partition(g, df, K, "population", 0.01)
        end
        tree_sg, marked_edges = BeanoInit.partition_to_tree_marked_edges(g, districts)
    else
        tree_sg      = initialization[1]
        marked_edges = initialization[2]
    end

    state     = from_simplegraph(g, tree_sg, marked_edges, K)
    pop_ideal = sum(df.population) / K
    score_cur = calculate_score(g, state.node_to_dist, K)

    # swap which energy function drives the chain by changing this line only
    # energy_fn = make_polsby_popper_energy(1000.0, areas, boundary_lengths, perim_dict, K)
    # energy_fn = make_cuts_energy(0.1, 600) 
    # energy_fn = make_county_splits_energy(0.4, df[!, "COUNTYFP"])
    energy_fn = make_combined_energy(0.1, df[!, "COUNTYFP"], 0.1, 10, 600)

    initial_partition = ntd_to_partition_dict(state.node_to_dist)
    current_ntd       = copy(state.node_to_dist)
    flips             = Vector{Dict{Int, Int}}()

    for _ in ProgressBar(1:N_ITERS)
        score_cur, accepted = mcmc_step!(
            state, df, EPSILON, energy_fn, score_cur, pop_ideal
        )
        flip = compute_flip_ntd(current_ntd, state.node_to_dist)
        push!(flips, flip)
        if accepted
            current_ntd = copy(state.node_to_dist)
        end
    end

    final_tree = lct_state_to_simplegraph(state)
    return initial_partition, flips, final_tree, copy(state.marked_edges)
end