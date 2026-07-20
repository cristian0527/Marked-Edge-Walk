
include("lct_mew.jl")


# Wrap the original init utilities to avoid name collisions with lct_mew.jl
module BeanoInit
    using ProgressBars
    include(joinpath(@__DIR__, "..", "Marked_edges", "beano2.2_WI.jl"))
    # include("Marked_edges/beano2.2_WI.jl") # this if can't find Beano
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

# __ MAIN __ 
const K          = 6 # Louisiana
const EPSILON    = 0.05 # 0.01 before
# within 5% balance of the population

const N_ITERS    = 4_000_000


# LOAD GRAPH AND DF FROM 
#function load_graph_and_df()
#    data  = JSON.parsefile("Minnesota/mn.json")
#    nodes = data["nodes"]
#    adj   = data["adjacency"]
#    println("Loaded graph with $(length(nodes)) nodes.")
#    g = Graphs.SimpleGraph(length(nodes))
#    for (i, neighbors) in enumerate(adj)
#        for nb in neighbors
#            j = nb["id"] + 1   # 0-based → 1-based
#            i < j && add_edge!(g, simple_edge(i, j))
#        end
#    end
#    println("Graph has $(ne(g)) edges.")
#    df_dict = Dict()
#   for node in nodes
#       for (key, value) in node
#            key == "id" && continue
#            if !haskey(df_dict, key)
#                df_dict[key] = []
#            end
#            push!(df_dict[key], value)
#        end
#    end
#    max_len = maximum(length(v) for v in values(df_dict))
#    for (_, vals) in df_dict
#        while length(vals) < max_len
#            push!(vals, missing)
#        end
#    end
#    df = DataFrame(df_dict)
#    insertcols!(df, 1, :id => 1:nrow(df))
#    return g, df
#end

# ntd_to_tree
# markov edge



# main function
function main(; initialization=nothing)

    # data  = JSON.parsefile("LA/LA_dual_graph_stripped_links.json") # , allownan=true)
    data = JSON.parsefile("LA/LA_processed_precincts_Julia.json")
    nodes = data["nodes"]
    # links = data["links"]
    # println(data["links"])
    adjacency = data["adjacency"]

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
    for i in 1:length(nodes)
        u = i
        for nbr in adjacency[i]
            v = nbr["id"] + 1
            e = simple_edge(u, v)
            add_edge!(g, e)
            perim_dict[e] = nbr["shared_perim"]
        end
    end
    
    #for link in links
    #    u, v = link["source"], link["target"]
    #    add_edge!(g, simple_edge(u + 1, v + 1))
    #    perim_dict[simple_edge(u + 1, v + 1)] = link["shared_perim"]
    # end
    
    # what was output  
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
    # county_ids = df[!, "COUNTYFP"]
    # Set(county_ids)

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

  
    county_ids = df[!, "COUNTYFP"]
    # energy_fn = make_combined_energy(0.75, df[!, "COUNTYFP"], 0.015, 11, 700) 
    # LA WINNER RAHHH NOW FOR PARTY SHARE

    targetsLA = [0.3, 0.6, 0, 0, 0, 0]
    usesLA = [:equal, :equal, :do_nothing, :do_nothing, :do_nothing, :do_nothing]

    energy_fn = make_combined_gaussian_party_energy(county_ids, 0.75, 0.015, 11, 700; df=df, k=K, targets=targetsLA, uses=usesLA, beta_voteshare=550, voteshare_slope_down=0.5, d_col="G24PREDHAR", r_col="G24PRERTRU")
    # function make_combined_gaussian_party_energy(
    #county_ids          :: Vector{Any},
    #beta_county_splits :: Number,
    #beta_cuts           :: Number,
    #cty_target          :: Number,
    #target_cuts         :: Number;
    #df                   = nothing,
    #k                    :: Int64   = 0,
    #targets             :: Vector{Float64},
    #uses                :: Vector{Symbol},
    #beta_voteshare      :: Number  = 0.0,
    #voteshare_slope_down :: Number = 0.5,
    #d_col                :: String  = "G24PREDHAR",
    #r_col                :: String  = "G24PRERTRU",


    # cut edges looks good 
    # 0.5 -> 1
    # 1/2 and 13
    # decrease county number 
    # increase cutedges

    # make super combined 
    # 30 or less first one / over 70 for the one 
    # below 40 / above 60 for the rest
    # 


# function make_party_combined_energy(
    #beta_county_splits :: Number,
    #county_ids :: Vector{Any},
    #beta_cuts :: Number,
    #cty_target :: Number,
    #target_cuts :: Number;
    #df                   = nothing,
    #k                    :: Int64   = 0,
    #beta_voteshare :: Number = 0.0,
    #n_top :: Number = 3,
    #voteshare_threshold :: Float64 = 0.45,
    #voteshare_slope_down :: Float64 = 0.5,
    #d_col                :: String  = "G24PREDHAR",
    #r_col                :: String  = "G24PRERTRU"

    # perim_dit :: Dict{SimpleEdge, Float64}

    
    # ATTICUS PARTY FUNCTION 
    # energy_fn = make_party_combined_energy(5, county_ids, 0.015, 9, 380; df, k=9, beta_voteshare=0.45)
    println(energy_fn)
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




function prepare_warm_start()
    data  = JSON.parsefile("LA/LA_processed_precincts_Julia.json")
    nodes = data["nodes"]
    # links = data["edges"]
    adjacency = data["adjacency"]

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
    # for link in links
        # u, v = link["source"], link["target"]
        # add_edge!(g, simple_edge(u + 1, v + 1))
        # perim_dict[simple_edge(u + 1, v + 1)] = link["shared_perim"]
    # end
    for i in 1:length(nodes)
        u = i
        for nbr in adjacency[i]
            v = nbr["id"] + 1
            e = simple_edge(u, v)
            add_edge!(g, e)
            perim_dict[e] = nbr["shared_perim"]
        end
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

    # county_ids = df[!, "COUNTYFP"]
    
    ##### now begins the new business #####

    seed_1 = JSON.parsefile("LA/seed_plan1.json") # needed for seed n!!
    # louisiana_plans = JSON.parsefile("LA/LA_processed_precincts_Julia.json")

    # con2026 = df[!, "CON"] # need this for CON PLAN
    # start from the enacted plan for batched run warm, rather than the seeds
    # unique(con2026)
    # ntd_to_tree
    # markov edge
    # ntd_to_tree ... find thoroughly

    seed1_int_int = Dict{Int, Int}(parse(Int, k) => v for (k, v) in pairs(seed_1)) # needed for seed n!!

    seed_1_ntd = [seed1_int_int[i - 1] + 1 for i in 1:length(nodes)] # needed for seed n!!

    # GEOIDs = [node["GEOID"] for node in nodes]
    # seed_1_ntd = [seed_1[GEOIDs[i]] + 1 for i in 1:length(nodes)]
 

    districts = [[i for i in 1:length(seed_1_ntd) if seed_1_ntd[i] == d] for d in unique(seed_1_ntd)] # needed for seed n!!
    #districts = [[i for i in 1:length(con2026) if con2026[i] == d] for d in unique(con2026)]


    t, m = BeanoInit.partition_to_tree_marked_edges(g, districts)

    return t, m
end

