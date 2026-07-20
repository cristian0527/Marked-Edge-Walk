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
const K          = 9
const EPSILON    = 0.05 # 0.01 before
# within 5% balance of the population

const N_ITERS    = 10_000

# main function
function main(; initialization=nothing)

    data  = JSON.parsefile("TN/TN_Processed_Precincts_w_CON22.json")
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
    """
    for link in links
        u, v = link["source"], link["target"]
        add_edge!(g, simple_edge(u + 1, v + 1))
        perim_dict[simple_edge(u + 1, v + 1)] = link["shared_perim"]
    end
    """
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

    # swap which energy function drives the chain by changing this line only
    # energy_fn = make_polsby_popper_energy(1000.0, areas, boundary_lengths, perim_dict, K)
    # energy_fn = make_cuts_energy(0.1, 600) 
    # energy_fn = make_county_splits_energy(0.4, df[!, "COUNTYFP"])
    # energy_fn = make_combined_energy(.1, df[!, "COUNTYFP"], 0.015, 10, 500)
    
    
    # make counties more important
    # energy_fn = make_combined_energy(.25, df[!, "COUNTYFP"], 0.015, 10, 500)
    # let's try 0.35
    
    county_ids = df[!, "COUNTYFP"]
    # energy_fn = make_combined_energy(5, df[!, "COUNTYFP"], 0.015, 10, 500) 
    # energy_fn = make_combined_energy(5, county_ids, 0.015, 9, 380)
    # THIS ONE ABOVE WORKS - 4 million JUL 10


    # 4 million walk make_combined_energy(5, county_ids, 0.015, 10, 380)

    # county -> 9 
    # relaxing coefficient a little (?)
    # 380 -> 400
    # 0.1 -> change to 1 
    # changed 1 - 0.8
    # tried 2 ->3.25

    # try 380 for burnin -> let's try this after this 100k walk(?)
    # or budgeting



    # recom runs look back 
    # beta county to zero to check if something is up and then increment slowly up...




    # county beta value that looks like a good distribution (unimodal) 
    # look at recom cut edges 

    # this is in the energy function 

    # interested in nudging the coefficient up and down 
    # above: once we find an ideal coefficient 
    # also for target, shift from up and down 
    # notes above 

    # trying this out
    # doing a version that just tries to minimize the county splits coefficients rather than the guassian 
    # it is possibile: replace cut edge target to minimize cut edges target 
    # modify the new and old 
    # changing the distribution by doing the summary plots and observing
    # do a warm start 
    

    # running 2 ensembles: 1 minimized cut edges / 1 normal
    # minimize: minus times value
    # running 4 ensembles: 2 of each
    
    # see the difference on the same parameters, just minimized
    # energy_fn_minimized = make_combined_energy_county_minimized(0.35, df[!, "COUNTYFP"], 0.015)

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
    energy_fn = make_party_combined_energy(5, county_ids, 0.015, 9, 380; df, k=9, beta_voteshare=100)
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


# ntd_to_tree
# markov edge
"""
# copy paste from seed file 
function prepare_warm_start()

    data  = JSON.parsefile("TN/TN_Processed_Precincts_w_CON22.json")
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
        
        # for link in links
            # u, v = link["source"], link["target"]
            # add_edge!(g, simple_edge(u + 1, v + 1))
            # perim_dict[simple_edge(u + 1, v + 1)] = link["shared_perim"]
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

    county_ids = df[!, "COUNTYFP"]

    ### where business starts ### 

    seed1 = JSON.parsefile("TN/seed_plan1.json")
    GEOIDs = [node["GEOID"] for node in nodes]

    # geoids 

    seed_1_ntd = [seed1[GEOIDs[i]] + 1 for i in 1:length(nodes)] # +1 for 1-based indexing in Julia
    # 0 to 1 indexing from python to Julia

    districts = [[i for i in 1:length(seed_1_ntd) if seed_1_ntd[i] == d] for d in unique(seed_1_ntd)]

    t, m = BeanoInit.partition_to_tree_marked_edges(g, districts)

    return t,m # run the MCMC with the seed plan as initialization 
end 

"""



function prepare_warm_start()
    data  = JSON.parsefile("TN/TN_Processed_Precincts_w_CON22.json")
    nodes = data["nodes"]
    # links = data["links"]
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

    county_ids = df[!, "COUNTYFP"]


    ##### now begins the new business #####

    seed_1 = JSON.parsefile("TN/seed_plan5.json")
    seed1_int_int = Dict{Int, Int}(parse(Int, k) => v for (k, v) in pairs(seed_1))

    seed_1_ntd = [seed1_int_int[i - 1] + 1 for i in 1:length(nodes)]

    # GEOIDs = [node["GEOID"] for node in nodes]
    # seed_1_ntd = [seed_1[GEOIDs[i]] + 1 for i in 1:length(nodes)]

    districts = [[i for i in 1:length(seed_1_ntd) if seed_1_ntd[i] == d] for d in unique(seed_1_ntd)]

    t, m = BeanoInit.partition_to_tree_marked_edges(g, districts)

    return t, m
end

