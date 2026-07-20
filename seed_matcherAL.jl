using JSON
include("link_cut_MEW/lct_run_AL.jl")
include("link_cut_MEW/lct_mew.jl")

data  = JSON.parsefile("AL/AL_processed_precincts_Julia.json")
nodes = data["nodes"]
links = data["links"]
    # println(data["links"])
# adjacency = data["adjacency"]

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
#for i in 1:length(nodes)
#    u = i
#    for nbr in adjacency[i]
#        v = nbr["id"] + 1
#        e = simple_edge(u, v)
#        add_edge!(g, e)
#        perim_dict[e] = nbr["shared_perim"]
#    end
#end

for link in links
    u, v = link["source"], link["target"]
    add_edge!(g, simple_edge(u + 1, v + 1))
       perim_dict[simple_edge(u + 1, v + 1)] = link["shared_perim"]
end
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

seed1 = JSON.parsefile("AL/AL_assignment_2_100000.json")

seed1
# GEOIDs = [node["GEOID"] for node in nodes]
seed1_int_int = Dict{Int, Int}(parse(Int, k) => v for (k, v) in pairs(seed1))
# geoids 

seed_1_ntd = [seed1_int_int[i - 1] + 1 for i in 1:length(nodes)]
# +1 for 1-based indexing in Julia
# 0 to 1 indexing from python to Julia

# maximum(keys(seed1_int_int))

module BeanoInit
    using ProgressBars
    include(joinpath(@__DIR__, "Marked_edges", "beano2.2_WI.jl"))
    # include("Marked_edges/beano2.2_WI.jl") # this if can't find Beano
end

districts = [[i for i in 1:length(seed_1_ntd) if seed_1_ntd[i] == d] for d in unique(seed_1_ntd)]

using Graphs
mini_gs = [induced_subgraph(g, district) for district in districts]

println([length(connected_components(g_[1])) for g_ in mini_gs])


county_splits(g, county_ids, seed_1_ntd)

county_splits(g, seed_1_ntd, county_ids)




at, m = BeanoInit.partition_to_tree_marked_edges(g, districts)
# [1, 1, 1, 1, 1, 1, 1, 1, 1]

main(; initialization = [t, m]) # run the MCMC with the seed plan as initialization

return t, m

# main(; initialization = [t, m]) # run the MCMC with the seed plan as initialization


# post def/n of prepare_warm_start testing #### 

include("link_cut_MEW/lct_run_TN.jl")
initialization = prepare_warm_start()

a, b, c, d = main(; initialization) # run the MCMC with the seed plan as initialization

partitions = BeanoInit.replay_partitions(a, b)
final_part = BeanoInit.replay_partitions(c, d)

@assert partitions[end] == final_part

ntds = [[partition[i] for i in 1:length(partition)] for partition in partitions]
cty_splits = county_splits.county_splits(ntds, county_ids)

# using Plots


# 200 k run on AL for the updated betas and targets 
# --> 


