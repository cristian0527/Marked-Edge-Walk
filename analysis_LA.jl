# analysis_TN.jl — mirrors analysis.jl, adapted for TN's adjacency-format dual graph JSON
# fix for Louisina 

include("link_cut_MEW/lct_run_TN.jl")
using Plots

init, flip, final_t, final_m = main()

partitions       = BeanoInit.replay_partitions(init, flip)
final_from_tree  = BeanoInit.partition(final_t, final_m)
@assert partitions[end] == final_from_tree "trajectory replay disagrees with final tree/marked edges"


data  = JSON.parsefile("TN/TN_Processed_Precincts_w_CON22.json")
nodes = data["nodes"]
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
"""
for link in links
    u, v = link["source"], link["target"]
    add_edge!(g, simple_edge(u + 1, v + 1))
    perim_dict[simple_edge(u + 1, v + 1)] = link["shared_perim"]
end
"""

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

cs = length.(BeanoInit.cut_edges.(partitions, Ref(g)))
plot(cs)

cty_splits = county_splits.(Ref(g), partitions, Ref(county_ids))
plot(cty_splits)


a, b, c, d = main(; initialization=[final_t, final_m])

partitions       = BeanoInit.replay_partitions(a, b)
final_from_tree  = BeanoInit.partition(c, d)
@assert partitions[end] == final_from_tree "trajectory replay disagrees with final tree/marked edges"


"""cs = length.(BeanoInit.cut_edges.(partitions, Ref(g)))
plot(cs)
histogram(cs)
"""
"""
pps = BeanoInit.polsby_popper.(Ref(areas), Ref(boundary_lengths), Ref(perim_dict), partitions, Ref(g))

plot(pps)

unique(partitions)
"""
"""
data  = JSON.parsefile("TN/TN_Processed_Precincts_w_CON22.json")
nodes = data["nodes"]
adjacency = data["adjacency"]
enacted_ntd = [node["CON"] for node in nodes]

enacted_cuts = BeanoInit.cut_edges(enacted_ntd, g)
println(length(enacted_cuts))
"""