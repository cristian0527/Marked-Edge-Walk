
######################### initial seed was population imbalanced. population tolerance changed to 5%, unstuck chain! ##########################



include("link_cut_MEW/lct_run_TN.jl")



init = prepare_warm_start()

part = BeanoInit.partition(init[1], init[2])


pops = BeanoInit.tally(df, "population", part)

BeanoInit.within_percent_of_ideal(pops, sum(pops)/9, 0.03) ####initial seed has population inbalance at epsilon ~ 0.2


a, b, c, d = main(; initialization=init)




partitions       = BeanoInit.replay_partitions(a, b)
final_from_tree  = BeanoInit.partition(c, d)
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


using Plots
cs = length.(BeanoInit.cut_edges.(partitions, Ref(g)))
plot(cs)

cty_splits = county_splits.(Ref(g), partitions, Ref(county_ids))
plot(cty_splits)


