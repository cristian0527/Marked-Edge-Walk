include("link_cut_MEW/lct_run_SC.jl")
using Plots

init, flip, final_t, final_m = main()

partitions       = BeanoInit.replay_partitions(init, flip)
final_from_tree  = BeanoInit.partition(final_t, final_m)
@assert partitions[end] == final_from_tree "trajectory replay disagrees with final tree/marked edges"


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

cs = length.(BeanoInit.cut_edges.(partitions, Ref(g)))

plot(cs)


a, b, c, d = main(; initialization=[final_t, final_m])

partitions       = BeanoInit.replay_partitions(a, b)
final_from_tree  = BeanoInit.partition(c, d)
@assert partitions[end] == final_from_tree "trajectory replay disagrees with final tree/marked edges"


cs = length.(BeanoInit.cut_edges.(partitions, Ref(g)))
plot(cs)
histogram(cs)

pps = BeanoInit.polsby_popper.(Ref(areas), Ref(boundary_lengths), Ref(perim_dict), partitions, Ref(g))

plot(pps)


unique(partitions)



using Plots

county_ids = df[!, "COUNTYFP"]

cty_splits = county_splits.(Ref(g), partitions, Ref(county_ids))
plot(cty_splits)



a, b, c, d = main(; initialization = [c, d])

partitions = BeanoInit.replay_partitions(a, b)
final_part = BeanoInit.partition(c, d)
@assert partitions[end] == final_part "trajectory replay disagrees with final tree/marked edges"

ntds = [[partition[i] for i in 1:length(partition)] for partition in partitions]

cty_splits = county_splits.(Ref(g), ntds, Ref(county_ids))

plot(cty_splits)

unique(ntds)

cs = length.(BeanoInit.cut_edges.(partitions, Ref(g)))

plot(cs)