using Serialization, Shapefile, DataFrames, Distributions, JSON, ProgressBars
using StatsBase, Graphs, SpecialFunctions
using Distributed, SharedArrays, JLD2

"""
using Pkg
Pkg.activate(".")     # activate the folder containing this Project.toml
Pkg.instantiate()      # installs everything already listed in Project.toml
Pkg.add(["JLD2", "Distributions"])
""" 

addprocs(min(40, Sys.CPU_THREADS - 1))

@everywhere begin
    using Serialization, Graphs, StatsBase, DataFrames, Shapefile, JSON, ProgressBars, Distributions, SpecialFunctions
    include("Marked_edges/beano2.2_WI.jl")
end

@everywhere g_ref  = Ref{Any}(nothing)
@everywhere df_ref = Ref{Any}(nothing)
@everywhere county_ids_ref = Ref{Any}(nothing)

data  = JSON.parsefile("TN/TN_Processed_Precincts_w_CON22.json")
nodes = data["nodes"]
# links = data["links"]
adjacency = data["adjacency"]
g = Graphs.SimpleGraph(length(nodes))


# println(data["adjacency"])
# for link in links
#    u, v = link["source"], link["target"]
#    add_edge!(g, simple_edge(u + 1, v + 1))
# end
for i in 1:length(nodes)
    u = i
    for nbr in adjacency[i]
        v = nbr["id"] + 1
        e = simple_edge(u, v)
        add_edge!(g, e)
        # perim_dict[e] = nbr["shared_perim"]
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
for pid in workers()
    remotecall_fetch(pid, g, df, county_ids) do g_local, df_local, county_ids_local
        g_ref[]          = g_local
        df_ref[]         = df_local
        county_ids_ref[] = county_ids_local
    end
end
@everywhere function sorted_percents_sc(df, ptition)
    tot = tally(df, "G24PREDHAR", ptition) + tally(df, "G24PRERTRU", ptition)
    dem = tally(df, "G24PREDHAR", ptition)
    percs = dem ./ tot
    return sort!(percs)
end

@everywhere function county_splits(ntd, county_ids)
    splits = 0
    for cty in Set(county_ids)
        cty_nodes = findall(county_ids .== cty)
        cty_districts = Set(ntd[cty_nodes])
        if length(cty_districts) > 1
            splits += 1
        end
    end
    return splits
end

#fix it / copy from what we have
# the one from the beano.22_WI.jl 

# @everywhere function calculate_cut_edges(partition, g)
#   cuts = []
#    for e in edges(g)
#        u, v = src(e), dst(e)
#        if partition[u] != partition[v]
#            push!(cuts, e)
#        end
#    end
#    return cuts
# end

@everywhere function calculate_cut_edges(partition, g)
    # cuts = Tuple{Int,Int}[]
    cuts = []
    for e in edges(g)
        u, v = src(e), dst(e)
        if partition[u] != partition[v]
            push!(cuts, simple_edge(u, v))
        end
    end
    return cuts
end


@everywhere function process_single_run(run_dir, i)
    df         = df_ref[]
    county_ids = county_ids_ref[]
    n          = nrow(df)
    rows       = Vector{Vector{Float16}}()
    splits     = Vector{Int16}()
    # make a vector for cut edges
    cut_edges = Vector{Int16}()
    # cut_edges = Vector{Vector{Tuple{Int,Int}}}()
    try
        c = open("$(run_dir)/run$(i).jls", "r") do io
            deserialize(io)
        end
        initial_partition = c[1]
        flips             = c[2]
        c                 = nothing
        current_partition = copy(initial_partition)
        ntd_vec = [current_partition[node] for node in 1:n]
        push!(rows, Vector{Float16}(sorted_percents_sc(df, current_partition)))
        push!(splits, county_splits(ntd_vec, county_ids))
        cuts = calculate_cut_edges(current_partition, g_ref[])
        # will this work? just counting edges cut...
        push!(cut_edges, Int16(length(cuts)))
        # something like this for cut edges / more arrays wanted
        for flip in flips
            if haskey(flip, 0)
                push!(rows, Vector{Float16}(sorted_percents_sc(df, current_partition)))
                push!(splits, county_splits(ntd_vec, county_ids))
                push!(cut_edges, Int16(length(cuts)))
                # add this (?)
                continue
            end
            for (node, new_part) in flip
                current_partition[node] = new_part
                ntd_vec[node]           = new_part
            end
            push!(rows, Vector{Float16}(sorted_percents_sc(df, current_partition)))
            push!(splits, county_splits(ntd_vec, county_ids))
            push!(cut_edges, Int16(length(cuts))) 
            # add this (?)
        end
        flips = nothing
        GC.gc(false)
    catch e
        println("Error processing $(run_dir)/run$(i).jls: $e")
        println(stacktrace(catch_backtrace()))
        return nothing
    end
    # return as (niters × ndistricts) matrix — no jagged arrays — plus the
    # per-iteration county-split counts
    return reduce(hcat, rows)', splits, cut_edges
end
# Process one run_dir, write result to a temp file, return filename
function process_one_dir(run_dir, num_runs, chain_idx, tmp_prefix="tmp_perc")
    println("Processing $run_dir ...")
    tasks = collect(1:num_runs)
    results = pmap(tasks; on_error=ex->nothing) do i
        process_single_run(run_dir, i)
    end
    valid = [r for r in results if r !== nothing]
    if isempty(valid)
        println("  WARNING: all runs failed for $run_dir")
        return nothing
    end
    println("  $(length(valid)) / $num_runs runs succeeded")
    perc_mats    = [r[1] for r in valid]
    split_vecs   = [r[2] for r in valid]
    edges_vecs = [r[3] for r in valid]
    # do this if you want to add more arrays to save

    # stack all runs along iterations axis: each mat is (niters × ndistricts)
    # cat along dim 1 gives (niters*num_runs × ndistricts)
    combined = reduce(vcat, perc_mats)       # (total_iters × ndistricts)
    arr = Array{Float16, 2}(undef, size(combined, 2), size(combined, 1))
    arr .= combined'                         # (ndistricts × total_iters)
    splits_arr = reduce(vcat, split_vecs)    # (total_iters,)
    edges_arr = reduce(vcat, edges_vecs)    # (total_iters,)
    tmp_file = "$(tmp_prefix)_chain$(chain_idx).jld2"
    jldsave(tmp_file; arr, splits_arr, edges_arr)
    println("  Written $tmp_file  ($(round(sizeof(arr)/1024^3, digits=3)) GB)")
    results    = nothing
    valid      = nothing
    perc_mats  = nothing
    split_vecs = nothing
    edges_vecs = nothing
    combined   = nothing
    arr        = nothing
    splits_arr = nothing
    edges_arr  = nothing
    GC.gc(true)
    return tmp_file
end
# Concatenate temp files into final output
function compile_tmp_files(tmp_files, out_file, splits_out_file, edges_out_file)
    println("Compiling $(length(tmp_files)) temp files into $out_file ...")
    # Read first to get dims
    ndistricts, niters = jldopen(tmp_files[1], "r") do jf
        size(jf["arr"])
    end
    nchains = length(tmp_files)
    println("Final array: ndistricts=$ndistricts, niters=$niters, nchains=$nchains")
    println("Estimated size: $(round(ndistricts * niters * nchains * 2 / 1024^3, digits=2)) GB (Float16)")
    final_arr        = Array{Float16, 3}(undef, ndistricts, niters, nchains)
    final_splits_arr = Array{Int16, 2}(undef, niters, nchains)
    final_edges_arr  = Array{Int16, 2}(undef, niters, nchains)
    for (ci, tmp_file) in enumerate(tmp_files)
        jldopen(tmp_file, "r") do jf
            final_arr[:, :, ci]     .= jf["arr"]
            final_splits_arr[:, ci] .= jf["splits_arr"]
            final_edges_arr[:, ci]  .= jf["edges_arr"]
        end
        println("  Loaded chain $ci / $nchains")
    end
    jldsave(out_file; perc_arr=final_arr, compress=true)
    println("Saved $out_file")
    jldsave(splits_out_file; county_splits_arr=final_splits_arr, compress=true)
    println("Saved $splits_out_file")
    jldsave(edges_out_file; cut_edges_arr=final_edges_arr, compress=true)
    println("Saved $edges_out_file")
    # do this if you want to add more arrays to save 
    # (!)(!)(!)(!)(!) add another variable for out to file (!)(!)(!)(!)(!)


    # Clean up temp files
    for tmp_file in tmp_files
        rm(tmp_file)
    end
    println("Temp files deleted.")
end
# --- Main ---
run_dirs = ["Cohn_combat_runs_TN/run1"]
num_runs = 100
tmp_files = String[]
for (ci, run_dir) in enumerate(run_dirs)
    tmp = process_one_dir(run_dir, num_runs, ci)
    tmp !== nothing && push!(tmp_files, tmp)
end
compile_tmp_files(tmp_files, "percents_TN1-100runs.jld2", "county_splits_TN1-100runs.jld2", "cut_edges_TN1-100runs.jld2")
exit()

### Loading the final array for analysis

percs_array = load("percents_TN1-100runs.jld2", "perc_arr")
# sorted dem vote shares
# sorted_dem_votes = size(percs_array, 2)

using Plots
# plot(percs_array)

for i in 1:9
    # Slice the i-th list: row i, all columns, first slice
    # Use vec() to flatten the slice into a 1D vector for clean plotting
    list_data = vec(percs_array[i, :, 1])
    
    # Create and display the plot
    p = plot(list_data, title="$i", xlabel="run", ylabel="percs")
    display(p) 
    
    # savefig("plot_list_$i.png")
end

county_array = load("county_splits_TN.jld2", "county_splits_arr")

plot(county_array)


