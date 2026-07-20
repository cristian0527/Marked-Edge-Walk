using Serialization, Shapefile, DataFrames, Distributions, JSON, ProgressBars
using StatsBase, Graphs, SpecialFunctions
using Distributed, SharedArrays, JLD2

addprocs(min(40, Sys.CPU_THREADS - 1))

@everywhere begin
    using Serialization, Graphs, StatsBase, DataFrames, Shapefile, JSON, ProgressBars, Distributions, SpecialFunctions
    include("Marked_edges/beano2.2_WI.jl")
end

@everywhere g_ref          = Ref{Any}(nothing)
@everywhere df_ref         = Ref{Any}(nothing)
@everywhere county_ids_ref = Ref{Any}(nothing)

data  = JSON.parsefile("TN/TN_Processed_Precincts_w_CON22.json")
# data  = JSON.parsefile("TN/seed_plan1.json")

nodes = data["nodes"]
# links = data["links"]
adjacency = data["adjacency"]
g     = Graphs.SimpleGraph(length(nodes))
#for link in links
#    u, v = link["source"], link["target"]
#    add_edge!(g, simple_edge(u + 1, v + 1))
#end

for i in 1:length(nodes)
    u = i
    for nbr in adjacency[i]
        v = nbr["id"] + 1
        e = simple_edge(u, v)
        add_edge!(g, e)
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

@everywhere function sorted_percents_black(df, ptition)
    tot_black = tally(df, "Total", ptition) # + tally(df, "G24PRERTRU", ptition)
    black = tally(df, "NH_Black", ptition)

    percs_black = black ./ tot_black

    return sort!(percs_black)
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

@everywhere function process_single_run(run_dir, i)
    g          = g_ref[]
    df         = df_ref[]
    county_ids = county_ids_ref[]
    n          = nrow(df)
    rows       = Vector{Vector{Float16}}()
    splits     = Vector{Int16}()
    cuts       = Vector{Int32}()
    black_rows = Vector{Vector{Float16}}()

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
        push!(cuts, Int32(length(cut_edges(current_partition, g))))
        push!(black_rows, Vector{Float16}(sorted_percents_black(df, current_partition)))

        for flip in flips
            if haskey(flip, 0)
                push!(rows, Vector{Float16}(sorted_percents_sc(df, current_partition)))
                push!(splits, county_splits(ntd_vec, county_ids))
                push!(cuts, Int32(length(cut_edges(current_partition, g))))
                push!(black_rows, Vector{Float16}(sorted_percents_black(df, current_partition)))
                continue
            end
            for (node, new_part) in flip
                current_partition[node] = new_part
                ntd_vec[node]           = new_part
            end
            push!(rows, Vector{Float16}(sorted_percents_sc(df, current_partition)))
            push!(splits, county_splits(ntd_vec, county_ids))
            push!(cuts, Int32(length(cut_edges(current_partition, g))))
            push!(black_rows, Vector{Float16}(sorted_percents_black(df, current_partition)))
        end

        flips = nothing
        GC.gc(false)

    catch e
        println("Error processing $(run_dir)/run$(i).jls: $e")
        println(stacktrace(catch_backtrace()))
        return nothing
    end

    # return as (niters × ndistricts) matrix — no jagged arrays — plus the
    # per-iteration county-split and cut-edge counts
    return reduce(hcat, rows)', splits, cuts, reduce(hcat, black_rows)'
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

    perc_mats  = [r[1] for r in valid]
    split_vecs = [r[2] for r in valid]
    cut_vecs   = [r[3] for r in valid]
    black_perc_mats = [r[4] for r in valid]

    # stack all runs along iterations axis: each mat is (niters × ndistricts)
    # cat along dim 1 gives (niters*num_runs × ndistricts)
    combined = reduce(vcat, perc_mats)       # (total_iters × ndistricts)
    arr = Array{Float16, 2}(undef, size(combined, 2), size(combined, 1))
    arr .= combined'                         # (ndistricts × total_iters)

    splits_arr = reduce(vcat, split_vecs)    # (total_iters,)
    cuts_arr   = reduce(vcat, cut_vecs)      # (total_iters,)

    # black percents 
    combined_black = reduce(vcat, black_perc_mats)       # (total_iters × ndistricts)
    arr_black = Array{Float16, 2}(undef, size(combined_black, 2), size(combined_black, 1))
    arr_black .= combined_black'                         # (ndistricts × total_iters)


    tmp_file = "$(tmp_prefix)_chain$(chain_idx).jld2"
    jldsave(tmp_file; arr, splits_arr, cuts_arr, arr_black)
    println("  Written $tmp_file  ($(round(sizeof(arr)/1024^3, digits=3)) GB)")

    results    = nothing
    valid      = nothing
    perc_mats  = nothing
    split_vecs = nothing
    cut_vecs   = nothing
    black_perc_mats = nothing
    combined   = nothing
    arr        = nothing
    splits_arr = nothing
    cuts_arr   = nothing
    combined_black = nothing
    arr_black = nothing
    GC.gc(true)

    return tmp_file
end

# Concatenate temp files into final output
function compile_tmp_files(tmp_files, out_file, splits_out_file, cuts_out_file, black_out_file)
    println("Compiling $(length(tmp_files)) temp files into $out_file ...")

    # Read first to get dims
    ndistricts, niters = jldopen(tmp_files[1], "r") do jf
        size(jf["arr"])
    end
    nchains = length(tmp_files)
     # there is an error with this line  # endndistricts ? maybe i accidentally copied down the variable / pasted next to end i think?
    # how did this error come? 


    println("Final array: ndistricts=$ndistricts, niters=$niters, nchains=$nchains")
    println("Estimated size: $(round(ndistricts * niters * nchains * 2 / 1024^3, digits=2)) GB (Float16)")

    # RESUME BACK HERE
    final_arr        = Array{Float16, 3}(undef, ndistricts, niters, nchains)
    final_splits_arr = Array{Int16, 2}(undef, niters, nchains)
    final_cuts_arr   = Array{Int32, 2}(undef, niters, nchains)
    final_arr_black  = Array{Float16, 3}(undef, ndistricts, niters, nchains)


    for (ci, tmp_file) in enumerate(tmp_files)
        jldopen(tmp_file, "r") do jf
            final_arr[:, :, ci]     .= jf["arr"]
            final_splits_arr[:, ci] .= jf["splits_arr"]
            final_cuts_arr[:, ci]   .= jf["cuts_arr"]
            final_arr_black[:, :, ci]     .= jf["arr_black"]

        end
        println("  Loaded chain $ci / $nchains")
    end

    jldsave(out_file; perc_arr=final_arr, compress=true)
    println("Saved $out_file")

    jldsave(splits_out_file; county_splits_arr=final_splits_arr, compress=true)
    println("Saved $splits_out_file")

    jldsave(cuts_out_file; cut_edges_arr=final_cuts_arr, compress=true)
    println("Saved $cuts_out_file")

    jldsave(black_out_file; blackperc_arr=final_arr_black, compress=true)
    println("Saved $black_out_file")

    # Clean up temp files
    for tmp_file in tmp_files
        rm(tmp_file)
    end
    println("Temp files deleted.")
end


# --- Main ---
run_dirs = ["Warm_runs_TN/run1"]
num_runs = 1

tmp_files = String[]
for (ci, run_dir) in enumerate(run_dirs)
    tmp = process_one_dir(run_dir, num_runs, ci)
    tmp !== nothing && push!(tmp_files, tmp)
end

# this should now work? 
compile_tmp_files(tmp_files, "warm_runs_percents_TN.jld2", "warm_runs_county_splits_TN.jld2", "warm_runs_cut_edges_TN.jld2", "warm_runs_black_percs_TN.jld2")
exit() # push_reader_claude_w_cuts

# 23







# visualization of county splits and cut edges
using Plots, JLD2

county_array_15steps = load("warm_runs_county_splits_TN.jld2", "county_splits_arr")
plot(county_array_15steps)


cut_array_15steps = load("warm_runs_cut_edges_TN.jld2", "cut_edges_arr")
plot(cut_array_15steps)



# PLOTS FOR THE NORMAL WARM 
# using Plot 

# percs_array = load("warm_runs_percents_TN.jld2", "perc_arr")
# sorted dem vote shares
# sorted_dem_votes = size(percs_array, 2)
# plot(percs_array)

#for i in 1:9
    # Slice the i-th list: row i, all columns, first slice
    # Use vec() to flatten the slice into a 1D vector for clean plotting
#    list_data = vec(percs_array[i, :, 1])
    
    # Create and display the plot
#    p = plot(list_data, title="$i", xlabel="run", ylabel="percs")
#    display(p) 
    
    # savefig("plot_list_$i.png")
#end

#p = plot(
#    vec(percs_array[1, :, 1]),
#    label="1",
#    xlabel="run",
#    ylabel="percs",
#    title="dem percents"
#)

# Add the rest
#for i in 1:9
#    list_data = vec(percs_array[i, :, 1])
#    plot!(p, list_data, label="$i")
#end

#display(p)


p = plot(title="% dem over run", xlabel="step", ylabel="% dem", legend=true)

for i in 1:9
    for j in 1:size(percs_array, 3)
        plot!(p, percs_array[i, :, j])
    end
end

display(p)




import Pkg
Pkg.add("StatsPlots")
using StatsPlots
# boplots 

percs_array = load("warm_runs_percents_TN.jld2", "perc_arr")

enacted = [0.20767780842194122, 0.3189138681245119, 0.32882478205896054, 0.361378679747658, 0.3782771292833578, 
                      0.3851191451702776, 0.3897077049211585, 0.38999520567046975, 0.399749312865921]





# note: running for 4 million steps

# plot the cs, ce, and the democractic percents in about 3-4 hours

# using StatsPlots

b = plot(title="% Dem Box Plots TN")

for i in 1:9
    boxplot!(b, [i], percs_array[i,:,1], label="", legend=false, outliers=false)
end

# stars aligned correctly
scatter!(b, 1:9, enacted,
         marker=:star5,
         markersize=8,
         color=:red,
         label="Enacted")

# add district labels
xticks!(b, 1:9, ["D$i" for i in 1:9])

display(b)



# black percents boxplot

black_percs_array = load("warm_runs_black_percs_TN.jld2", "blackperc_arr")
enacted_black = [0.030920502743701426, 0.07068214330038997, 0.10842391194669013, 
                       0.1123127079937296, 0.15509604891543902, 0.16308640810613034, 
                       0.28144667592498385, 0.2949367972807872, 0.31560996779663003]


c = plot(title="% Black Box Plots TN")

for i in 1:9
    boxplot!(c, [i], black_percs_array[i,:,1], label="", legend=false, outliers=false)
end

# stars aligned correctly
scatter!(c, 1:9, enacted_black,
         marker=:star5,
         markersize=8,
         color=:red,
         label="Enacted")

# add district labels
xticks!(c, 1:9, ["D$i" for i in 1:9])

display(c)












# Atticus next step:

# OLD SETUP 
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

# function _district_score(share::Float64;
#                         threshold::Float64=0.501, # 0.55
#                         slope_down::Float64=0.5)
#    share <= threshold ? share / threshold : 1.0 - slope_down * (share - threshold)
# end


"""
Only scores the `n` most-Republican districts, each on a standard
_district_score tent peaking at `threshold` (default 55% R). The other
k-n districts are ignored entirely.
Useful when you only care about packing/protecting a fixed number of
safe-R seats and don't want to penalize or reward the rest.
"""
#function rep_voteshare_score_topn(df::DataFrame, node_to_dist::Vector{Int64}, k::Int64, n::Int64;
#                                  d_col::String="PRE20D", r_col::String="PRE20R",
#                                  threshold::Float64=0.55, slope_down::Float64=0.5)
#   d = tally(df, d_col, node_to_dist, k)
#   r = tally(df, r_col, node_to_dist, k)
#   shares = Float64[]
#    for i in 1:k
#        tot = d[i] + r[i]
#       tot == 0 && continue
#        push!(shares, r[i] / tot)
#    end
#    sort!(shares; rev=false) # rev = true
#    score = 0.0
#    for i in 1:min(n, length(shares))
#        score += _district_score(shares[i]; threshold, slope_down)
#    end
#    return score
# end


# -------------------------- NEW SETUP --------------------------

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
                                  d_col::String="PRE20D", r_col::String="PRE20R",
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



# partisian terms in function 
# under 41% TN
#
