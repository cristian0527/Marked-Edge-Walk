# voteshare_viz.jl — district-level D voteshares over time (Colorado)
# Assumes `result` is in scope from a completed fool.jl run.
# Run: include("voteshare_viz.jl")

using CSV, JSON, DataFrames, Plots

# ── state config ──────────────────────────────────────────────────────────────

const VIZ_JSON  = "Colorado/co.json"
const VIZ_D_COL = "PRE20D"
const VIZ_R_COL = "PRE20R"

# ── reload df ─────────────────────────────────────────────────────────────────

let
data       = JSON.parsefile(VIZ_JSON)
nodes_json = data["nodes"]
df_dict    = Dict{String, Vector{Any}}()
for node in nodes_json
    for (key, value) in node
        key == "id" && continue
        v = get!(df_dict, key, Any[])
        push!(v, value)
    end
end
max_len = maximum(length(v) for v in values(df_dict))
for vals in values(df_dict)
    while length(vals) < max_len; push!(vals, missing); end
end
df_vis = DataFrame(df_dict)
insertcols!(df_vis, 1, :id => 1:nrow(df_vis))
println("Loaded df: $(nrow(df_vis)) precincts")

# ── rebuild initial ntd from result.partition ─────────────────────────────────

ntd = zeros(Int, nrow(df_vis))
for (node, dist) in result.partition
    ntd[node] = dist
end

# ── helper ────────────────────────────────────────────────────────────────────

function dem_voteshares_vec(df, ntd, k)
    d  = tally(df, VIZ_D_COL, ntd, k)
    r  = tally(df, VIZ_R_COL, ntd, k)
    vs = Vector{Float64}(undef, k)
    for i in 1:k
        tot   = d[i] + r[i]
        vs[i] = tot > 0 ? d[i] / tot : 0.5
    end
    return vs
end

# ── replay flips, snapshot voteshares at each trace step ──────────────────────

trace_steps = result.trace.step
n_trace     = length(trace_steps)
voteshares  = Matrix{Float64}(undef, n_trace, K)

best_idx   = argmax(result.trace.metric_value)
maxdem_idx = argmax(result.trace.dem_seats)
ntd_at_best   = zeros(Int, nrow(df_vis))
ntd_at_maxdem = zeros(Int, nrow(df_vis))

trace_idx = 1
if trace_steps[1] == 0
    voteshares[1, :] = dem_voteshares_vec(df_vis, ntd, K)
    1 == best_idx   && (ntd_at_best   = copy(ntd))
    1 == maxdem_idx && (ntd_at_maxdem = copy(ntd))
    trace_idx = 2
end

println("Replaying $(length(result.flips)) flips → $(n_trace) snapshots...")
for step in 1:length(result.flips)
    for (node, new_dist) in result.flips[step]
        ntd[node] = new_dist
    end
    if trace_idx <= n_trace && trace_steps[trace_idx] == step
        voteshares[trace_idx, :] = dem_voteshares_vec(df_vis, ntd, K)
        trace_idx == best_idx   && (ntd_at_best   = copy(ntd))
        trace_idx == maxdem_idx && (ntd_at_maxdem = copy(ntd))
        trace_idx += 1
    end
end
println("Done.")

# ── summary ───────────────────────────────────────────────────────────────────

best_step = trace_steps[best_idx]
best_vs   = voteshares[best_idx, :]
println("Best metric $(round(result.trace.metric_value[best_idx], digits=3)) at step $best_step")
println("  Dem seats at best: $(count(>(0.5), best_vs)) / $K")

maxdem_step = trace_steps[maxdem_idx]
maxdem_vs   = voteshares[maxdem_idx, :]
println("Most Dem districts $(result.trace.dem_seats[maxdem_idx]) at step $maxdem_step")

# ── export partition CSVs ─────────────────────────────────────────────────────

CSV.write("partition_best_metric.csv",
    DataFrame(node_id = 1:nrow(df_vis), district = ntd_at_best))
CSV.write("partition_max_dem.csv",
    DataFrame(node_id = 1:nrow(df_vis), district = ntd_at_maxdem))
println("Partitions saved → partition_best_metric.csv, partition_max_dem.csv")

# ── export sorted voteshares matrix (read by viz.py) ─────────────────────────

sorted_vs_time = mapslices(sort, voteshares; dims = 2)   # (n_trace × K), ascending
df_svs = DataFrame(sorted_vs_time, ["d$i" for i in 1:K])
insertcols!(df_svs, 1, :step => trace_steps)
CSV.write("sorted_voteshares.csv", df_svs)
println("Sorted voteshares saved → sorted_voteshares.csv")

# ── spaghetti plot ────────────────────────────────────────────────────────────

colors = [sorted_vs_time[best_idx, j] > 0.5 ? :blue : :red for j in 1:K]

spaghetti = plot(
    trace_steps, sorted_vs_time;
    color   = reshape(colors, 1, K),
    alpha   = 0.65,
    lw      = 1.2,
    xlabel  = "Step",
    ylabel  = "D two-party voteshare",
    title   = "Colorado — Sorted district D voteshares over time (rank traces)",
    legend  = false,
    size    = (900, 450),
)
hline!(spaghetti, [0.5]; color = :black, lw = 1.5, ls = :dash)
vline!(spaghetti, [best_step];   color = :darkgreen, lw = 2, label = "Best metric")
vline!(spaghetti, [maxdem_step]; color = :orange,    lw = 2, label = "Most Dem seats")

# ── bar chart at the best-scoring step ────────────────────────────────────────

sorted_best = sort(best_vs)
bar_colors  = [vs > 0.5 ? :blue : :red for vs in sorted_best]

best_bar = bar(
    1:K, sorted_best;
    color   = bar_colors,
    xlabel  = "District (sorted by D share)",
    ylabel  = "D two-party voteshare",
    title   = "Best step ($best_step): $(count(>(0.5), best_vs)) Dem / $(count(<=(0.5), best_vs)) Rep districts",
    legend  = false,
    xticks  = false,
    ylims   = (0.0, 1.0),
)
hline!(best_bar, [0.5]; color = :black, lw = 1.5, ls = :dash)

# ── bar chart at the most-democratic step ─────────────────────────────────────

sorted_maxdem = sort(maxdem_vs)
maxdem_colors = [vs > 0.5 ? :blue : :red for vs in sorted_maxdem]

maxdem_bar = bar(
    1:K, sorted_maxdem;
    color   = maxdem_colors,
    xlabel  = "District (sorted by D share)",
    ylabel  = "D two-party voteshare",
    title   = "Most Dem step ($maxdem_step): $(result.trace.dem_seats[maxdem_idx]) Dem / $(K - result.trace.dem_seats[maxdem_idx]) Rep districts",
    legend  = false,
    xticks  = false,
    ylims   = (0.0, 1.0),
)
hline!(maxdem_bar, [0.5]; color = :black, lw = 1.5, ls = :dash)

# ── combined figure ────────────────────────────────────────────────────────────

p = plot(spaghetti, best_bar, maxdem_bar; layout = (3, 1), size = (900, 1150))
display(p)
savefig(p, "voteshares_over_time.png")
println("Saved → voteshares_over_time.png")
end  # let