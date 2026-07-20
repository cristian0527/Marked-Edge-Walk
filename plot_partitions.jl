# plotting plots of the ensemble + partitions
using CairoMakie
using DataFrames
using GeoDataFrames
using GeoMakie

gdf = GeoDataFrames.read("LA/LA_Processed_Precincts_w_CON22_and24.shp") # this is the shp file with the enacted plan (2026) + CON 2022 and CON 2024


function add_partition_column!(
    gdf::AbstractDataFrame,
    partition,
    nodes;
    column::Symbol = :current,
    geography_id::Symbol = :PRECINCTID,
    graph_id::AbstractString = "GEOID",
)
    geography_id in propertynames(gdf) ||
        throw(ArgumentError("geographic data has no $(geography_id) column"))

    #nodes = JSON.parsefile(graph_json)["nodes"]
    assignment = Dict{String, Int}()

    for node in nodes
        # haskey(node, graph_id) ||
            # throw(ArgumentError("graph node has no id"))

        # Graph JSON IDs are zero-based; the sampler's graph vertices are one-based.
        sampler_id = Int(node["id"]) + 1
        haskey(partition, sampler_id) ||
            throw(ArgumentError("partition has no assignment for node $(sampler_id)"))
        assignment[string(node[graph_id])] = Int(partition[sampler_id])
    end

    geographic_ids = string.(gdf[!, geography_id])
    unmatched = unique([id for id in geographic_ids if !haskey(assignment, id)])
    isempty(unmatched) || throw(ArgumentError(
        "$(length(unmatched)) geographic unit(s) have no graph assignment; " *
        "first unmatched ID: $(first(unmatched))",
    ))

    gdf[!, column] = [assignment[id] for id in geographic_ids]
    return gdf
end






add_partition_column!(
        gdf,
        NTD_THING,
        nodes;
        column = :current,
        geography_id = "PRECINCTID",
        graph_id = "GEOID",
    )

districts = sort(unique(gdf[!, :current]))
    district_index = Dict(districts .=> eachindex(districts))
    colors = [district_index[d] for d in gdf[!, :current]]

    fig = Figure(size = (1000, 750))
    ax = Axis(fig[1, 1]; title = "State Name Goes Here", aspect = DataAspect())
    hidespines!(ax)
    hidedecorations!(ax)

    plot = poly!(
        ax,
        gdf.geometry;
        color = colors,
        colormap = :tab20,
        colorrange = (0.5, length(districts) + 0.5),
        strokecolor = (:white, 0.35),
        strokewidth = 0.25,
    )
    Colorbar(
        fig[1, 2],
        plot;
        label = "District",
        ticks = (eachindex(districts), string.(districts)),
    )

display(fig)