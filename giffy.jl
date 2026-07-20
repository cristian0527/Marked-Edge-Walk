import Pkg
Pkg.add(["CairoMakie", "DataFrames", "GeoDataFrames", "GeoMakie"])
using CairoMakie
using DataFrames
using GeoDataFrames
using GeoMakie
using JSON
using Serialization, Graphs

gdf = GeoDataFrames.read("AL/AL_Processed_Precincts_w_CON21_and23.shp")


data  = JSON.parsefile("AL/AL_processed_precincts_Julia.json")  # match your actual AL filename
nodes = data["nodes"]
# '/Users/cristiancastellanos/URSI 2026/Alabama/AL_Processed_Precincts_w_CON21_and23.shp'

# initial_partition, flips, final_tree, marked_edges = c

function final_partition_from_flips(c)
    initial_partition, flips = c[1], c[2]
    current = copy(initial_partition)
    for flip in flips
        for (node, new_part) in flip
            current[node] = new_part
        end
    end
    return current  # Dict{Int,Int}, exactly what add_partition_column! wants
end

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
        haskey(node, graph_id) ||
            throw(ArgumentError("graph node has no id"))

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

c = deserialize("Warm_runs_AL4million/run1/run1.jls")
initial_partition, flips = c[1], c[2]
# nodes[1]["PRECINCTID"] == gdf.PRECINCTID[1]
length(flips)

# 
n_frames  = 1548
step_size = max(1, length(flips) ÷ n_frames)

current = copy(initial_partition)
frame   = 0

for (it, flip) in enumerate(flips)
    for (node, new_part) in flip
        current[node] = new_part
    end

    if it % step_size == 0
        frame += 1

        # add_partition_column! mutates gdf[!, :current] in place each frame
        add_partition_column!(
            gdf, current, nodes;
            column = :current, geography_id = :PRECINCTID, graph_id = "PRECINCTID",
        )

        districts      = sort(unique(gdf[!, :current]))
        district_index = Dict(districts .=> eachindex(districts))
        colors         = [district_index[d] for d in gdf[!, :current]]

        fig = Figure(size = (1000, 750))
        ax  = Axis(fig[1, 1]; title = "AL 4 million OG - Step $it", aspect = DataAspect())
        hidespines!(ax)
        hidedecorations!(ax)

        plot = poly!(
            ax, gdf.geometry;
            color = colors, colormap = :tab20,
            colorrange = (0.5, length(districts) + 0.5),
            strokecolor = (:white, 0.35), strokewidth = 0.25,
        )
        Colorbar(fig[1, 2], plot; label = "District", ticks = (eachindex(districts), string.(districts)))

        save("AL_gif_4millionOG/step$frame.png", fig)

        frame >= n_frames && break
    end
end




#NTD_THING = final_partition_from_flips(c)

#add_partition_column!(
#        gdf,
#        NTD_THING,
#        nodes;
#        column = :current,
#        geography_id = "PRECINCTID",
#        graph_id = "GEOID",
#        )

#add_partition_column!(
#        gdf,
#        NTD_THING,
#        nodes;
#        column = :current,
#        geography_id = :PRECINCTID,
#        graph_id = "PRECINCTID",
#        )


#districts = sort(unique(gdf[!, :current]))
 #   district_index = Dict(districts .=> eachindex(districts))
  #  colors = [district_index[d] for d in gdf[!, :current]]

   # fig = Figure(size = (1000, 750))
    #ax = Axis(fig[1, 1]; title = "AL 4 million walk original", aspect = DataAspect())
    #hidespines!(ax)
    #hidedecorations!(ax)

    #plot = poly!(
    #    ax,
    #    gdf.geometry;
    #    color = colors,
    ##    colormap = :tab20,
     #   colorrange = (0.5, length(districts) + 0.5),
     #   strokecolor = (:white, 0.35),
     #   strokewidth = 0.25,
    #)
    #Colorbar(
    #    fig[1, 2],
    #    plot;
    #    label = "District",
    #    ticks = (eachindex(districts), string.(districts)),
    #


#To actually add the column:
#add_partition_column!(
#        gdf,
#        partition,
#        nodes;
#        column = :current,
#        geography_id = "PRECINCTID",
#        graph_id = "GEOID",
#    )

#giffication
# mkpath("AL_gif_4millionOG")

#for zz in 1:1548
##    fig = Figure(size = (1000, 750))
 #   ax = Axis(fig[1, 1]; title = "AL 4 million OG - Step $zz", aspect = DataAspect())
 #   hidespines!(ax)
  #  hidedecorations!(ax)

   # plot = poly!(
    #   ax,
   #     gdf.geometry;
     #   color = colors,
     #   colormap = :tab20,
     ##   colorrange = (0.5, length(districts) + 0.5),
     #   strokecolor = (:white, 0.35),
     #   strokewidth = 0.25,
    #)
   # Colorbar(
       # fig[1, 2],
     #   plot;
       # label = "District",
      #  ticks = (eachindex(districts), string.(districts)),
   #)

   # save("./AL_gif_4millionOG/step$zz.png",fig)
#end