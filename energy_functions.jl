# energy_functions.jl
# Extracted from choniverse.jl: only the reusable function definitions,
# with no top-level script code that depends on df/g/partitions existing yet.
# Safe to `include` before calling main().

function county_splits(g, ntd, county_ids)
    splits = 0
    for cty in Set(county_ids)
        cty_nodes = findall(county_ids .== cty)
        cty_districts = Set(ntd[i] for i in cty_nodes)
        if length(cty_districts) > 1
            splits += 1
        end
    end
    return splits
end



# K is the number of districts
function make_county_splits_energy(
    beta :: Float64,
    county_ids :: Vector{Any}
)
    return function (g, new_ntd, old_ntd, county_ids, beta)
        cnty_splits_new = county_splits(g, new_ntd, county_ids)
        cnty_splits_old = county_splits(g, old_ntd, county_ids)
        return -beta * (cnty_splits_new - cnty_splits_old)
    end
end

function make_combined_energy(
    beta_county_splits :: Number,
    county_ids :: Vector{Any},
    beta_cuts :: Number,
    cty_target :: Number,
    target_cuts :: Number
    # perim_dit :: Dict{SimpleEdge, Float64}
)
    return function combined_energy(g, new_ntd, old_ntd)
        cnty_splits_new = county_splits(g, new_ntd, county_ids)
        cnty_splits_old = county_splits(g, old_ntd, county_ids)

        cuts_new = length(BeanoInit.cut_edges(new_ntd, g))
        cuts_old = length(BeanoInit.cut_edges(old_ntd, g))

        # minimize: county splits minus times value 
        # to minimize the cut edges
        
        return -beta_county_splits * ((cnty_splits_new - cty_target)^2 - (cnty_splits_old - cty_target)^2) -
        beta_cuts * ((cuts_new - target_cuts)^2 - (cuts_old - target_cuts)^2)
    end
end



# county splits minimized 
# delete this 
function make_county_splits_energy_minimized(
    beta :: Float64,
    county_ids :: Vector{Any}
)
    return function (g, new_ntd, old_ntd, county_ids, beta)
        cnty_splits_new = county_splits(g, new_ntd, county_ids)
        cnty_splits_old = county_splits(g, old_ntd, county_ids)
        return -beta * (cnty_splits_new - cnty_splits_old)
    end
end


function make_combined_energy_county_minimized(
    beta_county_splits :: Number,
    county_ids :: Vector{Any},
    beta_cuts :: Number,
    cty_target :: Number,
    target_cuts :: Number
    # perim_dit :: Dict{SimpleEdge, Float64}
    
)
    return function combined_energy_county_minimized(g, new_ntd, old_ntd)
        cnty_splits_new = county_splits(g, new_ntd, county_ids)
        cnty_splits_old = county_splits(g, old_ntd, county_ids)

        cuts_new = length(BeanoInit.cut_edges(new_ntd, g))
        cuts_old = length(BeanoInit.cut_edges(old_ntd, g))

        # minimize: county splits minus times value 
        # to minimize the cut edges
        return -beta_county_splits * (cnty_splits_new - cnty_splits_old) -
        beta_cuts * ((cuts_new - target_cuts)^2 - (cuts_old - target_cuts)^2)
        # - beta so we can accept it

    end
end

# for the beta value



# ntd_con = df[!, "CON"]

# districts = [[i for i in 1:length(ntd_con) if ntd_con[i] == d] for d in unique(ntd_con)]

# t, m = BeanoInit.partition_to_tree_marked_edges(g, districts)

# @Angela Moon this should give you a tree, marked edges that you can start a run with: 
    #for district in districts
    #    g_sub = induced_subgraph(g, district)
    #    if !isconnected(g_sub)
    #        println("District $district is not connected.")
    #    end
    #end 
# main(; initialization = [t, m])






