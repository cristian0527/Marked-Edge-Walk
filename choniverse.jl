

# println(names(df))

county_ids = df[!, "COUNTYFP"]
    
# Set(county_ids)


function county_splits(g, ntd, county_ids)
    splits = 0
    for cty in Set(county_ids) # takes the unique county ids
        cty_nodes = findall(county_ids .== cty)
        cty_districts = Set(ntd[cty_nodes])
        # ntd is .assignment in python
        # 
        if length(cty_districts) > 1
            splits += 1
        end
    
    end
    return splits

end

df[!, "CON"] 

ntd_enacted = df[!, "CON"]

county_splits(g, ntd_enacted, county_ids)
# do this for SLDU and SLDL to see if any errors occur

# K is the numbher of districts
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

methods(make_county_splits_energy)

# el_tigre = make_county_splits_energy(0.4, county_ids)

# el_tigre(g, ntd_sh, ntd_enacted)
#ntd_sh is state house

a, b, c ,d = main() 
partitions = BeanoInit.replay_partitions(a, b)
final_part = BeanoInit.partition(c, d)
@assert partitions[end] == final_part "trajectory replay disagrees with final tree/marked edges"


ntds = [[partition[i] for i in 1:length(partition)] for partition in partitions]

cty_splits = county_splits.(Ref(g), partitions, Ref(county_ids))

# energy function attempts to minimize county splits
# sets up score function e^-b(x)
# give a high score for low number of county splits, and a low score for high number of county splits
# ratio of new thing / old thing (e^-b(x_new - x_old))
# score function into the energy function 
#  

using Plots 

plot(cnty_splits)

unique(ntds)

cs = length.(BeanoInit.cut_edges.(partitions, Ref(g)))
plot(cs)

# next step is to do cut edges and make 

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

        cuts_new = length(BeanoInit.cut_edges(g, new_ntd))
        cuts_old = length(BeanoInit.cut_edges(g, old_ntd))

        return -beta_county_splits * ((cnty_splits_new - cty_target)^2 - (cnty_splits_old - cty_target)^2) - 
        beta_cuts * ((cuts_new - target_cuts)^2 - (cuts_old - target_cuts)^2)
    end



end

# USE THE FUNCTION IN SC LC

a, b, c , d = main(; initialization=[c, d])

partitions = BeanoInit.replay_partitions(a, b)
