function rep_voteshare_score_vector(df::DataFrame, node_to_dist::Vector{Int64}, k::Int64, targets::Vector{Float64}, uses::Vector{Symbol}, slope_down::Number=0.5;
                                  d_col::String="G24PREDBID", r_col::String="G20PRERTRU")
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

    for i in 1:length(shares)
        if uses[i] == :do_nothing
            continue
        elseif uses[i] == :equal
            score += shares[i] <= targets[i] ? shares[i] / targets[i] : 1.0 - slope_down * (shares[i] - targets[i])
        elseif uses[i] == :less
            score += (1-shares[i]) <= (1-targets[i]) ? (1-shares[i]) / (1-targets[i]) : 1.0
        elseif uses[i] ==:greater
            score += shares[i] <= targets[i] ? shares[i] / targets[i] : 1.0
        end
    end
    return score
end

function make_combined_super_party_energy(
    county_ids          :: Vector{Any},
    beta_county_splits :: Number,
    beta_cuts           :: Number,
    cty_target          :: Number,
    target_cuts         :: Number;
    df                   = nothing,
    k                    :: Int64   = 0,
    targets             :: Vector{Float64},
    uses                :: Vector{Symbol},
    beta_voteshare      :: Number  = 0.0,
    voteshare_slope_down :: Number = 0.5,
    d_col                :: String  = "G20PREDBID",
    r_col                :: String  = "G20PRERTRU",
)
    if beta_voteshare != 0.0
        @assert !isnothing(df) && k > 0 "make_combined_energy: df and k are required when beta_voteshare != 0"
    end
    return function combined_energy(g, new_ntd, old_ntd)
        cty_splits_new = county_splits(g, new_ntd, county_ids)
        cty_splits_old = county_splits(g, old_ntd, county_ids)
        cuts_new = count(e -> new_ntd[src(e)] != new_ntd[dst(e)], edges(g))
        cuts_old = count(e -> old_ntd[src(e)] != old_ntd[dst(e)], edges(g))
        energy = -beta_county_splits * ((cty_splits_new - cty_target)^2 - (cty_splits_old - cty_target)^2) -
                 beta_cuts          * ((cuts_new - target_cuts)^2 - (cuts_old - target_cuts)^2)
        if beta_voteshare != 0.0
            score_new = rep_voteshare_score_vector(df, new_ntd, k, targets, uses, voteshare_slope_down;
                                                  d_col=d_col, r_col=r_col)
            score_old = rep_voteshare_score_vector(df, old_ntd, k, targets, uses, voteshare_slope_down;
                                                  d_col=d_col, r_col=r_col)
            energy += beta_voteshare * (score_new - score_old)
        end
        return energy
    end
end

function rep_voteshare_score_vector_gaussian(df::DataFrame, node_to_dist::Vector{Int64}, k::Int64, targets::Vector{Float64}, uses::Vector{Symbol}, slope_down::Number=0.5;
                                  d_col::String="PRE20D", r_col::String="PRE20R")
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

    for i in 1:length(shares)
        if uses[i] == :do_nothing
            continue
        elseif uses[i] == :equal
            score += shares[i] <= targets[i] ? shares[i] / targets[i] : 1.0 - slope_down * (shares[i] - targets[i])^2
        elseif uses[i] == :less
            continue
            # score += (1-shares[i]) <= (1-targets[i]) ? (1-shares[i]) / (1-targets[i]) : 1.0
        elseif uses[i] ==:greater
            continue
            # score += shares[i] <= targets[i] ? shares[i] / targets[i] : 1.0
        end
    end
    return score
end

function make_combined_gaussian_party_energy(
    county_ids          :: Vector{Any},
    beta_county_splits :: Number,
    beta_cuts           :: Number,
    cty_target          :: Number,
    target_cuts         :: Number;
    df                   = nothing,
    k                    :: Int64   = 0,
    targets             :: Vector{Float64},
    uses                :: Vector{Symbol},
    beta_voteshare      :: Number  = 0.0,
    voteshare_slope_down :: Number = 0.5,
    d_col                :: String  = "G20PREDBID",
    r_col                :: String  = "G20PRERTRU",
)
    if beta_voteshare != 0.0
        @assert !isnothing(df) && k > 0 "make_combined_energy: df and k are required when beta_voteshare != 0"
    end
    return function combined_energy(g, new_ntd, old_ntd)
        cty_splits_new = county_splits(g, new_ntd, county_ids)
        cty_splits_old = county_splits(g, old_ntd, county_ids)
        cuts_new = count(e -> new_ntd[src(e)] != new_ntd[dst(e)], edges(g))
        cuts_old = count(e -> old_ntd[src(e)] != old_ntd[dst(e)], edges(g))
        energy = -beta_county_splits * ((cty_splits_new - cty_target)^2 - (cty_splits_old - cty_target)^2) -
                 beta_cuts          * ((cuts_new - target_cuts)^2 - (cuts_old - target_cuts)^2)
        if beta_voteshare != 0.0
            score_new = rep_voteshare_score_vector_gaussian(df, new_ntd, k, targets, uses, voteshare_slope_down;
                                                  d_col=d_col, r_col=r_col)
            score_old = rep_voteshare_score_vector_gaussian(df, old_ntd, k, targets, uses, voteshare_slope_down;
                                                  d_col=d_col, r_col=r_col)
            energy -= beta_voteshare * (score_new - score_old)
        end
        return energy
    end
end
