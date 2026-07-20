# using Distributed
# addprocs(10)

# @everywhere begin
include("link_cut_MEW/lct_run_MS.jl")

function long_run(n, j)
    run_dir = "Warm_runs_MS/run$(j)"

    if isdir(run_dir)
        existing_files = filter(f -> endswith(f, ".jls"), readdir(run_dir))
        if !isempty(existing_files)
            run_numbers = [parse(Int, match(r"run(\d+)\.jls", f).captures[1]) for f in existing_files]
            last_run    = maximum(run_numbers)
            println("Found existing runs up to run$(last_run) in $(run_dir)")
            last_file      = "$(run_dir)/run$(last_run).jls"
            c_last         = deserialize(last_file)
            initialization = (c_last[3], c_last[4])  # (tree, marked_edges)
            start_i        = last_run + 1
        else
            start_i        = 1
            # c1 = main() need to cahnged this to get initialization from our warm start
            initialization = prepare_warm_start()
            # c is the vector 
            # d is the tree

            # main(; initialization=(c, d))

            #seed1 = JSON.parsefile("TN/seed_plan1.json")
            #GEOIDs = [node["GEOID"] for node in nodes]
            #seed_1_ntd = [seed1[GEOIDs[i]] + 1 for i in 1:length(nodes)]
            #districts = [[i for i in 1:length(seed_1_ntd) if seed_1_ntd[i] == d] for d in unique(seed_1_ntd)]
            


            #initialization = prepare_warm_start()
            # initialization = (c1[3], c1[4])
        end
    else
        println("Creating new directory: $(run_dir)")
        mkpath(run_dir)
        # c1             = main() initialization from our warm start
        # initialization = (c1[3], c1[4])
        initialization = prepare_warm_start()
        start_i        = 1
    end

    if start_i <= n
        for i in start_i:n
            c2             = main(; initialization) # main(; initialization)
            # main(;initialization=[c,d])
            serialize("$(run_dir)/run$(i).jls", c2)
            initialization = (c2[3], c2[4])
            GC.gc()
        end
    else
        println("All $(n) runs already completed for this parameter set.")
    end

    return nothing
end

# function run_wrapper(params)
#    j = params.j
#    println("Starting run $j")
#    try
#        long_run(num_iterations, j)
#    catch e
#        println("Error in run $j: $e")
#        println(stacktrace(catch_backtrace()))
#    end
#    return nothing
# end

num_iterations = 1
# it was 2, maybe because we had to debug
# maybe to 15? 
# bet     = 0.0    # set beta here
# m2      = 0      # set target_cuts here
#epsilon = 0.1   # set population tolerance here

# j_values          = 1:10
# param_combinations = vec([(j=j,) for j in j_values])
# end

# results = pmap(run_wrapper, param_combinations)


long_run(num_iterations, 1)

exit() # vital

# visuals
# main(;initialization=(c,d))


# main(;initialization=(c,d))







# scoreg(P) = e^Bs[ Bcs(xcs new - xcs old) + Bce(cutedges new x - cutedges old x)]
# scale this Bs temperature --> this tells us to make the value more important overtime 
# beta goes up to cool 

# seed one but target 360
# measure cut edges at seeds


