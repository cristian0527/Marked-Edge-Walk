# The Marked Edge Walk (MEW) for a Statistical Analysis of Minority Representation after Callais

This repository is the keepsake of Professor Daryl DeFord's Undergraduate Research Summer Institute group at Vassar College, which endeavors to analyze the future of Southern States post Louisiana v. Callais.

First of all, we would like to thank Atticus McWhorter for joining us in our summer research journey. Atticus, alongside Professor DeFord, guided us in implementing the Marked Edge Walk for our analysis and taught us the fundamentals in Julia.

### Overview

In late April 2026, the Supreme Court decided a Louisiana congressional map was unconstitutional due to race being a predominant factor in creating the districts. This case, Louisiana v. Callais, significantly reduces protections for minority representation through the Voting Rights Act. Following this case, there have been questions, concerns, and speculation surrounding how this case will impact minority representation throughout the case. One recent New York Times article written by Nate Cohn claimed that the Voting Rights Act did not create additional representation for minority voters and that race blind redistricting would maintain this representation.

This repository analyzes what redistricting would look like in those eight states Cohn argued in his NYT article through the Marked Edge Walk. These eight states are Alabama, Georgia, Florida, Louisiana, Mississippi, North Carolina, South Carolina, and Tennessee. In our analysis, we find sampling distributions that are similar to the current enacted congressional plan in each state, and then generate large ensembles of redistricting plans to grasp an understanding of minority representation across plans. We rely on the Marked Edge Walk to carry out our analysis since the algorithm has the useful property to sample from a targeted distribution.

For a thorough walkthrough of our analysis, please see the following report: [The Impact of the VRA on Minority Representation: What “Race-Blind” Redistricting Looks Like in Southern States](https://docs.google.com/document/d/1MZHpPIdR5elz5MzpLTaWqnESlYOkuR72L6Cy55GbeS4/edit?usp=sharing)

### Summary of the Marked Edge Walk

The Marked Edge Walk is a Markov Chain Monte Carlo (MCMC) algorithm that is implemented for generating ensembles of redistricting plans by sampling from a tuneable targeted distribution. This algorithm utilizes a spanning tree with marked edges, which moves between districting plans through small adjustments to the tree. In comparison to other redistricting algorithms, such as ReCom, the Marked Edge Walk is capable of calculating transition probabilities, making it more desirable in sampling plans from such targeted distributions.

For a thorough walkthrough on the Marked Edge Walk, please see the following paper by Atticus and Professor DeFord: https://arxiv.org/abs/2510.17714v2

## How to MEW (Marked Edge Walk)

From here on, we will explain how to run the Marked Edge Walk (MEW) to replicate results from our analysis. Or, if you're curious, you may use this guide to play around with the parameters or try the analysis on a new state of your choice (e.g. Texas)

The following is the overall workflow of our MEW process:

1. Prepare & Load State Data
	* Load dual graph with demographic & election data, alongside the .shp file and generated .json plans through ReCom.
2. Choose Constraints for the Energy Function
	* Decide on what parameters to be taken into account with the energy function. For our analysis, we set on compactness/cut-edges and county splits.
3. Select an Energy Function
	* Across our eight states, we used either the Gaussian approach or the exponential/minimized approach
4. Set Beta & Target Values
	* Based on the currently enacted plan statistics on compactness and county splits, tune the parameters so we may sample from a target distribution that's similar so we may draw baseline comparisons
5. Base Test & Calibrate
	* Tuning the parameters will come with many trials and error, so it is always good to rerun a test run to feel confident that we are sampling from a distribution that is mixing & as close as to the enacted plan.  
6. Add Partisan constraints
	* note....
8. Partisan Test & Calibrate
	* Once we add partisan, tune the partisanship parameter and run test trials to feel confident that we set appropriate parameters.
9. Run a Lengthy Ensemble
	* Once we are set on parameters, prepare to do a long ensemble run that is from the millions and onwards.

### Repository Files

Before walking through each step in depth, it is essential to first get familiar with the main files in this repository and what they do.

**Table of Content for Main Files**

| File | About |
|---|---|
| `post_do_a_quick_little_thing.ipynb` | Transforms a State's Processed Precincts `.json` file into a dual graph so we may carry out our analysis in Julia. |
| `link_cut_MEW/lct_mew.jl`| Functions for the Marked Edge Walk using Link-Cut Trees. These functions define the energy terms. `make_combined_energy_gaussian` combines county splits and compactness into a single energy function of the Gaussian form: `-beta*((new - target)^2 - (old - target)^2)`. `rep_voteshare_score_vector_gaussian` ranks districts by Republican vote share, from least to most Republican, and assigns each ranked district its own target voteshare. Each district earns a reward as its voteshare climbs toward its target, but gets penalized if it goes past that target. `make_combined_gaussian_party_energy` plugs this republican voteshare gaussian function and is introduced as a term into the energy function from before.
| `link_cut_MEW/lct_run_{STATE}`| These state-by-state files runs the MEW chain from `lct_mew.jl` for each state. Each state has its own version of this file, differing mainly in district count, population tolerance, data path, run length, and the energy function's targets and betas tuned to the state's currently enacted political map. `prepare_warm_start` helps prepare the chain to be ready to begin from either an enacted district plan or from a random seed map generated with ReCom.|
| `{STATE}/`| These state-by-state folders hold each state's data, such as the `.shp` shapefile, its corresponding dual graph, and a set of random seed `.json` plans generated with ReCom MCMC.|
| `batched_runner_warm{STATE}`| These state-by-state batched-runner files run the MEW algorithm for each state. |
| `push_reader_claude_w_cuts_{STATE}`| These state-by-state files turns the raw MEW ensemble file into accessible data. For every run, at each step the partition's data is recorded, such as sorted Democratic vote shares, sorted Black population shares, county split count, and cut-edges count. After the MEW ensemble data is saved, plotting exists that generates plots comparing Democratic vote share per district, as well as Black population vote share. Also code is present to plot compactness and cut-edges counts throughout the ensemble run.|
| `seed_matcher{STATE}`| This state-by-state file is a safety net to make sure that a seed is useable as a warm start for the Marked Edge Walk.| 

**Note on {STATE}:** {STATE} is a placeholder for a state's abbreviation (e.g. `AL`, `FL`, `GA`). Each state has its own copy of these files following the same naming convention.


### In-Depth Step-by-Step Workflow of MEW

1. Prepare & Load State Data

&emsp;Before running the Marked Edge Walk, the state's processed precinct data needs to be converted into a graph format Julia can work with. This is done with a small Python preprocessing script using the `networkx` package. We load in a state's processed precincts file (e.g., AL/AL_Processed_Precincts.json) as a `.json` file and rebuild it as a `networkx` graph object. Afterwards, we convert the string-labeled `GEOID` nodes into integers. Finally, the relabeled graph is exported in a node-link format. The resulting dual graph is what gets loaded into Julia (e.g. in `lct_run_AL` as `AL/AL_processed_precincts_Julia.json`) to use for our MEW ensemble proceedure.

&emsp; *A couple of pitfalls for this process:*
* If you are working with an updated version of `networkx`, the following piece of code in the script, `json_graph.node_link_data(g)`, may not label the `links` data type in the `.json` as "`links`." Rather, `edges` will be the key name of the "`links`" list in the dictionary. We handled this by renaming the `edges` key name to be `links` since our written Julia code used the `links` key name when we loaded in all of the states' data.
* For Louisiana, technical issues arose in attempting to convert the processed precinct data into the dual graph we want for the Marked Edge Walk in Julia. Since we could not assemble a dual graph for the state, this means we won't have the `links` data structure in our dictionary. To work about this, we loaded Louisiana’s processed precincts `.json` in and built the graph from its `adjacency` list in the file. Since each node’s `adjacency` list stores the information of the neighboring node IDS, we iterate through every node, and for each neighbor listed in `adjacency`, we add an edge to the graph and store that edge’s shared perimeter, `share_perim`,  in a `perim_dict`. This allows us to obtain the same dual graph that the python script `json_graph.node_link_data(g)` would have given us. 

2. Choose Constraints for the Energy Function
	* Based on the question you are asking, you may want to sit with thinking about which exact parameters you want to embed in your energy function. (we should describe the purpose of what these constraints do in the energy function.)

3. Select an Energy Function
	* (first write a bit of this section and ask for help on what may be improved) 

**ROUGH DRAFT**

&emsp;Once a set of constraints were chosen for approaching our analysis, step two, the next step is deciding how those constraints get integrated into the Marked Edge Walk as a guide to sample from plans we care about. This leads us into discussing the importance of the energy function’s role in the Marked Edge Walk. In essence, the energy function samples from a distribution to grab a districting plan and returns a score that represents how good the plan is relative to the constraints we care about. At each step, the Marked Edge Walk proposes a tiny change to the current plan (enacted / random seed from the init), and the energy function compares the score of the new plan against the old one to decide on how likely that change is to be accepted. The energy function is essentially the powerhouse of the Marked Edge Walk. 

Across our eight states, we used two approaches for how we structure the energy function, which are the Exponential/Minimized approach and a Gaussian approach.   

&emsp;Table of States' Energy Function Approach 
| Energy Function Approach | States |
|---|---|
| Exponential / Minimized| Florida, Georgia, North Carolina, South Carolina|
| Gaussian | Alabama, Louisiana, Mississippi, Tennessee

* Exponential/Minimized: 

* Gaussian: In the Gaussian approach, each constraint is scored based on its square distance from a target value. 

&emsp;`-β county_splits * ((cnty_splits_new - cty_target)^2 - (cnty_splits_old - cty_target)^2) - β_cuts * ((cuts_new - target_cuts)^2 - (cuts_old - target_cuts)^2)`

&emsp;Inside each energy function, we are able to control how strongly each constraint is (weighted?) in the energy function. This allows us to communicate with the Marked Edge Walk on how much to care about a given constraint relative to the others. 
	
&emsp;Scaling Betas
* Lowering a beta loosens the constraints' influence, allowing the chain more freedom to drift from the target value in favor of exploring other parts of the state space or satisfying other constraints. On the other hand, toggling the beta towards a higher value makes the constraint play a more important part to the energy function. This means that during the walk, the plans that keep that constraint close to its target are favored. Understanding the power(?) scaling(?) Betas are crucial when it comes to tuning our parameters (see Step 4 & Step 5).


4. Set Beta & Target Values

&emsp;After you settle on the constraints for the energy function, as well as choosing the kind of energy function (Exponential & Gaussian) for ..., we will set our parameters for our beta and target constraints. For our analysis, since we are using cut-edges as a compactness measure in our energy function, as well as the total number of county splits across a state, that will be the focus of our *How to Mew* repository. For this section on how to use the Marked Edge Walk, there will be many trials and errors, and most likely you will not find the optimal parameters in the first shot. This step goes in hand with step five, so this section is primarily a forewarning of how lengthy this process may be and guide you on how to figure out the initial parameters.

&emsp;A good start for finding an initial set of parameters for the energy function will be to use each state's geographic measures statistics on the currently enacted districting plan as a starting point; so for our analysis, we will be observing where each state's county splits and cut-edges sit around. These statistics are important for tuning the target parameters, so to go about with county splits and cut-edges, set the target parameters for county splits at the enacted statistic and for cut-edges approximately about where the enacted plan lies at. For convention, we set the beta parameters to be at 0.2 for county splits and 0.010 for cut edges. The number of steps to have the Marked Edge Walk to do, to get a rough idea of how optimal our parameters are for our analysis, a 100,000 through 500,000 step walk will be sufficient.

&emsp;Since in this section we introduced on how to get started with finding these parameters, the following section will continue with this section, elaborating what's the next step after doing an initial run on test parameters and how to identify if the parameters are valid, or how to approach towards optimal parameters. 
	
&emsp;*Acquiring Geographic Statistics:*

* For each state, we acquired our statistics through Professor DeFord's [CISER GerryChain](https://github.com/drdeford/CISER_GerryChain/blob/main/2_Compactness.ipynb) repository. An example script is present in our repository that streamlines this process of acquiring the necessary statistics for our analysis, but for a thorough overview, please refer to the provided repository by Professor DeFord.


5. Base Test & Calibrate

&emsp;(Very rough rough draft) Now that we ran our first short walk ensemble with these intial parameters for our betas and target values, we will discuss about the process to finding parameters that will be optimal for each state. This will be a burdensome stress, since not only do we need to grasp an understanding of how our beta values function in relation to one another, but we also have to come to understand some other underlaying factors, such as the epsilon (population balance) for the state and where the Marked Edge Walk starts exploring from. (note, not really burdensome, but just a bit fustrating understanding what went exactly is not going smoothly) (needs to be better written lol but just write : note to self). 

Now that we ran our first short walk ensemble with these intial parameters for our betas and target values, we will discuss about the process to finding parameters that will be optimal for each state. This portion of the journey will be a bit lackluster due to the fact many Marked Edge Walk ensembles will have to be executed before finally finding ideal parameters. Do not worry though as we will describe in this paper in our best abilities how to proceeed (redundant get to the point; just for my writing sake). Visualization plots will come into play. We must construct trace plots of our county splits and cut-edges values throughout the ensemble. Ideally, for our trace plots for both of our parameters, they should be tracing(?) converging(?) towards the enacted map geographic statistics. If on the first shot it appears that the statistics are converging, congratulations! That's fantastic. Now we must...

Continued on google, will report back to here

Issues with this:
population balance, where the MEW starts from, whether it looks like the trace plots for Democratic % are really moving. 
# mention the 





### notes

* **for step 4 & 5 (probably step 5) mention the use of starting the MEW from either the currently enacted districting plan or a randomly generated seed from ReCom algorithm**


* **low amount of county splits compared to recom in MEW 
	but for party aware experiment, we turned it off so we would have better mixing.**

* **define what are county splits / cut edges (and why we chose these constraints)**

* **go back and mention the exact files / scrips to run and how exactly**





### Repository Structure

Maybe not important, throwout if not necessary [keep right now though]
```text
.
├── README.md
├── Project.toml
├── State
│   ├── dual graph
│   ├── shp
│   ├── seed json
├── post_do_a_quick_little_thing.ipynb
├── link_cut_MEW
│   ├── lct_run_{STATE}.jl
│   ├── lct_MEW.jl
└── Marked_edges
    └── beano2.2_WI.jl
```

Equations

Gaussian

-β<sub>county_splits</sub> • ((county_splits<sub>new</sub> - target_county)<sup>2</sup> - (county_splits<sub>old</sub> - target_county)<sup>2</sup>) - β<sub>cuts</sub> • ((cuts_edges<sub>new</sub> - target_cuts)<sup>2</sup> - (cuts_edges<sub>old</sub> - target_cuts)<sup>2</sup>)

Exponential: 

Have a clear runnable version so if someone stumbles upon the repo, they just want to :heckit



&emsp;Table of States' Energy Function Approach 
| Energy Function Approach | States |
|---|---|
| Exponential / Minimized| Florida, Georgia, North Carolina, South Carolina|
| Gaussian | Alabama, Louisiana, Mississippi, Tennessee




