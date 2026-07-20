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
3. Choose Constraints for the Energy Function
	* Decide on what parameters to be taken into account with the energy function. For our analysis, we set on compactness/cut-edges and county splits.
5. Select an Energy Function
	* Across our eight states, we either were set on the Gaussian approach or the exponential/minimized approach.
7. Set Beta & Target Values
	* Based on the currently enacted plan statistics on compactness and county splits, tune the parameters so we may sample from a target distribution that's similar so we may draw baseline comparisons
9. Test & Calibrate
	* Tuning the parameters will come with many trials and error, so it is always good to rerun a test run to feel confident that we are sampling from a distribution that is mixing & as close as to the enacted plan.  
11. Add Partisan constraints
13. Test & Calibrate
	* Once we add partisan, tune the partisanship parameter and run test trials to feel confident that we set appropriate parameters.
15. Run a Lengthy Ensemble
	* Once we are set on parameters, prepare to do a long ensemble run that is from the millions and onwards.

### Repository Files

Before walking through each step in depth, it is essential to first get familiar with the main files in this repository and what they do.

### Table of Files

| File | About |
|---|---|
| `post_do_a_quick_little_thing.ipynb` | Transforms a State's Processed Precincts `.json` file into a dual graph so we may carry out our analysis in Julia|
| `link_cut_MEW/lct_mew.jl`| Functions for the Marked Edge Walk using Link-Cut Trees. These functions define the energy terms. `make_combined_energy` combines county splits and compactness into a single energy function of the Gaussian form: `-beta*((new - target)^2 - (old - target)^2)`. This function is the original energy function variant. `rep_voteshare_score_vector_gaussian` ... [finish mew] |
| `link_cut_MEW/lct_run_{STATE}`| mew |
| `Marked_edges/beano2.2_WI.jl`| mew |
| `{STATE}/`| mew |
| `batched_runner_warm{STATE}`| mew |
| `push_reader_claude_w_cuts_{STATE}`| mew |
| `seed_matcher{STATE}`| mew |

**STATE:** placeholder - explain what this is 




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





