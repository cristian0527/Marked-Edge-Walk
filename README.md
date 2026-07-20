# The Marked Edge Walk (MEW) for a Statistical Analysis of Minority Representation after Callais

This repository is the keepsake of Professor Daryl DeFord's Undergraduate Research Summer Institute group at Vassar College, which endeavors to analyze the future of Southern States post Louisiana v. Callais.

First of all, we would like to thank Atticus McWhorter for joining us in our summer research journey. Atticus, alongside Professor DeFord, guided us in implementing the Marked Edge Walk for our analysis and taught us how to program in Julia.

### Overview

[quick overview of what our purpose is with our analysis, why we are studying this, and discuss why we are using MEW.] 

In late April 2026, the Supreme Court decided a Louisiana congressional map was unconstitutional due to race being a predominant factor in creating the districts. This case, Louisiana v. Callais, significantly reduces protections for minority representation through the Voting Rights Act. Following this case, there have been questions, concerns, and speculation surrounding how this case will impact minority representation throughout the case. One recent New York Times article written by Nate Cohn claimed that the Voting Rights Act did not create additional representation for minority voters and that race blind redistricting would maintain this representation.

This repository analyzes what redistricting would look like in those eight states Cohn argued in his NYT article through the Marked Edge Walk. These eight states are Alabama, Georgia, Florida, Louisiana, Mississippi, North Carolina, South Carolina, and Tennessee. In our analysis, we find sampling distributions that are similar to the current enacted congressional plan in each state, and then generate large ensembles of redistricting plans to grasp an understanding of minority representation across plans. We rely on the Marked Edge Walk to carry out our analysis since the algorithm has the useful property to sample from a targeted distribution.

For a thorough walk of our analysis, please see the following report of our analysis: __link__

### Summary of the Marked Edge Walk

The Marked Edge Walk is a Markov Chain Monte Carlo (MCMC) algorithm that is implemented for generating ensembles of redistricting plans by sampling from a tuneable targeted distribution. This algorithm utilizes a spanning tree with marked edges, which moves between districting plans through small adjustments to the tree. In comparison to other redistricting algorithms, such as ReCom, the Marked Edge Walk is capable of calculating transition probabilities, making it more desirable in sampling plans from such targeted distributions.

For a thorough walkthrough on the Marked Edge Walk, please see the following paper by Professor DeFord and Atticus: https://arxiv.org/abs/2510.17714v2


## How to MEW (Marked Edge Walk)


### Repository Structure

the following is a quick example (consult with what to add and such...) 

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





