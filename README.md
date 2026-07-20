# The Marked Edge Walk (MEW) for a Statistical Analysis of Minority Representation after Callais

This repository is the keepsake of Professor Daryl DeFord's Undergraduate Research Summer Institute group at Vassar College, which endeavors to analyze the future of Southern States post Louisiana v. Callais.

First of all, we would like to thank Atticus McWhorter for joining us in our summer research journey. Atticus, alongside Professor DeFord, guided us in implementing the Marked Edge Walk for our analysis and teaching us how to program in Julia.

The Marked Edge Walk is a Markov Chain Monte Carlo (MCMC) algorithm that is implemented for generating ensembles of redistricting plans by sampling from a tuneable targeted distribution. This algorithm utilizes a spanning tree with marked edges, which moves between districting plans through small adjustments to the tree. In comparison to other redistricting algorithms, such as ReCom, the Marked Edge Walk is capable of calculating transition probabilities, making it more desirable in sampling plans from such targeted distributions.

For a thorough walkthrough on the Mark Edge Walk, please see the following paper by Professor DeFord and Atticus: https://arxiv.org/abs/2510.17714v2



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

