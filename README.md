BacArena
========

BacArena is an agent based modeling framework for microbial communities.

See <a href="./doc/documentation.pdf">Documentation</a>


Microbial communities are essential for global ecosystems and human health.
Computational modeling of microbial consortia is thus a major goal in systems biology and microbial ecology. 

<i>BacArena</i> is a project to simulate bacterial behaviour in communities. A lot of progress is done in the last years to gain genome wide metabolic reconstructions of certain organisms, which open a wide field of mathematical analysis [1], [2]. <br>
One of this new methods is flux balanced analysis (fba) [3] to estimate optimal metabolic fluxes under certain constraints. By this work advanced models are possible, which are available in a defined, exchangeable format (<i>SBML</i>) [5]. The idea of this project is to use this existing reconstructions and put them in a spatial and temporal environment to study their possible interactions.
This is achieved by the combination of agent based modeling [4] with fba. Each bacterium is considered as an agent with individual states, own properties and rules to act. Agents are located on a grid where they can move and interact via metabolic exchanges computed by fba.

The starting point for our project is curiosity of what could be done with this huge models. We just throw those models into an arena to see what kind of actions will evolve.


[1] Lewis et. al ,,Constraining the metabolic genotype–phenotype relationship using a phylogeny of in silico methods'' Nature Reviews Microbiology 10, 291-305 (April 2012)<br>
[2] Feist et. al ,,Available predictive genome-scale metabolic network reconstructions'' http://systemsbiology.ucsd.edu/InSilicoOrganisms/OtherOrganisms <br>
[3] http://en.wikipedia.org/wiki/Flux_balance_analysis <br>
[4] http://en.wikipedia.org/wiki/Agent-based_model <br>
[5] http://en.wikipedia.org/wiki/SBML <br>
