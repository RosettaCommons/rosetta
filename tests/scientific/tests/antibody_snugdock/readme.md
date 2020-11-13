## AUTHOR AND DATE
SnugDock scientific test first implemented June, 2019 (revised Oct, 2020) by Jeliazko Jeliazkov (Gray Lab, JHU, jeliazkov@jhu.edu).

## PURPOSE OF THE TEST
This test seeks to evalute the executable snugdock. Briefly, this application simulates antibody--antigen interactions by local docking with CDR loop refinement. A simulation consists of fifty cycles of randomly selected moves, including, re-docking the relative orientation of the VH--VL and Ab--Ag are refined, remodeling the CDR H3/H2 loops in the context of the Ag interaction, side-chain repacking, and minimization of all relevant DOF (Ab--Ag interface, VH--VL interface, CDR loops, and side chains).

Here we test the protocol's ability to predict Ab--Ag interactions starting from ensembles of homology models and the unbound antigen crystal structure (when possible). We evaluate success/failure based on the <a href=https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3076516/>CAPRI quality criteria</a> of the lowest energy (based on interface score) models. 

## BENCHMARK DATASET
The dataset comprises six antibody targets from Sircar and Gray (<a href=https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000644>PLoS Comp. Bio. (2010)</a>). The targets are of varying difficulty (1AHW, 1JPS, 1MLC, and 1ZTX should be relatively easy - either sampling near-native states or having a funnel - whereas 2AEP and 2JEL should be challenging - no funnel, no near-native sampling). Native structures are crystals. Inputs for modeling are RosettaAntibody homology models and relaxed (unbound when possible) antigen.

## PROTOCOL
The Rosetta SnugDock protocol is described in our publication (Weitzner, Jeliazkov, Lyskov, et al., <a href=https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5739521> Nat. Protocol. 2017</a>). The benchmark is rather time consuming, taking ~1 hour per model. We aim to produce 500 models for 6 antibodies, so the runtime is 3000 CPU hours.

## PERFORMANCE METRICS
To consider a simulation successful, we would like to sample native states (low rmsd) and identify those states based on low interface energy. For this test to be passed, all simulations must produce at least one model of acceptable quality (according to the CAPRI criteria). Additionally, for the four easy targets, we expect low-energy models to have medium or acceptable quality. Thus, 1ahw, 1jps, and 1ztx must either produce at least 4 acceptable quality models or at least 3 medium quality models in the 10 lowest energy models (by interface score). 2jel, a slightly more challenging target must produce at least 4 acceptable quality models.

## KEY RESULTS
Antibody--antigen bound structure prediction is a challenging task. The current test assesses only six of a possible fifteen targets, which should display a breadth of behavior. For the easier targets (1AHW, 1JPS, 1MLC, and 1ZTX) we expect to observe either near-native (sub 2-Angstrom) sampling of the interface or a good funnel (negative discrimintation). For the harder targets (2AEP and 2JEL) we will not observe either.

## LIMITATIONS
The full assessment for Rosetta SnugDock typically covers the entire protocol from homology modeling to H3 modeling to antigen docking and consists of 15 targets. Due to time considerations, we obviously do not test the full protocol over 15 antibody-antigen complexes here. Instead we focus on six representative targets.
