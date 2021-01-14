## AUTHOR AND DATE
#### Who set up the benchmark? Please add name, email, PI, month and year

Adapted for the current benchmarking framework by Shourya S. Roy Burman (ssrburman@gmail.com; Gray Lab), Sep 2019 

## PURPOSE OF THE TEST
#### What does the benchmark test and why?

This benchmark is meant to test how well we discriminate native protein-protein binding orientations from decoys based on the interface score term by performing standard protein-protein docking experiments across a diverse set of protein-protein complexes. 

## BENCHMARK DATASET
#### How many proteins are in the set?
#### What dataset are you using? Is it published? If yes, please add a citation.
#### What are the input files? How were the they created?

The dataset consists of 10 protein-protein complexes extracted from the Docking Benchmark 5.0 (Vreven, T. et al. J. Mol. Biol., 2015). The set contains 4 rigid (conformational change < 1.5 Ang), 4 medium-flexible (conformational change between 1.5 and 2.2 Ang), and 2 difficult (conformational change > 2.2 Ang).

Structure preparation:
Proteins were extracted directly from the PDB according to the PDB ID of the native structure. They were then cleaned to ensure that both the bound and the unbound structures have the same residues. The unbound structures were superimposed on the bound and then the smaller partner (ligand) was moved away by 15 Ang and rotated by 60 degrees to scramble the interface

## PROTOCOL
#### State and briefly describe the protocol.
#### Is there a publication that describes the protocol?
#### How many CPU hours does this benchmark take approximately?

Protocol (as per the publication below): "RosettaDock is a Monte Carlo-plus-minimization algorithm consisting of a low-resolution stage, which simulates conformer selection during the formation of the encounter complex, followed by a high-resolution stage, which simulates induced fit in the bound complex. To produce a variety of starting states for the different trajectories, the ligand (the smaller protein) is first randomly rotated and translated about the receptor (the larger protein). In the low-resolution stage, side chains are replaced by coarse-grained “pseudoatoms”, allowing the ligand to efficiently sample the interface by rigid-body movements in a smoothened energy landscape. These rigid-body moves are coupled with backbone conformation swaps where the current backbone conformations of the ligand and the receptor are swapped with different ones from a pre-generated ensemble of conformations. In the high-resolution stage, the side chains are reintroduced to the putative encounter complex and those at the interface are packed for tight binding. There is minimal rigid-body motion in this second stage."

Publication:
The methodological details of the protocol, RosettaDock 4.0 and the performance on the benchmark have been thoroughly discussed in Marze, N. A., Roy Burman, S. S. et al. <a href="https://doi.org/10.1093/bioinformatics/bty355">Bioinfo., 2018</a>. 

CPU hours:
The benchmark takes ~150 hours. The debug mode takes ~4 mins.

## PERFORMANCE METRICS
#### What are the performance metrics used and why were they chosen?
#### How do you define a pass/fail for this test?
#### How were any cutoffs defined?

To assess the performance of a docking simulation, we assess the number of structures with a CAPRI-acceptable or better ranking, which the protocol auto-generates from the native structure. After bootstrap resampling, we classify complexes having N5>=3 as successful.

CAPRI model rankings are based on a combination of factors like fraction of native contacts, interface RMSD, and ligand RMSD.

## KEY RESULTS
#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
#### Describe outliers in the dataset. 

Unbound structures are docked and compared to the bound native structure.

## DEFINITIONS AND COMMENTS
#### State anything you think is important for someone else to replicate your results. 

For best results, use ensembles.

## LIMITATIONS
#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc. 
#### How could the benchmark be improved?
#### What goals should be hit to make this a "good" benchmark?

This is a small part of a much larger standard docking benchmark used in the docking community (Vreven, T. et al. J. Mol. Biol., 2015). More targets from that benchmark will help.
