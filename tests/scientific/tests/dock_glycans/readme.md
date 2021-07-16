## AUTHOR AND DATE

THIS TEST WILL NO LONGER BE RUN ON THE TEST SERVER BECAUSE IT HAS BEEN SUPERSEDED BY THE GLYCAN_DOCK TEST!!!

Labonte &lt;JWLabonte@jhu.edu&gt;, Gray Lab, 2019.06

## PURPOSE OF THE TEST
This benchmark confirms that the dock_glycans protocol is able to determine high-quality glycoligand models by CAPRI metrics in a bound&ndash;unbound local, flexible docking run.

## BENCHMARK DATASET
There are 4 protein&ndash;glycoligand complexes in the benchmark set, which were selected from a much larger set found in [Anita K. Nivedha, <i>et al.</i> <cite>J. Chem. Theory Comput.</cite> <b>2016</b>, <i>12</i> 2892&ndash;901]. The four selected complexes were initially chosen for their variety of size and success in the protocol. 

The input files are relaxed, native <code>.pdb</code> files.

## PROTOCOL
The initial protocol being run here is the dock_glycans protocol published with [Labonte, J.W.; Adolf-Bryfogle, J.; Schief, W.R.; Gray, J.J. “Residue-Centric Modeling and Design of Saccharide and Glycoconjugate Structures.” <cite>J. Comput. Chem.</cite> <b>2017</b>, <i>38</i>, 276–287]. 500 decoys are generated for each of the four input structures. An individual constraint file with a single site constraint makes sure that the sugar doesn't move out of the binding pocket. 

## PERFORMANCE METRICS
Interface energy is plotted <i>vs.</i> ligand RMSD. Pass/fail is currently determined by having at least 30% of the decoys below the mean for both ligand RMSD and interface energy, i.e. the lower left quadrant. 

## KEY RESULTS

## DEFINITIONS AND COMMENTS

## LIMITATIONS
After we publish our results with an enhanced GlycanDock protocol on the entire benchmark, we will expand the test here. Ideally, cutoffs should be determined by decoys considered to be of high-quality by CAPRI metrics.
