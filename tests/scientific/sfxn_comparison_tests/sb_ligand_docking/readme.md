## AUTHOR AND DATE
#### Who set up the benchmark? Please add name, email, PI, month and year

Adapted for the current benchmarking framework by Shannon Smith (shannon.t.smith.1@vanderbilt.edu; Meiler Lab), Feb 2019 

## PURPOSE OF THE TEST
#### What does the benchmark test and why?

This benchmark is meant to test how well we discriminate native small molecule binding orientations from decoys based on the interface score term by performing standard ligand docking experiments across 50 diverse protein-ligand complexes. 

## BENCHMARK DATASET
#### How many proteins are in the set?
#### What dataset are you using? Is it published? If yes, please add a citation.
#### What are the input files? How were the they created?

The dataset consists of 50 protein-ligand complexes extracted from the Diverse Platinum Dataset (Friedrich, N. et al. J. Chem. Inf. Model, 2017). In addition to the filters described thoroughly in Friedrich et al., the maximum resolution allowed was reduced to 2.0A, drug-like ligands were chosen using standard Lipinski Rules (Lipinski, C.A., et al. Advanced Drug Delivery Reviews, 1997), and visually inspected for unrealistic orientations and crystallographic artefact. This dataset was also filtered to eliminate cofactors or multiple ligands witin the binding site. The full list of protein-ligand complexes is given in the 1.submit.py file by PDB ID.

## PROTOCOL
#### State and briefly describe the protocol.
#### Is there a publication that describes the protocol?
#### How many CPU hours does this benchmark take approximately?

Structure preparation:
Proteins were extracted directly from the PDB according to their PDBIDs, underwent minimal cleaning to remove the ligand. Note that no prior relax or other structural manipulations were performed prior to beginning the docking run. 

Ligand preparation:
Initial ligand structures were downloaded from the Protein Data Bank using the ligand ideal SDF. Note that this is not the same as the structure from the input cocrystal structure in order to minimize conformational bias towards a particular pre-generated pose. Ligand files were cleaned using OpenBabel (J. Cheminf. 2011, 3, 33), run through BCL Conformer Generator (Kothiwale, Meiler. J. Cheminform., 2015) and generated Rosetta-readable parameter files according to the following scripts:

BCL Conformer Generator
bcl.exe molecule:ConformerGenerator -rotamer_library cod -top_models 100 -ensemble_filenames NAME.sdf -conformers_single_file NAME_conformers.sdf -conformation_comparer 'Dihedral(method=Max)' 30 -max_iterations 1000 

Rosetta Parameter File Generation:
/programs/x86_64-linux/rosetta/3.8/main/source/scripts/python/public/molfile_to_params.py -n $NAME -p $NAME --mm-as-virt --conformers-in-one-file NAME_conformers.sdf --chain X

This protocol currently uses 50 protein-ligand complexes, each containing the input file containing the cleaned protein + input ligand PDB (target_input.pdb), the native protein-ligand complex for RMSD calculations (target_native.pdb), the ligand params file (target_ligand.params) and the ligand conformer library file pointed to by the params file (target_ligand_conformers.pdb). 

NOTE: the protein and ligand preparation steps were performed previously and are given as the input in the data/ directory. In other words, these steps are not performed each time this benchmark is run.

Each test takes ~8 CPU hours x 50 tests = ~400 CPU hours.

## PERFORMANCE METRICS
#### What are the performance metrics used and why were they chosen?
#### How do you define a pass/fail for this test?
#### How were any cutoffs defined?

**Need to run several tests to see what we think should be the cutoff.
The big question that this benchmark intends to test is how well we discriminate native versus non-native binding poses. A run is determined successful if there is an near-native (<2A) structure within the top 1percent of models based on the interface_delta_X score. 

Sampling failure is defined as having no sub-2A output structures. 
Scoring failure is defined as not having a sub-2A output structure within the top 10% ranked by interface score (this is a pretty large margin and may adjust accordingly).

## KEY RESULTS
#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
#### Describe outliers in the dataset. 

This benchmark is meant to be a longitudinal test to determine how changes in the scorefunction impact small-molecule docking performance. 

Performance varies greatly across different test cases, so I am not sure how to go about grading overall performance. This also makes it difficult to define a binary pass/fail criteria to the entire benchmark.

## DEFINITIONS AND COMMENTS
#### State anything you think is important for someone else to replicate your results. 

There is a bit of noise in the runs, so you may need to compare several runs to get better sense of how different scorefunctions compare.

## LIMITATIONS
#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc. 
#### How could the benchmark be improved?
#### What goals should be hit to make this a "good" benchmark?

The run-to-run variablility is rather high, which might indicate that the benchmark could be improved by running several output structures for each input, rather than just the current one.

This benchmark currently utilizes the ligand scorefunction (an off-shoot of the score12 environment), as this is currently the best protocol for ligand docking in Rosetta. Could be updated for talaris2014, but performance has been shown to deteriorate in the newer scoring environments, notably ref2015 and later. 
