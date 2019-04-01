## AUTHOR AND DATE
#### Who set up the benchmark? Please add name, email, PI, month and year
The benchmark was set up by Julia Koehler Leman (julia.koehler.leman@gmail.com) in March 2019.
The research was performed in Jeff Gray's lab, but Julia Koehler Leman is now in the lab of Richard Bonneau. 

## PURPOSE OF THE TEST
#### What does the benchmark test and why?
We test whether we can recapitulate binding interfaces of proteins inside the membrane. For this, we run local docking with the membrane scorefunction from RosettaMP. 

## BENCHMARK DATASET
#### How many proteins are in the set?
#### What dataset are you using? Is it published? If yes, please add a citation.
#### What are the input files? How were the they created?
The benchmark set contains 10 test proteins that were run for mp_dock from the RosettaMP framework paper (Alford, Koehler Leman et. al, PlosCompBio, 2015). Note that the mp_dock application runs local docking (NOT global docking) and the mp_dock application could be further optimized. The paper describes mp_dock and other applications as a proof-of-concept!!!

Input PDBs were downloaded from the PDBTM database (which have better membrane embedding than PDBs from the OPM database). This means that the proteins are transformed into membrane coordinates and have the correct biomolecular assembly, even though we are only looking at dimers for this application. PDBs were cleaned using clean_pdb.py from the  Rosetta/tools/protein_tools/scripts directory. Spanfiles were created using the mp_span_from_pdb application in Rosetta. 

The PDB files were used as natives for RMSD calculations. Additionally, we created prepacked files using the docking_prepack application in Rosetta, picking the top-scoring model from 10 prepacked decoys. Detail command lines are described in the Supplement to (Alford, Koehler Leman et. al, PlosCompBio, 2015).

## PROTOCOL
#### State and briefly describe the protocol.
#### Is there a publication that describes the protocol?
#### How many CPU hours does this benchmark take approximately?

The mp_dock application uses the docking_protocol but adjusted with the membrane framework RosettaMP. mp_dock runs local docking with a docking perturbation of 3A translation and 8 degrees of rotation. We create 1000 models per protein to sample the interface. 

[Note that the natives shouldn't have a membrane residue (MEM), otherwise Rosetta crashes during RMSD calculation.]

On average for this benchmark set, mp_dock creates a decoy in ~70s; that is dependent on size of the protein. That makes 70s x 10 proteins x 1000 decoys = 700,000s which are <200 CPU hours. 

## PERFORMANCE METRICS
#### What are the performance metrics used and why were they chosen?
#### How do you define a pass/fail for this test?
#### How were any cutoffs defined?

We use the interface score vs. RMSD. Interface score has been found to be more predictive than the total Rosetta score. Alternatively, one can also look at the fraction of native contacts (Fnat). Various docking metrics (also the interface RMSD) are output into the score file. The thresholds for various docking metrics can be found in CAPRI papers. 

We define the modeling a success if the the top10 models (out of a 1000) are below a defined RMSD threshold. The RMSD threshold was defined by taking the mean RMSD of the top10 models from the first run and adding 2A. 

## KEY RESULTS
#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
#### Describe outliers in the dataset. 

We use crystal structures as natives for comparison, another general method for MP docking does not really exist - except maybe a specialized tool. 

Docking with this proof-of-concept application is easier for smaller systems and more difficult for larger protein or smaller interfaces. Soluble chains contribute to the interface, which can make docking easier, i.e. give more of a funnel. 

1afo and 2k9j are single TM helices, they are easier to dock
2l35 is a double helix docking to a single TM helix. There is a disulfide bond connecting the N-termini of the two partners which leads to incredibly low I_sc and large RMSDs 
2nq2 actually has 4 chains, for this applications, chains A+C and chains B+D were combined to A and B, respectively. Original chains C and D are the soluble parts of this dimer which create a larger docking interface and for docking results for the paper gave a really nice funnel. Removing chains C and D and therefore the soluble interface might make it more difficult to find a good docking interface, but this hypothesis was not tested.
2qi9 is a multi-helix homodimer
3din has a partner with a long loop that interacts with the top of the other partner
3odu GPCR with a small interface
3rlf very intertwined interface, which is hard to find, therefore RMSDs are large
3tui is a V-shaped transporter with an open and closed conformation. The Xray structure is in the open conformation but the closed conformation has been found to score better in Rosetta, presumably because it has a larger interface
4djh another GPCR with a small interface


## DEFINITIONS AND COMMENTS
#### State anything you think is important for someone else to replicate your results. 

## LIMITATIONS
#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc. 
#### How could the benchmark be improved?
#### What goals should be hit to make this a "good" benchmark?

The mp_dock app needs to be improved and to benchmark that properly, a larger dataset is needed. From that, a collection of good, intermediate and outlier proteins should be picked for continuous scientific benchmarking. The benchmark could further be improved by using ensemble docking or considering flexibility or using homology models instead of crystal structures. 

TODO: look at Fnat (=fraction of native contacts): for 3tui, the V-shaped transporter, the best RMSD structure is around 10A, which looks like a bad prediction. However, the Fnat for this is around 80%, which means most of the native interface has been found. 

We can only dock dimers, not multi-mers. Multi-body docking is currently not implemented. 
