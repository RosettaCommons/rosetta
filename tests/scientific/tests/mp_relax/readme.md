## AUTHOR AND DATE
#### Who set up the benchmark? Please add name, email, PI, month and year
The benchmark was set up by Julia Koehler Leman (julia.koehler.leman@gmail.com) in March 2019.
The research was performed in Jeff Gray's lab, but Julia Koehler Leman is now in the lab of Richard Bonneau. 

## PURPOSE OF THE TEST
#### What does the benchmark test and why?
We test whether how much flexibility we can sample with mp_relax. The mp_relax app runs FastRelax under the hood but by using the membrane framework RosettaMP. 

## BENCHMARK DATASET
#### How many proteins are in the set?
#### What dataset are you using? Is it published? If yes, please add a citation.
#### What are the input files? How were the they created?
The benchmark set contains 4 test proteins that were run for mp_relax from the RosettaMP framework paper (Alford, Koehler Leman et. al, PlosCompBio, 2015). Note that the mp_relax application could be further optimized or compared to the newer protocol mp_range_relax that hasn't been published yet. The paper describes mp_relax and other applications as a proof-of-concept.

Input PDBs were downloaded from the PDBTM database (which have better membrane embedding than PDBs from the OPM database). This means that the proteins are transformed into membrane coordinates and have the correct biomolecular assembly, even though we are only looking at dimers for this application. PDBs were cleaned using clean_pdb.py from the  Rosetta/tools/protein_tools/scripts directory. Spanfiles were created using the mp_span_from_pdb application in Rosetta. 

The PDB files were used as natives for RMSD calculations. Detail command lines are described in the Supplement to (Alford, Koehler Leman et. al, PlosCompBio, 2015).

## PROTOCOL
#### State and briefly describe the protocol.
#### Is there a publication that describes the protocol?
#### How many CPU hours does this benchmark take approximately?

The mp_relax application runs in RosettaScript and uses the FastRelax protocol as a mover but adjusted with the membrane framework RosettaMP. We create 100 models per protein to sample protein conformations. We use the mp_framework_smooth_2012 scorefunction for scoring until something better becomes available. 

[Note that the natives should have a membrane residue (MEM), otherwise Rosetta crashes during RMSD calculation.]

On average for this benchmark set, mp_relax creates a decoy in ~350s; that is dependent on size of the protein. That makes 350s x 4 proteins x 100 decoys = 140,000s which are <40 CPU hours. 

## PERFORMANCE METRICS
#### What are the performance metrics used and why were they chosen?
#### How do you define a pass/fail for this test?
#### How were any cutoffs defined?

We use the total Rosetta score and RMSD to the native structure as we are interested in the sampling range and the score ranges in sampling. A passed test constitutes 90% of the decoys generated have a smaller than the defined score and RMSD cutoffs. The cutoffs were defined by looking at these ranges in a single run without outlier decoys. Specifically, for RMSD we use the maximum RMSD + stdev. For score, we use the maximum score + stdev. 

## KEY RESULTS
#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
#### Describe outliers in the dataset. 

We use crystal and structures as natives for comparison. The sampling range is TO EDIT

2kse is an NMR structure of a histidine kinase, which normally has 2 chains but is here modeled as a single chain. That's bad. Not sure why this was chosen. 
2leg is the NMR structure of DsbB, not a great structure and has a long loop which flops around quite a bit
3pxo is the crystal structure of metarhodopsin, a decent structure
4a2n is a crystal sructure with relatively long loops which don't seem to be super floppy - Rosetta predicts a pretty narrow range

## DEFINITIONS AND COMMENTS
#### State anything you think is important for someone else to replicate your results. 

## LIMITATIONS
#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc. 
#### How could the benchmark be improved?
#### What goals should be hit to make this a "good" benchmark?

The mp_relax app needs to be improved (or moved to mp_range_relax) and to benchmark that properly, a larger dataset is needed. From that, a collection of good, intermediate and outlier proteins should be picked for continuous scientific benchmarking.

The protein embedding in the membrane should be optimized during relax. Since the FoldTree in mp_relax is such that the MEM is relaxed around the protein, doing that and superimposing the structures leads to MEM residues all over the place. The protein set is not very well picked. We need a broader set including beta-barrels and larger ranges of protein size with high-quality structures. There are NMR structures in here, which we might or might not want to avoid. We need to test multi-chain proteins in relax. And we need to test against range_relax and mp_range_relax which runs much master and provides better results. 
