## AUTHOR AND DATE
#### Who set up the benchmark? Please add name, email, PI, month and year
Configured by Ziv Ben-Aharon @ FurmanLab. June 2019. PI: Ora Furman-Schueler.
Contact info: ziv.benaharon1@mail.huji.ac.il, oraf@ekmd.huji.ac.il

## PURPOSE OF THE TEST
#### What does the benchmark test and why?
This test showcases the execution of the FlexPepDock refine protocol on 2 example targets.

Results compare executions starting from different receptor conformations:
(1) the bound (native) receptor conformation, 
(2) the unbound (free) receptor conformation,
and (3) the unbound receptor conformation, including receptor backbone minimization. 
The same native peptide starting conformation is used for all three simulations.

## BENCHMARK DATASET
#### How many proteins are in the set?
#### What dataset are you using? Is it published? If yes, please add a citation.
#### What are the input files? How were the they created?
The dataset includes:
1CZY_bound 		- the bound (native) conformation. Receptor chain C, peptide E.
1CZY_ub 		- this is receptor 1CA4 chain C with peptide from 1CZY:E
1CZY_ub_bb		- same as above but with backbone minimization protocl.

2FOJ_bound 		- the bound (native) conformation. Receptor chain A, peptide B.
2FOJ_ub 		- this is receptor 2FOJ chain A with peptide from 2FOJ:B
2FOJ_ub_bb		- same as above but with backbone minimization protocl.

## PROTOCOL
#### State and briefly describe the protocol.
#### Is there a publication that describes the protocol?
#### How many CPU hours does this benchmark take approximately?
Rosetta FlexPepDock is a high-resolution peptide-protein docking protocol that is able to refine a coarse starting structure of a peptide-protein complex, to a near-native model of the interaction. The full degrees of freedom of the peptide are optimized (rigid body orientation, peptide backbone and side chains). Optionally, the receptor backbone can be minimized during optimization.

For more information, read: 
1) Raveh, B., London, N. & Schueler-Furman, O.: Sub-angstrom modeling of complexes between flexible peptides and globular proteins. Proteins (2010). https://onlinelibrary.wiley.com/doi/full/10.1002/prot.22716, and 
2) Alam, N. & Schueler-Furman, O. Modeling peptide-protein structure and binding using monte carlo sampling approaches: Rosetta flexpepdock and flexpepbind. in Methods in Molecular Biology (2017).

Running this benchmark takes approximately 70 CPU hours.

## PERFORMANCE METRICS
#### What are the performance metrics used and why were they chosen?
#### How do you define a pass/fail for this test?
#### How were any cutoffs defined?

A passing test means that at least 8 out of the top10-scoring points are below the interface RMSD cutoff. The cutoff was determined as the maximum of the top10-scoring points in revision 61451 from 2020-Oct-12. The interface RMSD is represented by the rmsBB_if tag in the score file. 

## KEY RESULTS
#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
#### Describe outliers in the dataset. 
The default implementation of FlexPepDock can reliably refine peptide conformations to near-native resolution for starting structures that are up to 4-5 Angstrom away from the native conformation. An extended range can be obtained by including an ab initio search of the peptide backbone conformation (see Raveh et al. PLoSOne 2011), or by coupling the FlexPepDock refinement step to a low-resolution fast rigid body search using other approaches (e.g. FFT, as implemented in  PIPER-FlexPepDock, Alam PlosCB 2017).

## DEFINITIONS AND COMMENTS
#### State anything you think is important for someone else to replicate your results. 
We use the reweighted_score to rank different models: this score term is the sum of total score, interface score and peptide score, providing more weight to the energy terms contributed by the peptide, compared to the energy of the full complex. Alternatively, interface score (I_Sc) may be used. Both terms will reduce the influence of possible conformational changes far away from the binding site that introduce noise.

## LIMITATIONS
#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc. 
#### How could the benchmark be improved?
#### What goals should be hit to make this a "good" benchmark?
The protocol needs an approximate starting conformation of the peptide. This starting structure may be obtained from structures of homologous complexes, or from low-resolution docking protocols.

