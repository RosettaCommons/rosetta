## AUTHOR AND DATE
Rebecca F. Alford (ralford3@jhu.edu)
PI: Jeffrey J. Gray (Johns Hopkins ChemBE)
Test created 6/6/19

## PURPOSE OF THE TEST
The purpose of this test is to evaluate the scientific performance of franklin2019, the default energy function for membrane protein structure prediction and design. Specifcally, this test evaluates the ability of the energy function to reproduce experimentally measured ddG of mutation values. 

## BENCHMARK DATASET
The dataset includes two sets of 19 ddG measurements taken in two corresponding scaffolds: the Outer Membrane Protein Phospholipiase A (OmpLA; 1qd6) [1] and the Palmitoyl transferase PagP (3gp6) [2]. Each set of ddG measurements capture the change in free energy upon mutation from alanine to one of the 19 canonical amino acids at a lipid facing site on the protien. The references for each set of measurements are given here: 

	[1] Moon CP, Fleming KG (2011) "Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers" Proc Natl Acad Sci 108(25):10174-7.
	[2] Marx DC, Fleming KG (2017) "Influence of Protein Scaffold on Side-Chain Transfer Free Energies" Biophysical Journal 113(3):597-604

## PROTOCOL
This test uses the ddG prediction protocol described in Alford & Koehler Leman et al. [3]. Here, a mutation is introduced at the host site and the side chains are optimized within 8Å of the mutated residue. Then, the ddG of mutation is computed as the difference in energy between the mutant and native conformations. 

The ddG of mutation protocol is described in: 

	(Alford RF, Koehler Leman J, Weitzner BD, Duran Am, Tiley DC, Elazar A, Gray JJ (2015) "An integrated framework advancing membrane protein modeling and design" PLoS Comput. Biol. 11(9):e1004398.)

The established test results for franklin2019 are described in: 

	(Alford RF, Fleming PJ, Fleming KG, Gray JJ (2019) "Protein structure prediction and design in a biologically-realistic implicit membrane" Submitted.)

The test should take approximately 1 CPU hour. 

## PERFORMANCE METRICS
The performance metric for this test is the Pearson correlation coefficient between the experimentally measured and predicted ddG of mutation values. Pass/fail is defined by the MAE falling within +/- 0.5 REU of the values established in Alford et al. 2019. 

## KEY RESULTS
The pearson correlation coefficient for Franklin2019 ddG predictions is R = 0.85 for the OmpLA adataset and R = 0.83 for the PagP dataset. This result is described in detail in Alford et al. 2019. 

## DEFINITIONS AND COMMENTS
The starting structures were downloaded from the Orientations in Membrane proteins database and the DLPC lipid composition parameters were used for calculations. 

## LIMITATIONS
It should be noted that these ddG measurements are high quality data. This is because the measurements were taken in a reversibly folding scaffold and the lipid composition between the experiment and similation are consistent. 


