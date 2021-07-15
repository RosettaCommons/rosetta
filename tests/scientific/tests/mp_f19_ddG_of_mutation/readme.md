## AUTHOR AND DATE
Rebecca F. Alford (ralford3@jhu.edu)
Rituparna Samanta (rsamant2@jhu.edu)
PI: Jeffrey J. Gray (Johns Hopkins ChemBE)
Test created 6/6/19

## PURPOSE OF THE TEST
The purpose of this test is to evaluate the scientific performance of franklin2019, the default energy function for membrane protein structure prediction and design. Specifically, this test evaluates the ability of the energy function to reproduce experimentally measured ddG of mutation values.

## BENCHMARK DATASET
The dataset includes three sets of ddG measurements taken in two protein scaffolds: the Outer Membrane Protein Phospholipiase A (OmpLA; 1qd6) [1] and the Palmitoyl transferase PagP (3gp6) [3]. The ddG measurements for OmpLA and PagP capture the change in free energy upon mutation from alanine to one of the 19 canonical amino acids at a lipid facing site on the protein. The ddG measurements for OmpLA_aro capture the change in free energy upon mutation from alanine to the aromatic amino acids Trp, Tyr, and Phe [2]. 

The references for each set of measurements are given here:

	[1] Moon CP, Fleming KG (2011) "Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers" Proc Natl Acad Sci 108(25):10174-7.
	[2] McDonald, S. K. & Fleming, K. G. Aromatic Side Chain Water-to-Lipid Transfer Free Energies Show a Depth Dependence across the Membrane Normal. J. Am. Chem. Soc. 138, 7946–7950 (2016).
	[3] Marx DC, Fleming KG (2017) "Influence of Protein Scaffold on Side-Chain Transfer Free Energies" Biophysical Journal 113(3):597-604

## PROTOCOL
This test uses the ddG prediction protocol described in Alford & Koehler Leman et al. [3]. Here, a mutation is introduced at the host site and the side chains are optimized within 8Å of the mutated residue. Then, the ddG of mutation is computed as the difference in energy between the mutant and native conformations. In this test, the fa_water_to_bilayer weight value of 1.5, instead of 0.5 used for other tests to improve the correlations for OmpLA and PagP. 

The ddG of mutation protocol is described in:

	Alford RF, Koehler Leman J, Weitzner BD, Duran AM, Tiley DC, Elazar A, Gray JJ (2015) "An integrated framework advancing membrane protein modeling and design" PLoS Comput. Biol. 11(9):e1004398.

The test results for OmpLA and PagP for franklin2019 are described in:

	Alford RF, Fleming PJ, Fleming KG, Gray JJ (2020) "Protein structure prediction and design in a biologically-realistic implicit membrane" Biophys J 180(8): 2042-2055
	
The test results for OmpLA_aro for franklin2019 are described in:

	Alford RF, Samanta R, Gray JJ (2021) "Diverse scientific benchmarks for implicit membrane energy functions" Journal of Chemical Theory and Computation (under review)

The test should take approximately 1 CPU hour.

## PERFORMANCE METRICS
The performance metric for this test is the Pearson correlation coefficient between the experimentally measured and predicted ddG of mutation values. Pass/fail is defined by the correlation coefficient and the fitting parameters (slope and intercept) falling within +/- 0.05  and +/-3 of the values established in Alford et al. 2019. 

## KEY RESULTS
The Pearson correlation coefficients for Franklin2019 ddG predictions are R = 0.922 / 0.158 / 0.919 for OmpLA, OmpLA_aro and PagP, respectively. These results are described in detail in the benchmark paper Alford et al. 2019. The correlation coefficient for OmpLA is high due to the fact that the water to bilayer transfer energy of franklin2019 was caliberated on the basis of its experimental values [1]. The encouraging results for PagP is partly due to the chemically identical regions of the bilayer to that of OmpLA around the substituing site as discussed in [3]. Contrarily, the low value of pearson correlation coefficient for OmpLA_aro are due to the high energy values. The large contributions are from rotamer energies, suggesting steric clashes between the guest side chain and the neighbouring side chains. 

## DEFINITIONS AND COMMENTS
The starting structures were downloaded from the Orientations in Membrane proteins database and the DLPC lipid composition parameters were used for calculations.

## LIMITATIONS
It should be noted that these ddG measurements are high quality data. This is because the measurements were taken in a reversibly folding scaffold and the lipid composition between the experiment and simulation are consistent. Possible improvements are to use the MAE (mean absolute error) of the ddG values directly as a quality metric instead of just for the correlation. Further possible improvements are to include electrostatic terms in franklin2019, to capture the effect of pH and the interaction of peptides with the lipid head groups. 
