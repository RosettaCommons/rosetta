## AUTHOR AND DATE
Rebecca F. Alford (ralford3@jhu.edu)
PI: Jeffrey J. Gray (Johns Hopkins ChemBE)
Test created 6/6/19

## PURPOSE OF THE TEST
The purpose of this test is to evaluate the scientific performance of franklin2019, the default energy function for membrane protein structure prediction and design.

## BENCHMARK DATASET
The benchmark dataset includes five transmembrane peptides. Each peptide has an NMR structure and a tilt angle measured by solid-state NMR. The first four targets (2nr1, 1a11, 1mp6, and 1pje) were taken from [1]. The fifth target, WALP3. was taken from [2].

	1. Ulmschneider MB, Sansom MSP, Di Nola A (2005) "Evaluating tilt angles of membrane-associated helices: comparison of computational and NMR techniques" Biophys J 90(5): 1650-1660
	2. Holt A, Koehorst RBM, Rutters-Meijneke T, Gelb MH, Rijkers DTS, Hemminga MA, Antoinette-Killian J (2009) "Tilt and rotation angles of a transmembrane model peptide as studied by fluorescence spectroscopy" Biophys J 97(8): 2258-2266

The input files are the PDB coordinates for each peptide downloaded from the Protein Databank. Instead of spanfiles, we specify the "single_TM_mode" option, indicating that the protein is a single-TM segment. The PDB coordinate files were cleaned using the `clean_pdb.py` script.

## PROTOCOL
This scientific benchmark identifies low-free energy peptide orientation by calculating an energy landscape: a mapping between all possible orientations relative to the membrane and their energies. Orientation is defined by two coordinates: distance between the center of the membrane and center-of-mass of the peptide, and (2) angle between the membrane normal and helical axis. The protocol first applies side-chain packing and minimization to resolve steric clashes in the peptide structure. Then, rigid body moves are applied to sample all combinations of angle and depth values. Membrane depth values are sampled between -60Å and 60Å with a 1Å step size and tilt angles are sampled between 0-360 degree with a 1 degree step size.

For each peptide, the protocol takes approximately 1-2 CPU hours (for a total 10 CPU hours).

The membrane energy landscape sampling protocol is described in:

	(Alford RF, Fleming PJ, Fleming KG, Gray JJ (2019) "Protein structure prediction and design in a biologically-realistic implicit membrane" Under Review)

## PERFORMANCE METRICS
There are two key performance metrics for this test. The first is the water-to-bilayer partitioning energy of the peptide, calculated as the energy difference between a peptide submerged in water and a peptide oriented vertically in the membrane. The test passes if the calculated partitioning energy is within +/- 1REU of the previously calculated value. The second metric is the minimum energy tilt angle and membrane depth. The test passes if the calculated orientation is within +/- 1Å and +/- 1˚ of the previously calculated orientation. In addition, the energy landscape mapping will be plotted and checked through visual inspection.

## KEY RESULTS
The water-to-bilayer partitioning energy, delta G should be less than 0 for all peptides. Further, the low-energy peptide orientations should be between 0-45 degrees tilted relative to the membrane normal and no more than +/- 5Å from the membrane center. The results are compared to two prior membrane energy functions, a knowledge-based model and an organic-sovlent based model in Alford et al. 2019. 

## DEFINITIONS AND COMMENTS
The membrane normal and center, in addition to the coordinate frame setup are described in:

	Alford RF, Koehler Leman J, Weitzner BD, Duran Am, Tiley DC, Elazar A, Gray JJ (2015) "An integrated framework advancing membrane protein modeling and design" PLoS Comput. Biol. 11(9):e1004398.

## LIMITATIONS
I am working on better quantification metrics for the full landscape.
