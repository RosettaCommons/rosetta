## AUTHOR AND DATE
#### Who set up the benchmark? Please add name, email, PI, month and year
Johanna Tiemann: johanna.tiemann@gmail.com, Stein lab, July 2021
Hope Woods: hope.woods@vanderbilt.edu, Meiler Lab
Julia Koehler Leman: julia.koehler.leman@gmail.com, Bonneau lab
 
## PURPOSE OF THE TEST
#### What does the benchmark test and why?
This benchmark compares ddG predictions for single point mutations in alpha helical membrane proteins to experimental data. It also compares the results from four different score functions: mpframework_smooth_fa_2012, ref2015_memb, franklin2019, and, as a baseline, ref2015. Here, we follow the protocol by Park, et al. DiMaio (2016) with adaptions by Frenz, et al. Song (2020), which calculates the ddgs in cartesian space.

## BENCHMARK DATASET
#### How many proteins are in the set?
There are 3 proteins in this dataset, covering 150 unique mutations in total. 

#### What dataset are you using? Is it published? If yes, please add a citation.
The data is pooled together in (Kroncke, B. M., Duran, A. M., Mendenhall, J. L., Meiler, J., Blume, J. D. & Sanders, C. R. Documentation of an Imperative to Improve Methods for Predicting Membrane Protein Stability. Biochemistry 55, 5002-5009 (2016)). From the whole dataset, 
we removed Glycophorin A because they look at interface mutations. We also removed beta-barrels from the dataset because the experiments are slightly different, contextually. Meaning for alpha-helical proteins, the experimental ddG measurements unfold the tertiary, but not all secondary structure, while for beta-barrels the experimental measurements completely unfold the barrel. Our dataset consists of:

bacteriorhodopsin: PDBID 1py6, 88 mutations, ddGs from
	reference1 (3 variants) - Cao, Z., Bowie, J. U., Proceedings of the National Academy of Sciences. 2012. doi:10.1073/pnas.1201298109 (no change)
	reference2 (24 variants) - Faham S., et al. Journal of molecular biology. 2004. doi:10.1016/j.jmb.2003.10.041 (multiplied by -1)
	reference3 (16 variants) - Joh, N. H., et al. Nature. 2008. doi:10.1038/nature06977 (multiplied by -1)
	reference4 (6 variants) - Joh, N. H., et al. Journal of the American Chemical Society. 2009. doi:10.1021/ja904711k (multiplied by -1)
	reference5 (16 variants) - Schlebach, J. P., Journal of the American Chemical Society. 2014. doi:10.1021/ja508359n (ddG_unf kin; multiplied by -1)
	reference6 (2 variants) - Yohannan, S., et al. PNAS. 2004. doi:10.1073/pnas.0306077101 (multiplied by -1)
	reference7 (13 variants) - Cao, Z. et al. Bowie, J. U. BBA Biomembranes. 2012 doi:10.1016/j.bbamem.2011.08.019 (ddG bRf-tobRu at Xsds=0.673, multiplied by -1)
	reference8 (8 variants) - Yohannan, S. et al. Bowie, J. U. JMB. 2004 doi:10.1016/j.jmb.2004.06.025 (multiplied by -1)

dsbB: PDBID 2k74, 18 mutations, ddGs from
	reference1 (12 variants) - Otzen, D.E., Protein Engineering, Design, and Selection. 2011. doi:10.1093/protein/gzq079 (ddG_D-N, oxidizing)
	reference2 (6 variants) - Otzen, D.E., Protein Engineering, Design, and Selection. 2011. doi:10.1093/protein/gzq079 (ddG_D-N, reducing)

glpG: PDBID 2xov, 89 mutations, ddGs from
	reference1 (17 variants) - Baker, R., Urban, S. Nat Chem Biol. 2012. doi:10.1038/nchembio.1021 (multiplied by -1)
	reference2 (69 variants) - Paslawski. W., et al. Proceedings of the National Academy of Sciences. 2015. doi:10.1073/pnas.1424751112 (no changes)
	reference3 (3 variants) Gaffney, K. A., Journal of General Physiology. 2019. doi:10.1085/jgp.201812047 (multiplied by -1)

The datasets are analyzed combined and separately. Only the combined datasets trigger a failure of the test.

#### What are the input files? How were they created?
First, structures are downloaded from the OPM website (https://opm.phar.umich.edu/). The PDB is cleaned with Rosetta/tools/protein_tools/scripts/clean_pdb.py and a spanfile is generated with the Rosetta/main/source/bin/mp_span_from_pdb executable. The cleaned pdb along with the spanfile are used as inputs to FastRelax with constraints to starting coordinates, outputting 100 models (in cartesian space). The lowest scoring model is transformed into membrane coordinates, just in case it had shifted during relax.

Mutfiles were created for each mutation following the format online https://new.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer#input-files.

## PROTOCOL
#### State and briefly describe the protocol.
The protocol takes an in cartesian space relaxed structure as input and executes the cartesian_ddG function as described here (https://new.rosettacommons.org/docs/latest/cartesian-ddG). Shortly, the algorithm allows backbone repacking in sequencial space of the mutant of 1 AA and sidechaing relax within 9 A of the mutant site and performs minimization on the whole structure. 
The whole protocol is repeated 5 times and the score with the lowest energy is selected.
For the membrane sfxns, the Membrane Framework is initiated. For comparison, we do not initiate the membrane framework for ref2015.

#### Is there a publication that describes the protocol?
The protocol is described and benchmarked in detail in protocol:
Park, H., Bradley, P., Greisen, P., Jr, Liu, Y., Mulligan, V. K., Kim, D. E., Baker, D., & DiMaio, F. Simultaneous Optimization of Biomolecular Energy Functions on Features from Small Molecules and Macromolecules. Journal of chemical theory and computation, 12(12), 6201–6212. (2016)
It was further optimized by
Frenz, B., Lewis, S. M., King, I., DiMaio, F., Park, H., & Song, Y. Prediction of Protein Mutational Free Energy: Benchmark and Sampling Improvements Increase Classification Accuracy. Frontiers in bioengineering and biotechnology, 8, 558247. (2020)
#### How many CPU hours does this benchmark take approximately?
The protocol runs about (0.5 hrs per mutant) x (150 mutants) = 75 CPU hours.

## PERFORMANCE METRICS
#### What are the performance metrics used and why were they chosen?
The performance metric is the correlation coefficient between measured ddG values and Rosetta predicted ddG values. 

#### How do you define a pass/fail for this test?
The datasets are analyzed combined and separately. Only the combined datasets trigger a failure of the test. A test is passed if the Pearson and Spearman correlation coefficient is larger than the reference one.

#### How were any cutoffs defined?
As pass/fail cutoff we use the Pearson correlation coefficient of the first run minus 0.1.

## KEY RESULTS
#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
We're comparing ot experimentally derived ddG values. 

#### Describe outliers in the dataset.
The protocol currently predicts artificially high values for most mutations to proline.

## DEFINITIONS AND COMMENTS
#### State anything you think is important for someone else to replicate your results.
Relaxation does result in slightly different models, which could lead to deviation of the single ddG values.

## LIMITATIONS
#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc.
The dataset is very small and biases towards mutations to small residues. The correlation with all score functions is quite low (~0.3) compared to the correlation for soluble proteins (~0.6).

#### How could the benchmark be improved?
#### What goals should be hit to make this a "good" benchmark?
A goal for this benchmark (or rather, protocol development) would be to improve the correlation coefficient to closer to that of soluble proteins. While increasing the size and diversity of the dataset would improve this benchmark, measuring ddG values for membrane proteins is quite difficult, so that will depend on available data.
