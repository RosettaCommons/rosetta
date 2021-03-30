## AUTHOR AND DATE
#### Who set up the benchmark? Please add name, email, PI, month and year
Hope Woods: hope.woods@vanderbilt.edu, Meiler Lab, March 2021
Julia Koehler Leman: julia.koehler.leman@gmail.com, Bonneau lab
​Johanna Tiemann, Amelie Stein lab
 
## PURPOSE OF THE TEST
#### What does the benchmark test and why?
This benchmark compares ddG predictions for single point mutations in alpha helical membrane proteins to experimental data. It also compares the results from four different score functions: mpframework_smooth_fa_2012, ref2015_memb, franklin2019, and, as a baseline, ref2015.
​
## BENCHMARK DATASET
#### How many proteins are in the set?
#### What dataset are you using? Is it published? If yes, please add a citation.
#### What are the input files? How were the they created?
There are 3 proteins in this dataset, covering 150 mutations in total. The data is pooled together in (Kroncke, B. M., Duran, A. M., Mendenhall, J. L., Meiler, J., Blume, J. D. & Sanders, C. R. Documentation of an Imperative to Improve Methods for Predicting Membrane Protein Stability. Biochemistry 55, 5002–5009 (2016)). From the whole dataset, 
we removed Glycophorin A because they look at interface mutations. We also removed beta-barrels from the dataset because the experiments are slightly different, contextually. Meaning for alpha-helical proteins, the experimental ddG measurements unfold the tertiary, but not all secondary structure, while for beta-barrels the experimental measurements completely unfold the barrel. Our dataset consists of:
					
bacteriorhodopsin: PDBID 1py6, 67 mutations, ddGs from
	Cao, Z., Bowie, J. U., Proceedings of the National Academy of Sciences. 2012.
	Faham S., et al. Journal of molecular biology. 2004.
	Joh, N. H., et al. Nature. 2008.
	Joh, N. H., et al. Journal of the American Chemical Society. 2009.
	Schlebach, J. P., Journal of the American Chemical Society. 2014.
	Yohannan, S., et al. PNAS. 2004.
​
dsbB: PDBID 2k74, 12 mutations, ddGs from
	Otzen, D.E., Protein Engineering, Design, and Selection. 2011.
​
glpG: PDBID 2xov, 71 mutations, ddGs from
	Baker, R., Urban, S. Nat Chem Biol. 2012.
	Paslawski. W., et al. Proceedings of the National Academy of Sciences. 2015.
	Gaffney, K. A., Journal of General Physiology. 2019.
 
Input Files:
​
First, structures are downloaded from the OPM website (https://opm.phar.umich.edu/). The PDB is cleaned with Rosetta/tools/protein_tools/scripts/clean_pdb.py and a spanfile is generated with the Rosetta/main/source/bin/mp_span_from_pdb executable. The cleaned pdb along with the spanfile are used as inputs to FastRelax with constraints to starting coordinates, outputting 100 models. The lowest scoring model is transformed into membrane coordinates, just in case it had shifted during relax.
​
Resfile were created for each mutation following the format online https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/resfiles.
​
The protocol itself is implemented in RosettaScripts and XML scripts to create mutant and wildtype models that are used to calculate the ddG are provided.
 
## PROTOCOL
#### State and briefly describe the protocol.
#### Is there a publication that describes the protocol?
#### How many CPU hours does this benchmark take approximately?
​
The protocol takes a relaxed structure as input into two different scripts - one for the wildtype and one for the mutant. One introduces a mutation, does backrub sampling on residues within 8 A of the mutation site, then iterates through repacking within 8 A of the mutation site and minimization on the whole structure four times. The other script does the same sampling without introducing the mutation. The scores of the three lowest-scoring models are averaged for mutant and WT and the average score from WT is subtracted from the average score of the mutant.
​
The protocol is not published, but is based on the FlexddG protocol from Barlow, K. A., et al. The Journal of Physical Chemistry B. 2018.
​
The protocol runs about (12 hrs per mutant) x (150 mutants) = 1800 CPU hours.
​
## PERFORMANCE METRICS
#### What are the performance metrics used and why were they chosen?
#### How do you define a pass/fail for this test?
#### How were any cutoffs defined?
The performance metric is the correlation coefficient between measured ddG values and Rosetta predicted ddG values.
​
As pass/fail cutoff we use the Pearson correlation coefficient of the first run minus 0.1.
​​
## KEY RESULTS
#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
#### Describe outliers in the dataset.
​
The protocol currently predicts artificially high values for most mutations to proline.
​
## DEFINITIONS AND COMMENTS
#### State anything you think is important for someone else to replicate your results.
​
## LIMITATIONS
#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc.
#### How could the benchmark be improved?
#### What goals should be hit to make this a "good" benchmark?
The dataset is very small and biases towards mutations to small residues. The correlation with all score functions is quite low (~0.3) compared to the correlation for soluble proteins (~0.6). A goal for this benchmark (or rather, protocol development) would be to improve the correlation coefficient to closer to that of soluble proteins. While increasing the size and diversity of the dataset would improve this benchmark, measuring ddG values for membrane proteins is quite difficult, so that will depend on available data.

