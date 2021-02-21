## AUTHOR AND DATE

#### Who set up the benchmark? Please add name, email, PI, month and year
Shannon Smith, shannon.t.smith.1@vanderbilt.edu, PI: Jens Meiler, June 2019. 

## PURPOSE OF THE TEST

#### What does the benchmark test and why?
This benchmark tests our ability to correctly correlate experimental binding affinities and computed Rosetta protein-ligand interface scores. This test is meant to use the same dataset and success metrics as those in the Comparative Assessment of Scoring Functions 2016 (CASF). This is a diverse dataset that allows us to compare Rosetta protocols directly to other software suites on the same tests.
This benchmark is meant to be different than the previously-implemented "ligand_docking" benchmark, as this is only looking at the score function instead of performing a docking protocol with full sampling. 

## BENCHMARK DATASET

#### How many proteins are in the set? 
285 native, co-crystal structures of 57 protein-ligand complexes, i.e. 57 proteins (=targets) with varying ligands. The full list with associated logKa values and target information can be found in CoreSet.dat file.

#### What dataset are you using? Is it published? If yes, please add a citation.
CASF-2016 dataset: Comparative Assessment of Scoring Functions: The CASF-2016 Update. (J Chem Inf Model. 2019 Feb 25;59(2):895-913. doi: 10.1021/acs.jcim.8b00545.)

#### What are the input files? How were the they created?
Input files are those directly taken from the CASF dataset. These are the co-crystal structures downloaded directly from the RCSB and ligand params files are generated using molfile_to_params.py

## PROTOCOL

#### State and briefly describe the protocol.
This protocol is implemented in RosettaScripts and takes native protein-ligand complexes, runs a minimization at the interface using the FinalMinimizer mover, then calculates the score using the InterfaceScoreCalculator mover.
Output scores for analysis are taken from the 'interface_delta_X' term without any additional filters/calculations. 

#### Is there a publication that describes the protocol?
A publication describing this benchmark and protocol is available at (Smith Shannon T, Meiler Jens: “Assessing multiple score functions in Rosetta for drug discovery”, PlosOne, 2020, https://doi.org/10.1371/journal.pone.0240450)

#### How many CPU hours does this benchmark take approximately?
This is a fast test. Only takes ~30 seconds max / structure * 285 ~ 2 hours. 

## PERFORMANCE METRICS

#### What are the performance metrics used and why were they chosen?
This test has two components: scoring and ranking. Scoring power refers to a score functions' ability to linearly correlate experimental binding affinities and calculated interface energies in Rosetta and uses a Pearson correlation to determine success. The ranking power test refers to a score functions' ability to correctly rank compounds against a single target and uses a Spearman correlation to determine success.
These metrics were chosen simply by using the same analysis scheme as CASF making these results directly comparable to protocols/score functions used in other softwares. 

#### How do you define a pass/fail for this test?
Values for correlation can range between -1 and 1 where 1 is a perfect correlation. 
We ran this test-set twice to obtain an average Spearman correlation for each target. In each test, we now compare the Spearman correlations of each target to the corresponding value. We are looking for major changes in correlation between runs, for better or worse. Success or failure is based on how many target cases are we deviating more than 0.2 from the standard for each respective target. Overall success is defined as 75% of targets, we maintain within 0.2 of original value. Change these values as necessary. 

#### How were any cutoffs defined?
Cutoffs were defined in the first test run. I should probably do more thorough analysis on previous runs, but looking at the most recent two, correlations between most targets has a standard deviation <0.2.  

## KEY RESULTS

In previous runs, ~75% of targets had Spearman correlations greater than 0.25 - this is not defined as a pass/fail but simply stated.

#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
Previous Rosetta protocols. 

#### Describe outliers in the dataset. 

## DEFINITIONS AND COMMENTS

#### State anything you think is important for someone else to replicate your results. 
Analysis scripts were distributed by CASF, which made it easy to implement here (seriously a lot of copy/paste) meaning the analyses used here can be directly compared to others running the CASF benchmark. 
The given XML uses old ligand score functions (which do work better than the newer ones across all tests), but to assess score functions that are currently being developed, I will likely move this to the new ones in the near future.

## LIMITATIONS

#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc. 
The scoring test is rather vague as this takes the scores of all 285 complexes and determines how linearly correlated these scores are to experimental values. In such a diverse dataset, it is very difficult to capture information about specific interactions that might be outliers. 
In the ranking test, there are clusters where ligands are in different binding pockets. Typically in a drug discovery campaign, groups acquire compound series that are directed towards the same pocket and there are slight chemical changes to tease out SAR information. In some of the ranking tests, those with lower correlations can be attributed to this. It would be nice to have a ranking test where we simply look at similar compound series that greater mimic hit-to-lead optimization efforts. 
The dataset relies strictly on cocrystal structures, which can present issues when working with ligands due to strange geometries, non-ideal bonds or connectivities, etc. 

#### How could the benchmark be improved?
More specific analysis, but probably not here. AKA Which cases do we see consistent good or poor results? Types of interactions that are characterized poorly? Ligand side: functional groups, descriptor classification? Protein side: particular problem residues, protein families, dynamics? It will be good to have this running on a consistent basis to see if these results are consistent. Currently, we only output 1 minimized structure (hopefully at the local minimum). In what I have seen before, results do not change much when we output more decoys, but this is something to keep in mind when looking at the results over time. 
A technical thing: to make sure our analyses are consistent with CASF-published results, I used their scripts with no modification aside from making this compatible with the benchmark server setup. This means that I needed to include several python packages (pandas, sklearn and scipy), which I hope does not cause problems on the server. 

#### What goals should be hit to make this a "good" benchmark?
This is a fast and easy benchmark that can yield information about specific protein-ligand interactions that the score function mischaracterizes. I hope that by running these tests consistently, we can better tease more specific information about where our score functions need improvement. 
Again, I really want to emphasize that provided here are the datasets and analysis scripts that were used in the latest CASF assessment where we can compare these results to >30 other software suites that are popular in computational drug discovery. 
