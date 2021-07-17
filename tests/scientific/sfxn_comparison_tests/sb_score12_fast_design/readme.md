## AUTHOR AND DATE
#### Who set up the benchmark? Please add name, email, PI, month and year
This benchmark was set up by Julia Koehler Leman (julia.koehler.leman@gmail.com), PI Richard Bonneau, March 2020
Input data and command lines are from Chris Bahl and Jack Maguire. 

## PURPOSE OF THE TEST
#### What does the benchmark test and why?
The benchmark tests how well FastDesign can recover native sequences on the benchmark set.

## BENCHMARK DATASET
#### How many proteins are in the set?
#### What dataset are you using? Is it published? If yes, please add a citation.
#### What are the input files? How were the they created?
The benchmark set contains 48 proteins between 102 and 176 residues, originally used by Frank DiMaio for his improvements to the energy function. The set covers alpha-helical bundles, beta-sheets proteins and mixed alpha/beta folds.

## PROTOCOL
#### State and briefly describe the protocol.
#### Is there a publication that describes the protocol?
#### How many CPU hours does this benchmark take approximately?
The protocol runs FastDesign in RosettaScripts currently with 1 iteration, nstruct 100, no extrachi. Probably should try 5 iterations as originally suggested. 1 iteration generates a decoy in about 2000 seconds. This makes this protocol run for about 48 x 100 x 2000 / 3600 = 2666 CPU hours.

## PERFORMANCE METRICS
#### What are the performance metrics used and why were they chosen?
#### How do you define a pass/fail for this test?
#### How were any cutoffs defined?
We use sequence recovery between the native and the design computed via SimpleMetrics in RosettaScripts. The cutoffs were defined for sequence recovery, for each protein take the minimum minus 2 stdev. For the score, per protein take the maximum plus 5 stdev. For the score12 version we use 90% of the score values below the cutoff because there are occasional outliers with high scores. 

## KEY RESULTS
#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
#### Describe outliers in the dataset. 
The sequence recovery metric has been used for many years to benchmark design applications. Historically, sequence recoveries are somewhere between 30% and 60% at the maximum. It is difficult for the scorefunction to recapitulate native sequences accurately. It might be worth noting that we do not expect 100% sequence recovery even with a "perfect" energy function and "perfect" optimizer, since evolution optimizes proteins for marginal stability (to allow for degradation) and for other things (function, genetic code, amino acid costs/abundances), while we're trying to optimize for high stability (and maximize the stability of the designed state, without knowing what we're doing to the stability of alternative conformations).

## DEFINITIONS AND COMMENTS
#### State anything you think is important for someone else to replicate your results. 

## LIMITATIONS
#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc. 
#### How could the benchmark be improved?
#### What goals should be hit to make this a "good" benchmark?
The benchmark set only consists of small, soluble proteins. It would be good to know how design performs on larger proteins and more complex folds. For the quality metrics, rotamer recovery could be considered as well. 
