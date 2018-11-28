## AUTHOR AND DATE
#### Andrew Watkins (andy.watkins2@gmail.com); Rhiju Das (October 2018)

## PURPOSE OF THE TEST
#### This test judges the stability of the stepwise Monte Carlo algorithm's performance.

## BENCHMARK DATASET
#### The benchmark set contains the 12 loop modeling problems from favorites.txt, a benchmark first established in Watkins et al., Science Advances 2018. The input files are ideal A-form RNA helices, fasta files, and portions of crystal structures.

## PROTOCOL
#### See Watkins et al., Science Advances 2018 for a description of the SWM protocol. The benchmark takes 6 hours on the test server (~240 CPU-hours).

## PERFORMANCE METRICS
#### Benchmark performance is based on the minimum RMSD sampled and minimum energy sampled; each case has a separate threshold. Since sampling is limited for each test, cutoffs are determined conservatively through repeated runs..

## KEY RESULTS
#### n/a

## DEFINITIONS AND COMMENTS
#### n/a

## LIMITATIONS
#### The benchamrk could be expanded to include favorites2.txt, challenges.txt, followups.txt, and more of the benchmarking challenges already developed for stepwise Monte Carlo. More CPU power would permit better sampling and therefore more precise cutoffs.

