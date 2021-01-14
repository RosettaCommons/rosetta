## AUTHOR AND DATE
#### Who set up the benchmark? Please add name, email, PI, month and year

Adapted for the current benchmarking framework by Rocco Moretti (rmorettiase@gmail.com; Meiler Lab), Sep 2018 

## PURPOSE OF THE TEST
#### What does the benchmark test and why?

This benchmark tests how well the enzyme design code is able to recapitulate native-like sequences when run over a set of cocrystal structures of small-molecule binding proteins with their native substrates.

## BENCHMARK DATASET
#### How many proteins are in the set?
#### What dataset are you using? Is it published? If yes, please add a citation.
#### What are the input files? How were the they created?

There are 50 proteins in this set, chosen for being a high-quality structures of proteins binding to their native substrates. This benchmark set (and the basic protocol) is partly described in <a href="https://doi.org/10.1002/prot.24463">Nivon et al. (2014)</a> "Automating human intuition for protein design." The input PDBs are from the previous benchmark tests, so their provenance is not 100% clear, but I believe that they have been downloaded from the RCSBV, minimally cleaned, and subjected to the all-atom constrainted relax protocol of <a href="https://doi.org/10.1371/journal.pone.0059004">Nivon et al. (2013)</a> (probably under the score12/enzdes scorefunction).

## PROTOCOL
#### State and briefly describe the protocol.
#### Is there a publication that describes the protocol?
#### How many CPU hours does this benchmark take approximately?

The protocol follows more-or-less that of <a href="https://doi.org/10.1002/prot.24463">Nivon et al. (2014)</a>, updated for RosettaScripts XML. Briefly, the residues surrounding the ligand (design within 6 Ang (8 if pointed toward ligand) and repack within 10 (12) Ang) are subjected to 2 cycles of softpack/hardmin followed by 1 cycle of hardpack/hardmin. The protein-ligand interactions are upweighted by 1.8-fold for those residues being designed.

One big change from the publication is that the scorefunction being used is the (currently default) REF2015, rather than the older scorefunction used in the paper/in previous benchmarks. (It's the intention that the score function be updated based on whatever the current default is.)

Currently we are only running one output structure for each input, which results in the test taking somewhere around 25-50 CPU hours.

## PERFORMANCE METRICS
#### What are the performance metrics used and why were they chosen?
#### How do you define a pass/fail for this test?
#### How were any cutoffs defined?

The original test only looked at percent sequence recovery at designed positions; whether or not the design recaputulated the identical residue type as the input. ("seqrec")
The concept is that, as the structures are native binders, their current sequence should be close to optimal for binding the ligand.

This reimplementation added two new metrics, based on matching the design result to a PSSM of the input protein. 
This PSSM was generated with the 2.6.0 version of psiblast, using the BLAST nr database from 21May2014
(yes, both the psiblast and nr database were old when this was done in late 2018) with the following command:

    psiblast -query 1A99.fasta -db /path/to/db.21May2014/nr -out_pssm 1A99.chk -out_ascii_pssm 1A99.pssm -save_pssm_after_last_round

The first PSSM-based metric ("pssm_seqrec") is a percent recovery metric, but instead of attempting to match the input sequence exactly, it counts as a "success" matching any amino acid which has a favorable score in the PSSM. This metric is described in <a href="https://doi.org/10.1021/bi200664b">DeLuca et al. 2011</a>. This is a 0-100% normalized scale which better allows for modest changes (e.g. S->T) which are permitted evolutionarily.

The second PSSM-based metric ("pssm_delta_seqrec") looks at the per-residue change in the PSSM score compared to the "native" input sequence. (Where more positive is a "better" design.) This attempts to capture the magnitude of the mutational change.

The benchmark set is summarized by averaging the scores for each protein. (This average is on a per-protein basis, and is not weighted by the number of positions mutated.)

The current tested cutoffs is based only on the match/fail seqrec metric, and is taken directly from the previous iteration of this benchmark. This value is somewhat arbitrary, and was set based on what didn't give overly noisy results when the benchmark was being run.

## KEY RESULTS
#### What is the baseline to compare things to - experimental data or a previous Rosetta protocol?
#### Describe outliers in the dataset. 

I am unaware of any objective level one could compare the results to. Generally you'd be limited to a comparative performance based on prior results.

There's currently no analysis of outliers or per-protein performance, over and above anything discussed in <a href="https://doi.org/10.1002/prot.24463">Nivon et al. (2014)</a>.

## DEFINITIONS AND COMMENTS
#### State anything you think is important for someone else to replicate your results. 

Comparing different runs and different scorefunctions, you're probably looking for something which is consistenty higher.

There's a bit of noise in the runs, so you may need to compare several runs to get better sense of how different scorefunctions compare.

## LIMITATIONS
#### What are the limitations of the benchmark? Consider dataset, quality measures, protocol etc. 
#### How could the benchmark be improved?
#### What goals should be hit to make this a "good" benchmark?

While the dataset attempted to be comprehensive (all structures which fit the criteria) when it was produced, the quality/breadth of structures these days may be better.

The preparation of the input structures may be another issue. The input structures may have been minimized/relax with an older version of the scorefunction. Updating the structure preparation may improve benchmark performance.

The run-to-run variablility is rather high, which might indicate that the benchmark could be improved by running several output structures for each input, rather than just the current one.

The exact protocol being used may not be optimal for ligand binding design, and a different packing/minimization scheme might improve performance.

Finally, keep in mind the intrinsic limitations of sequence-recovery based protocols. While low values are likely bad, you wouldn't necessarily expect the metric to saturate near a "perfect" score, as native proteins are optimized for more than just ligand-binding affinity, so even family PSSM profiles may not match the ideal design that Rosetta would be aiming for.
