Tests GenMeanFieldMover fixed-backbone packing use case 

Aliza Rubenstein, aliza.rubenstein@gmail.com

The GenMeanFieldMover predicts an amino acid probability distribution (design) or a rotamer probability distribution (sc) for a given set of designed/packed residues (as determined by TaskOperations).  This use case is on a single backbone with packing.

The test inputs include a single protease-peptide complex.  Two residues are packed and then output log, which includes per-rotamer probabilities, is compared for each run. 
