# Integration test "sample\_seq\_from\_probs"

##Author

Moritz Ertelt, Institute for Drug Discovery, University of Leipzig

## Description

This integration test confirms that Rosetta can load and save probabilities to/from a weights file, and sample sequences from the loaded probabilities.

## Expected output
A weights file that matches the input weights file. Also should change only the first residue of 1ubq to ALA instead of MET.
