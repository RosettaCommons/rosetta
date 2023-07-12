# Integration test "esm\_model\_perplexity"

##Author

Moritz Ertelt, Institute for Drug Discovery, University of Leipzig

## Description

This integration test confirms that Rosetta can load the ESM Tensorflow model, predict sequence probabilities and calculate the pseudo-perplexit from them. 

## Expected output

We expect an predicted pseudo-perplexity of ~10. Also outputs an position-specific-scoring-matrix (PSSM) in psi-blast format with the predicted logits.
