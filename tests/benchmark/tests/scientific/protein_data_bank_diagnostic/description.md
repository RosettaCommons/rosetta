This test takes the whole PDB (all 140K+ files) and tries to load them into Rosetta.

The default "fast" version of the test is JUST checking if they read in successfully.  Fancier versions will check if they can be scored, packed, and minimized.

If the run does not crash - great!  If the run crashes - the test saves the log file and attempts to classify the error based on the text of the error message, so that the most common errors can be identified.  

Once the test has run on the test server, a results page index.html will show up (that's not a link) which will give more details on interpretation.
