Steven Lewis, smlewi@gmail.com, summer 2018

## PURPOSE OF THE TEST
This benchmark runs the entire PDB database, all however many thousands of structures, through a simple Rosetta job.  The fast version just tries to load the PDB or CIF file, essentially exercising just the PDB/CIF reader and SFR->Pose machinery.  The purpose is to check what raw PDBs Rosetta can and cannot handle, and classify what errors it encounters.  This lets us ensure the PDB reading machinery is as robust as possible.

## BENCHMARK DATASET
The raw, whole PDB.  It automatically rsyncs from the PDB before the test tries to run.  Note these are intentionally the raw files.

## PROTOCOL
The protocol itself is trivial, it just loads a PDB and exits.  If if exits cleanly, that PDB passes.  Otherwise, the python wrapper examines the log file and attempts to classify what error occurred.  Errors are grouped in the report to make it easy to find sample errors for bugfixing purposes.

The underlying executable has tools to also exercise simple Rosetta machinery like scoring, repacking, minimization, etc - but those are off in the common test modes.  Once we hit zero errors in the simple "just read the PDB" test we can turn those on; for now they are off because they cost CPU time.

## PERFORMANCE METRICS
For PDBs that crash Rosetta, the test will attempt to classify why they failed.  It does so with string parsing on the log file.  The test is designed to be easy to extend with more error classifications as new bugs crop up, or as the "unknown" catchall category gets too big.

We don't yet expect every PDB to load successfully - a lot of them do stupid stuff.  Accordingly, the test has a manually maintained expected results list and blacklist.  The blacklist is for the case of "the stupid thing this PDB does is that it's huuuuuuuge" and we just skip checking it for time purposes.  The expected results list lets us check how a known-to-fail PDB is expected to fail.

The global test passes or fails based on whether any PDBs failed in a way that doesn't match an expected failure.  So 1000 PDBs crashing is not necessarily a failure if those are all crashes we expect, but 1 unexpected failure will turn the test red.

The blacklist and expected failures lists are meant to be updated manually.  Instructions are at the bottom of the results page where links to the updated copies are found.

## KEY RESULTS
The test sorts out by classification which PDBs failed and why.  You can view lists of all PDBs that failed in a particular category (good for getting samples if you wish to fix the error) and also get log files for each failed case (good for theorizing on the source of the bug).

## LIMITATIONS
This test is not doing any science - it just examines whether or not a PDB could be loaded.
