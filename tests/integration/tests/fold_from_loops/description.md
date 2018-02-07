# FoldFromLoops Test

This is an integration test for Rosetta's **FoldFromLoops** RS protocol.   
This test gets a small protein with a small binder, selects two regions as _motif_ and re-inserts them into itself.  
Ideally, the folding should be around the same and the disposition of the _motif_ regions in correlation to the binder should not change.  
The test, thus checks that folding multi-segment motifs with the binder (the most complex job type for FFL) actually keeps on working as expected.
Results from this test are not biologically relevant and it only concerns the folding part of the protocol, no design or anything else is carried out.