Tests FilterReportAsPoseExtraScoresMover

Steven Lewis, Cyrus Biotechnology, smlewi@gmail.com

The Mover's purpose is to allow mid-stream reporting of a Filter's report_sm value, since for some reason those are otherwise only reported at the end of a trajectory.

The test has a pair of residues.  It scores them, moves them apart, then scores them again.

Notice the reported value for just_elec (the filter alone), which matches elec_separated score even though the filter was originally evaluated at the elec_complex state.
