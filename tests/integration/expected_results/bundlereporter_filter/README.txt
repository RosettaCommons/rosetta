Integration test for protocols::helical_bundle::BundleReporterFilter filter.  This also tests the RosettaScripts connection. The test creates a mixed alpha/beta structure consisting of an alpha-helical bundle surrounded by an antiparallel beta-barrel "fence".  The reporter filter is then applied to report the helical bundle paramters.  The structure is more or less meaningless, but it confirms the proper functioning of the filter.

Changes to the logfile or to the output pdb mean that the filter has changed somehow.

The three versions of the test check the filter in three different modes: always returning true, always returning false, and filtering by energy.  Only the first should pass and return output.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

