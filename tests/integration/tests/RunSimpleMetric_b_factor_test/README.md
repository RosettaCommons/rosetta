# Integration test "RunSimpleMetric_b_factor_test"

##Author

Moritz Ertelt, Institute for Drug Discovery, Leipzig University

## Description

Tests the metric_to_bfactor option of the RunSimpleMetricsMover, which should set an arbitrary PerResidueRealMetric as temperature factors, which later end up in the b factor column of the pdb.

## Expected output

A line printing whether the output pdb has the right b factor column values.

