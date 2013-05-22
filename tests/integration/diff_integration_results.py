#!/usr/bin/env python
from integration_test_support import *

import argparse
parser = argparse.ArgumentParser(description="Diff integration tests.")

parser.add_argument("reference_results", help="Reference result directory.")
parser.add_argument("new_results", help="Test result directory.")
parser.add_argument("outdir", help="Output summary directory.")

options = parser.parse_args()

analysis = IntegrationAnalysis(options.reference_results, options.new_results)
analysis.write_test_report(options.outdir)
