// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/analysis/simple_metrics/ConstraintsMetricTests.cxxtest.hh
/// @brief  Unit tests for the ConstraintsMetric, which dumps out the constraints in a pose or selection.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/analysis/simple_metrics/ConstraintsMetric.hh>

// Protols Headers
#include <protocols/constraint_movers/ConstraintSetMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/import_pose/import_pose.hh>

// Utility Headers
#include <utility/string_util.hh>

// Basic Headers
#include <basic/Tracer.hh>

// C++ Headers
#include <sstream>

static basic::Tracer TR("ConstraintsMetricTests");


class ConstraintsMetricTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}

	void tearDown() {
	}

	/// @brief Are two values equal within error?
	bool equal_within_thresh( core::Real const val1, core::Real const val2 ) {
		return std::abs( val1-val2 ) < 1.0e-4;
	}

	/// @brief Print a failure message and return false.
	bool failure( std::string const & line1, std::string const & line2 ) {
		TR << "\"" << line1 << "\" != \"" << line2 << "\"" << std::endl;
		return false;
	}

	/// @brief Check that (whitespace aside) two lines match.
	bool lines_match(
		std::string const & line1,
		std::string const & line2
	) {
		utility::vector1< std::string > const line1_split( utility::split_whitespace(line1) );
		utility::vector1< std::string > const line2_split( utility::split_whitespace(line2) );
		if ( line1_split.size() != line2_split.size() ) return failure(line1, line2);

		for ( core::Size i(1), imax(line1_split.size()); i<=imax; ++i ) {
			{ //Try floats first:
				std::istringstream ss1( line1_split[i] ), ss2( line2_split[i] );
				core::Real val1, val2;
				ss1 >> val1;
				ss2 >> val2;
				if ( !ss1.fail() ) {
					if ( ss2.fail() || !equal_within_thresh( val1, val2 ) ) return failure( line1, line2 );
					continue;
				} else {
					if ( !ss2.fail() ) return failure( line1, line2 );
				}
			}

			{ //Try ints first:
				std::istringstream ss1( line1_split[i] ), ss2( line2_split[i] );
				signed long val1, val2;
				ss1 >> val1;
				ss2 >> val2;
				if ( !ss1.fail() ) {
					if ( ss2.fail() || (val1 != val2) ) return failure( line1, line2 );
					continue;
				} else {
					if ( !ss2.fail() ) return failure( line1, line2 );
				}
			}

			{ //Just compare strings next:
				if ( utility::upper(line1_split[i]) == "NONE" ) {
					if ( utility::upper(line2_split[i]) != "NONE" ) return failure( line1, line2 );
				} else {
					if ( line1_split[i] != line2_split[i] ) return failure( line1, line2 );
				}
			}
		}

		return true;
	}

	/// @brief Remove empty strings from a list of strings.
	void remove_empty(
		utility::vector1< std::string > & list
	) {
		TR << "----------" << std::endl;
		for ( utility::vector1<std::string>::iterator it(list.begin()); it != list.end(); ) {
			(*it) = utility::strip_whitespace( *it );
			TR << "Processing " << *it;
			if ( it->empty() ) {
				TR << "\t-- REMOVED" << std::endl;
				it = list.erase( it );
			} else {
				TR << "\t-- KEPT" << std::endl;
				++it;
			}
		}
		TR << "----------" << std::endl;
	}

	/// @brief Read in a constraints file, apply it to a pose, and then apply the ConstraintsMetric.
	/// Check whether the output from the metric is the same as the original constraints file.
	void test_read_and_summarize() {
		std::string const file_contents( utility::file_contents( "protocols/analysis/simple_metrics/trp_cage_constraints.cst" ) );
		utility::vector1< std::string > input_lines( utility::string_split( file_contents, '\n' ) );
		TR << "Read input file:\n" << file_contents << std::endl;
		TR << "Detected " << input_lines.size() << " input lines (including 1 blank)." << std::endl;
		core::pose::Pose pose( fullatom_pose_from_string( trp_cage_ideal() ) );
		protocols::constraint_movers::ConstraintSetMover cst_adder;
		cst_adder.set_cst_fa_file( "protocols/analysis/simple_metrics/trp_cage_constraints.cst" );
		cst_adder.apply( pose );
		TS_ASSERT( pose.constraint_set()->get_all_constraints().size() == input_lines.size() - 1 );

		//Apply the metric:
		protocols::analysis::simple_metrics::ConstraintsMetric cst_metric;
		std::string const metric_output( cst_metric.calculate(pose) );
		TR << "Metric output was:\n" << metric_output << std::endl;

		utility::vector1< std::string > output_lines( utility::string_split( metric_output, '\n' ) );
		TR << "Detected " << output_lines.size() << " output lines (including 2 blank)." << std::endl;
		remove_empty( input_lines );
		remove_empty( output_lines );

		TS_ASSERT_EQUALS( input_lines.size(), output_lines.size() );
		for ( core::Size i(1); i<=std::min(input_lines.size(), output_lines.size()); ++i ) {
			TS_ASSERT( lines_match( input_lines[i], output_lines[i] ) );
		}
	}


};
