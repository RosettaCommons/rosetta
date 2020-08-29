// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/PeptideInternalHbondsMetricTests.cxxtest.hh
/// @brief  Unit tests for the PeptideInternalHbondsMetric and its wrapper filter.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/cyclic_peptide/PeptideInternalHbondsMetric.hh>
#include <protocols/cyclic_peptide/PeptideInternalHbondsFilter.hh>

// Protocols Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("PeptideInternalHbondsMetricTests");


class PeptideInternalHbondsMetricTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	/// @brief Helper function.
	bool
	list_has_value(
		std::list< core::Size > const & mylist,
		core::Size const myval
	) const {
		return ( std::find( mylist.begin(), mylist.end(), myval ) != mylist.end() );
	}

	void test_peptide_internal_hbonds_filter() {
		TR << "Starting PeptideInternalHbondsMetricTests:test_peptide_internal_hbonds_filter." << std::endl;
		using namespace protocols::cyclic_peptide;

		core::pose::PoseOP pose( core::import_pose::pose_from_file( "protocols/cyclic_peptide/sample_cycpep_8mer.pdb" ) );
		DeclareBond declbond;
		declbond.set( 8, "C", 1, "N", false );
		declbond.apply(*pose);

		PeptideInternalHbondsFilter filter;

		filter.set_exclusion_distance(0);
		filter.set_hbond_cutoff(5);
		TS_ASSERT_EQUALS( filter.apply(*pose), true );
		TS_ASSERT_EQUALS( filter.report_sm(*pose), 5.0 );
		filter.set_hbond_cutoff(6);
		TS_ASSERT_EQUALS( filter.apply(*pose), false );
		TS_ASSERT_EQUALS( filter.report_sm(*pose), 5.0 );

		filter.set_exclusion_distance(1);
		filter.set_hbond_cutoff(5);
		TS_ASSERT_EQUALS( filter.apply(*pose), true );
		TS_ASSERT_EQUALS( filter.report_sm(*pose), 5.0 );
		filter.set_hbond_cutoff(6);
		TS_ASSERT_EQUALS( filter.apply(*pose), false );
		TS_ASSERT_EQUALS( filter.report_sm(*pose), 5.0 );

		filter.set_exclusion_distance(2);
		filter.set_hbond_cutoff(3);
		TS_ASSERT_EQUALS( filter.apply(*pose), true );
		TS_ASSERT_EQUALS( filter.report_sm(*pose), 3.0 );
		filter.set_hbond_cutoff(4);
		TS_ASSERT_EQUALS( filter.apply(*pose), false );
		TS_ASSERT_EQUALS( filter.report_sm(*pose), 3.0 );

		TR << "Finished PeptideInternalHbondsMetricTests:test_peptide_internal_hbonds_filter." << std::endl;
	}

	void test_peptide_internal_hbonds_metric() {
		TR << "Starting PeptideInternalHbondsMetricTests:test_peptide_internal_hbonds_metric." << std::endl;
		using namespace protocols::cyclic_peptide;

		core::pose::PoseOP pose( core::import_pose::pose_from_file( "protocols/cyclic_peptide/sample_cycpep_8mer.pdb" ) );
		DeclareBond declbond;
		declbond.set( 8, "C", 1, "N", false );
		declbond.apply(*pose);

		PeptideInternalHbondsMetric metric;

		metric.set_exclusion_distance(0);
		TS_ASSERT_EQUALS( metric.calculate(*pose), 5 );

		metric.set_exclusion_distance(1);
		TS_ASSERT_EQUALS( metric.calculate(*pose), 5 );

		metric.set_exclusion_distance(2);
		TS_ASSERT_EQUALS( metric.calculate(*pose), 3 );
		TR << "Finished PeptideInternalHbondsMetricTests:test_peptide_internal_hbonds_metric." << std::endl;
	}

	void test_generate_allowed_partners(){
		TR << "Starting PeptideInternalHbondsMetricTests:test_generate_allowed_partners." << std::endl;
		using namespace protocols::cyclic_peptide;
		PeptideStubMover builder;
		builder.add_residue( "Append", "ALA", 1, true, "", 0, 0, nullptr, "" );
		builder.add_residue( "Append", "ALA", 2, false, "", 0, 0, nullptr, "" );
		builder.add_residue( "Append", "ALA", 3, false, "", 0, 0, nullptr, "" );
		builder.add_residue( "Append", "ALA", 4, false, "", 0, 0, nullptr, "" );
		builder.add_residue( "Append", "ALA", 5, false, "", 0, 0, nullptr, "" );
		builder.add_residue( "Append", "ALA", 6, false, "", 0, 0, nullptr, "" );
		builder.add_residue( "Append", "ALA", 7, false, "", 0, 0, nullptr, "" );
		builder.add_residue( "Append", "ALA", 8, false, "", 0, 0, nullptr, "" );
		builder.add_residue( "Append", "ALA", 9, false, "", 0, 0, nullptr, "" );
		builder.add_residue( "Append", "ALA", 10, false, "", 0, 0, nullptr, "" );
		core::pose::Pose pose;
		builder.apply(pose);

		DeclareBond declbond;
		declbond.set( 10, "C", 1, "N", false );
		declbond.apply(pose);

		TS_ASSERT_EQUALS( pose.total_residue(), 10 );

		//Test the generate_allowed_partners() function with exclusion distance of zero:
		PeptideInternalHbondsMetric metric;
		metric.set_exclusion_distance(0);
		utility::vector1< core::Size > allowed_res(10);
		for ( core::Size i(1); i<=10; ++i ) { allowed_res[i] = i; }
		utility::vector1< std::list< core::Size > > const allowed0( metric.generate_allowed_partners( pose, allowed_res ) );
		for ( core::Size i(1); i<=10; ++i ) {
			TS_ASSERT_EQUALS( allowed0[i].size(), 9 );
			for ( core::Size j(1); j<=10; ++j ) {
				if ( i == j ) {
					TS_ASSERT( !list_has_value(allowed0[i], j) );
				} else {
					TS_ASSERT( list_has_value(allowed0[i], j) );
				}
			}
		}

		//Test the generate_allowed_partners() function with exclusion distance of one:
		metric.set_exclusion_distance(1);
		utility::vector1< std::list< core::Size > > const allowed1( metric.generate_allowed_partners( pose, allowed_res ) );
		for ( core::Size i(1); i<=10; ++i ) {
			TS_ASSERT_EQUALS( allowed1[i].size(), 7 );
			for ( core::Size j(1); j<=10; ++j ) {
				if ( i == j || i == (j-1) || i == (j+1) || (i == 1 && j == 10) || (i == 10 && j == 1) ) {
					TS_ASSERT( !list_has_value(allowed1[i], j) );
				} else {
					TS_ASSERT( list_has_value(allowed1[i], j) );
				}
			}
		}

		//Test the generate_allowed_partners() function with exclusion distance of two:
		metric.set_exclusion_distance(2);
		utility::vector1< std::list< core::Size > > const allowed2( metric.generate_allowed_partners( pose, allowed_res ) );
		for ( core::Size i(1); i<=10; ++i ) {
			TS_ASSERT_EQUALS( allowed2[i].size(), 5 );
			for ( core::Size j(1); j<=10; ++j ) {
				if ( i == j || i == (j-1) || i == (j-2) || i == (j+1) || i == (j+2) ||
						(i == 1 && (j == 10 || j == 9 )) || (i == 10 && (j == 1 || j == 2)) ||
						(i == 2 && (j == 10 || j == 1 )) || (i == 9 && (j == 10 || j == 1))
						) {
					TS_ASSERT( !list_has_value(allowed2[i], j) );
				} else {
					TS_ASSERT( list_has_value(allowed2[i], j) );
				}
			}
		}
		TR << "Finished PeptideInternalHbondsMetricTests:test_generate_allowed_partners." << std::endl;
	}



};
