// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_ddg/DdGScan.cxxtest.hh
/// @brief  test for DdGScan filter
/// @author Kyle Barlow (kb@kylebarlow.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <protocols/simple_ddg/DdGScan.hh>
#include <protocols/simple_ddg/ddG.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/file/FileName.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.simple_filters.DdGScan.cxxtest.hh");

// ------------ Mock ddGMover Class --------- //

class MockddG : public protocols::simple_ddg::ddG {

public:
	MockddG() {
		sf_ = StubMultiFilterOP( new StubMultiFilter( false ) );
		sf_->push_back( 7.0 );
		sf_->push_back( 107.0 );
	}

	void scorefxn( core::scoring::ScoreFunctionCOP /* scorefxn_in */ ) {}

	void calculate( Pose const & /* pose_in */ ) {}

	Real sum_ddG() const {
		return sf_->get_next_value();
	}

	void report_ddG( std::ostream & /* out */ ) const {}

private:

	StubMultiFilterOP sf_;

};

// --------------- Test Class --------------- //

class DdGScan : public CxxTest::TestSuite {

private:

	core::pose::PoseOP test_dimer_pose_;
	core::pose::PoseOP test_monomer_pose_;
	core::scoring::ScoreFunctionOP scorefxn_;
	protocols::simple_ddg::ddGOP ddg_mover_;
	protocols::simple_ddg::DdGScanOP ddg_scan_mover_;
	core::pack::task::TaskFactoryOP tf_;

public:

	void setUp() {
		core_init();

		test_dimer_pose_ = create_2res_1ten_2res_trp_cage_poseop(); //dimer structure
		// test_monomer_pose_ = create_twores_1ubq_poseop();
		scorefxn_ = scorefxn_ = core::scoring::get_score_function();

		ddg_mover_ = protocols::simple_ddg::ddGOP( new MockddG() );
		ddg_scan_mover_ = protocols::simple_ddg::DdGScanOP( new protocols::simple_ddg::DdGScan() );
		ddg_scan_mover_->scorefxn( scorefxn_ );
		ddg_scan_mover_->repeats( 1 );

		tf_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory );

		core::pack::task::operation::PreventRepackingOP prevent_repack_taskop(
			new core::pack::task::operation::PreventRepacking()
		);

		core::pack::task::operation::RestrictAbsentCanonicalAASOP design_taskop(
			new core::pack::task::operation::RestrictAbsentCanonicalAAS()
		);

		for ( core::Size i = 1; i <= test_dimer_pose_->size(); ++i ) {
			if ( i == 2 ) {
				design_taskop->include_residue(i);
				design_taskop->keep_aas("A");
			} else {
				prevent_repack_taskop->include_residue( i );
			}
		}

		tf_->push_back( prevent_repack_taskop );
		tf_->push_back( design_taskop );

	}

	void tearDown() {
	}

	void test_calculate() {
		(*scorefxn_)(*test_dimer_pose_);

		ddg_scan_mover_->task_factory( tf_ );
		ddg_scan_mover_->ddG_mover( ddg_mover_ );
		utility::vector1< ddG_data_tuple > ddG_data = ddg_scan_mover_->calculate( TR, *test_dimer_pose_ );

		core::Size resNum; std::string resname; core::Real ddG_value;
		for ( utility::vector1<ddG_data_tuple>::const_iterator iter = ddG_data.begin(), iter_end = ddG_data.end() ; iter != iter_end ; ++iter ) {
			boost::tie(resNum, resname, ddG_value) = *iter;
			TS_ASSERT_EQUALS( resNum, 2 );
			TS_ASSERT_EQUALS( resname, "ALA" );
			TS_ASSERT_EQUALS( ddG_value, 100 );
		}
	}

	void test_filter_parsing() {
		basic::datacache::DataMap data;
		Filters_map filters;
		Movers_map movers;

		prime_Data( data );
		movers["ddg_mover_name"] = ddg_mover_;

		(*scorefxn_)(*test_dimer_pose_);

		TagCOP tag = tagptr_from_string("<DdGScan name=test ddG_mover=ddg_mover_name />\n");
		ddg_scan_mover_->parse_my_tag( tag, data, filters, movers, *test_dimer_pose_ );
		ddg_scan_mover_->task_factory( tf_ );

		utility::vector1< ddG_data_tuple > ddG_data = ddg_scan_mover_->calculate( TR, *test_dimer_pose_ );
		core::Size resNum; std::string resname; core::Real ddG_value;
		for ( utility::vector1<ddG_data_tuple>::const_iterator iter = ddG_data.begin(), iter_end = ddG_data.end() ; iter != iter_end ; ++iter ) {
			boost::tie(resNum, resname, ddG_value) = *iter;
			TS_ASSERT_EQUALS( resNum, 2 );
			TS_ASSERT_EQUALS( resname, "ALA" );
			TS_ASSERT_EQUALS( ddG_value, 100 );
		}
	}

};
