// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/SwitchChainOrderMover.cxxtest.hh
/// @brief  test suite for SwitchChainOrderMover
/// @author Brian Coventry (bcov@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <core/types.hh>

#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/rosetta_scripts/XmlObjects.hh>
#include <protocols/simple_moves/SwitchChainOrderMover.hh>

#include <basic/datacache/DataMap.hh>

#include <utility/tag/Tag.hh>


#include <test/UTracer.hh>
#include <test/core/init_util.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

using namespace core;
using namespace protocols::simple_moves;


static basic::Tracer TR("protocols.simple_moves.SwitchChainOrderMover.cxxtest");

class SwitchChainOrderMoverTests : public CxxTest::TestSuite
{
public:
	void setUp()
	{
		core_init();
		import_pose::pose_from_file( pose_, "protocols/membrane/1AFO_AB.pdb" );
	}



	void test_pdbinfo_labels(){
		pose::Pose pose = pose_;

		conformation::Conformation const & conf = pose.conformation();

		TS_ASSERT( conf.num_chains() == 2 );
		TS_ASSERT( conf.chain_end(1) == 40 );
		TS_ASSERT( conf.chain_end(2) == 80 );

		pose.pdb_info()->add_reslabel(1,  "RES1");
		pose.pdb_info()->add_reslabel(2,  "RES2");
		pose.pdb_info()->add_reslabel(40, "RES40");
		pose.pdb_info()->add_reslabel(41, "RES41");
		pose.pdb_info()->add_reslabel(42, "RES42");
		pose.pdb_info()->add_reslabel(80, "RES80");


		// It grabs one of the default scorefunctions so we have to instantiate it this way
		protocols::moves::MoverOP switch_chain = protocols::rosetta_scripts::XmlObjects::static_get_mover(
			"<SwitchChainOrder name=\"swap\" chain_order=\"21\" />"
		);
		// std::string tag_string = ;
		// std::stringstream ss( tag_string );
		// utility::tag::TagOP tag( new utility::tag::Tag() );
		// tag->read( ss );
		// basic::datacache::DataMap dm;
		// protocols::filters::Filters_map fm;
		// protocols::moves::Movers_map mm;

		// SwitchChainOrderMover switch_chain;
		// switch_chain.parse_my_tag( tag, dm, fm, mm, pose );

		switch_chain->apply( pose );

		TS_ASSERT( pose.pdb_info()->res_haslabel( 1,  "RES41" ));
		TS_ASSERT( pose.pdb_info()->res_haslabel( 2,  "RES42" ));
		TS_ASSERT( pose.pdb_info()->res_haslabel( 40, "RES80" ));
		TS_ASSERT( pose.pdb_info()->res_haslabel( 41, "RES1" ));
		TS_ASSERT( pose.pdb_info()->res_haslabel( 42, "RES2" ));
		TS_ASSERT( pose.pdb_info()->res_haslabel( 80, "RES40" ));

	}


private:
	pose::Pose pose_;


};
