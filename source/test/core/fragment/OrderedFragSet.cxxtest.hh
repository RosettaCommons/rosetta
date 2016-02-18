// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/RotamerTrials.cxxtest.hh
/// @brief  test suite for rotamer_trials
/// @author Phil Bradley
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/util.hh>


#include <core/types.hh>

#include <basic/Tracer.hh>
#include <numeric/numeric.functions.hh>

//Auto Headers
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/id/AtomID_Mask.hh>
#include <utility/fix_boinc_read.hh>
#include <utility/vector1.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer tr("core.fragment.OrderedFragSet.cxxtest");
//MY_TRACERS("core.fragment.ConstantLengthFragments.cxxtest")

using namespace core;
using namespace fragment;
using namespace ObjexxFCL;

// hacky test code ---  should live somewhere else
// as defined in ConstantLengthFragments.cxxtest.hh
//extern void steal_constant_length_frag_set_from_pose ( pose::Pose const& pose, ConstantLengthFragSet& fragset );

class OrderedFragSetTest : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCOP residue_set_;
	pose::Pose pose_random_, pose_;
public:
	OrderedFragSetTest() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
		residue_set_ = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID );

		//core::import_pose::pose_from_file( pose_, "core/pack/test_in.pdb" , core::import_pose::PDB_file);
		pose_ = create_test_in_pdb_pose();

	}
	//helpers
	void generate_random_pose();

	// tests
	void test_frag_iterator();

	// Shared finalization goes here.
	void tearDown() {
	}

private:

};

void OrderedFragSetTest::generate_random_pose() {
	using namespace chemical;
	using namespace conformation;
	std::string sequence = pose_.sequence();
	//create extended pose from sequence
	for ( Size pos = 1; pos <= sequence.size(); pos++ ) {
		chemical::AA aa = aa_from_oneletter_code( sequence[ pos-1 ] );
		// Really, this could be changed to get_representative_type
		// But let's keep the unit test as-is
		ResidueTypeCOPs res_list = ResidueTypeFinder( *residue_set_ ).name3( name_from_aa ( aa ) ).get_all_possible_residue_types();
		ResidueOP new_rsd( ResidueFactory::create_residue( * ( res_list[ 1 ] ) ) );
		pose_random_.append_residue_by_bond( *new_rsd );
	}
	io::pdb::dump_pdb( pose_random_, "random_chain_pose" );
}


void OrderedFragSetTest::test_frag_iterator() {
	using namespace pose;
	using namespace fragment;
	kinematics::MoveMap movemap;
	movemap.set_bb( true );

	Size len9(9);
	ConstantLengthFragSet fragset9( len9 );
	steal_constant_length_frag_set_from_pose ( pose_, fragset9 );

	ConstantLengthFragSet fragset3( 3 );
	steal_constant_length_frag_set_from_pose ( pose_, fragset3 );

	OrderedFragSet fragset;

	// consolidate both fragsets in new OrderedFragSet
	for ( ConstFrameIterator it = fragset3.begin(), eit = fragset3.end();
			it!=eit; ++it ) {
		fragset.add( *it );
	}

	for ( ConstFrameIterator it = fragset9.begin(), eit = fragset9.end();
			it!=eit; ++it ) {
		fragset.add( *it );
	}

	// now there should be two frames with the same starting positions until we run out of
	// 9mers...  let's check this.
	Size ct = 1;
	for ( ConstFrameIterator it = fragset.begin(), eit = fragset.end();
			it!=eit; ++it ) {
		//    tr.Info << " ct: " << ct << " " << (*it)->start() << "\n";
		TS_ASSERT( (*it)->start() == ct );
		if ( ct <= pose_.total_residue() - 8 ) ++it;
		TS_ASSERT( (*it)->start() == ct );
		ct++;
	}

	// a call to region should return us two Frames per position
	FrameList frames;
	fragset.region( movemap, 5, 5, 0, 0, frames);
	TS_ASSERT_EQUALS( frames.size() , 2 );

	// a call to region should return us two Frames per position
	frames.clear();
	fragset.region( movemap, 5, 6, 0, 0, frames);

	TS_ASSERT_EQUALS( frames.size(), 4 );

	//let's check the clone function
	FragSetOP cloned = fragset.clone();
	ConstFrameIterator cit = cloned->begin();
	ConstFrameIterator ecit = cloned->end();
	for ( ConstFrameIterator it = fragset.begin(), eit = fragset.end();
			it!=eit; ++it,++cit ) {
		//    tr.Info << " ct: " << ct << " " << (*it)->start() << "\n";
		TS_ASSERT( cit != ecit );
		if ( !( cit!=ecit ) ) break; //avoid run-time errors
		TS_ASSERT( (*it)->start() == (*cit)->start() );
	}
}
