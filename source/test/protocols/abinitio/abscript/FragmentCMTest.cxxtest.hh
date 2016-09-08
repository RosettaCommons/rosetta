// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <core/id/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/kinematics/Jump.hh>

#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/BBTorsionSRFD.hh>

#include <protocols/environment/Environment.hh>
#include <protocols/environment/EnvExcn.hh>

#include <protocols/filters/Filter.hh>

#include <protocols/simple_moves/FragmentMover.hh>

#include <protocols/abinitio/abscript/FragmentCM.hh>

#include <test/core/init_util.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <utility/tag/Tag.hh>
// #include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <iostream>

std::string const FRAGFILE_LOCATION = "protocols/abinitio/abscript/one_frag3_per_pos";

using namespace core;
using namespace protocols::environment;

class FragmentCMTest : public CxxTest::TestSuite {
public:

	// Shared data elements go here.
	core::pose::Pose pose;

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();

		core::pose::make_pose_from_sequence(pose, "FRIENDLYFRIENDS", "fa_standard");

		for ( Size i = 1; i <= pose.size(); ++i ) {
			pose.set_phi( i, -65 );
			pose.set_psi( i, -41 );
		}

	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_fragment_insertion(){

		core::fragment::FragmentIO frag_io( 1, 1, true );


		protocols::abinitio::abscript::FragmentCMOP fragmover( new protocols::abinitio::abscript::FragmentCM( protocols::simple_moves::FragmentMoverOP( new protocols::simple_moves::ClassicFragmentMover( frag_io.read_data( FRAGFILE_LOCATION ) ) ) ) );

		EnvironmentOP env_op( new Environment( "env" ) );
		Environment & env = *env_op;
		env.register_mover( fragmover );

		core::pose::Pose ppose;
		TS_ASSERT_THROWS_NOTHING( ppose = env.start( pose ) );

		for ( core::Size i = 1; i <= ppose.size(); ++i ) {
			TS_ASSERT_DIFFERS( ppose.phi( i ), pose.phi( i ) );
			TS_ASSERT_DIFFERS( ppose.psi( i ), pose.psi( i ) );
			TS_ASSERT_DIFFERS( ppose.omega( i ) , pose.omega( i ) );
		}

		core::pose::Pose end_pose;
		TS_ASSERT_THROWS_NOTHING( end_pose = env.end( ppose ) );

		for ( core::Size i = 1; i <= ppose.size(); ++i ) {
			if ( i != 1 ) {
				TS_ASSERT_DELTA( ppose.phi( i ), end_pose.phi( i ), 1e-10 );
			}
			if ( i != ppose.size() ) {
				TS_ASSERT_DELTA( ppose.psi( i ), end_pose.psi( i ), 1e-10 );
				TS_ASSERT_DELTA( ppose.omega( i ) , end_pose.omega( i ), 1e-10 );
			}
		}
	}

	void test_chain2_fragment_init( core::pose::Pose & pose ){

		core::select::residue_selector::ChainSelectorOP selector( new core::select::residue_selector::ChainSelector() );
		core::select::residue_selector::ChainSelectorOP bad_selector( new core::select::residue_selector::ChainSelector() );
		selector->set_chain_strings( utility::vector1< std::string >( 1, "2" ) );
		bad_selector->set_chain_strings( utility::vector1< std::string >( 1, "1" ) );

		core::fragment::FragmentIO frag_io( 1, 1, true );

		protocols::abinitio::abscript::FragmentCMOP fragmover( new protocols::abinitio::abscript::FragmentCM( protocols::simple_moves::FragmentMoverOP( new protocols::simple_moves::ClassicFragmentMover( frag_io.read_data( "protocols/abinitio/abscript/one_frag3_per_pos" ) ) ) ) );

		fragmover->set_selector( selector );

		EnvironmentOP env_op( new Environment( "env" ) );
		Environment & env = *env_op;
		env.register_mover( fragmover );

		utility::vector1< bool > selection = selector->apply( pose );
		core::Size const ch2_begin = selection.index( true );
		TS_ASSERT( std::find( selection.begin(), selection.end(), true ) != selection.end() );
		core::Size const ch2_end = pose.size();
		TS_ASSERT( selection[ ch2_end ] );

		core::pose::Pose ppose;
		{
			TS_ASSERT_THROWS_NOTHING( ppose = env.start( pose ) );
			TS_ASSERT_THROWS_ANYTHING( fragmover->set_selector( bad_selector ) );

			for ( core::Size i = 1; i < ch2_begin; ++i ) {
				TS_ASSERT_DELTA( ppose.phi( i ), pose.phi( i ), 1e-10 );
				if ( i != ch2_begin-1 ) {
					TS_ASSERT_DELTA( ppose.psi( i ), pose.psi( i ), 1e-10 );
					TS_ASSERT_DELTA( ppose.omega( i ) , pose.omega( i ), 1e-10 );
				}
			}

			for ( core::Size i = ch2_begin; i <= ch2_end; ++i ) {
				if ( i != ch2_begin ) {
					TS_ASSERT_DIFFERS( ppose.phi( i ), pose.phi( i ) );
				}
				if ( i != ch2_end ) {
					TS_ASSERT_DIFFERS( ppose.psi( i ), pose.psi( i ) );
					TS_ASSERT_DIFFERS( ppose.omega( i ) , pose.omega( i ) );
				}
			}
		}

		core::pose::Pose end_pose;
		TS_ASSERT_THROWS_NOTHING( end_pose = env.end( ppose ) );

		// Verify no change
		for ( core::Size i = 1; i <= end_pose.size(); ++i ) {
			if ( i != 1 && i != ch2_begin ) {
				TS_ASSERT_DELTA( ppose.phi( i ), end_pose.phi( i ), 1e-10 );
			}
			if ( i != ppose.size() && i != ch2_begin-1 ) {
				TS_ASSERT_DELTA( ppose.psi( i ), end_pose.psi( i ), 1e-10 );
				TS_ASSERT_DELTA( ppose.omega( i ) , end_pose.omega( i ), 1e-10 );
			}
		}
	}

	void test_fragment_apply( core::pose::Pose const& pose ) {

		using namespace core::fragment;

		core::Size const frag_start_residue = 1;

		FragmentIO frag_io( 1, 1, true );
		FragSetOP fragset = frag_io.read_data( "protocols/abinitio/abscript/one_frag3_per_pos" );
		FragSetOP onefrag_fragset = fragset->empty_clone();
		onefrag_fragset->add( *(fragset->begin()) );
		TS_ASSERT_EQUALS( onefrag_fragset->size(), frag_start_residue );

		BBTorsionSRFDCOP r1_sfrd = utility::pointer::static_pointer_cast< BBTorsionSRFD const >( onefrag_fragset->begin()->fragment(1).get_residue(1) );

		protocols::abinitio::abscript::FragmentCMOP fragmover(
			new protocols::abinitio::abscript::FragmentCM(
			protocols::simple_moves::FragmentMoverOP( new protocols::simple_moves::ClassicFragmentMover( onefrag_fragset ) ) ) );
		TS_ASSERT_THROWS_NOTHING( fragmover->initialize( false ) );

		protocols::environment::Environment env( "env" );
		env.register_mover( fragmover );

		core::pose::Pose ppose;
		{
			TS_ASSERT_THROWS_NOTHING( ppose = env.start( pose ) );

			// Verify no initialization
			for ( core::Size i = 1; i <= ppose.size(); ++i ) {
				if ( i != 1 ) {
					TS_ASSERT_DELTA( ppose.phi( i ), pose.phi( i ), 1e-10 );
				}
				if ( i != ppose.size() ) {
					TS_ASSERT_DELTA( ppose.psi( i ), pose.psi( i ), 1e-10 );
					TS_ASSERT_DELTA( ppose.omega( i ) , pose.omega( i ), 1e-10 );
				}
			}

			fragmover->apply( ppose );

			for ( core::Size i = onefrag_fragset->min_pos(); i <= onefrag_fragset->max_pos(); ++i ) {
				BBTorsionSRFDCOP sfrd = utility::pointer::static_pointer_cast< BBTorsionSRFD const >( onefrag_fragset->begin()->fragment(1).get_residue( i ) );

				TS_ASSERT_DELTA( ppose.phi( i ), sfrd->torsion( 1 ), 1e-10 );
				TS_ASSERT_DELTA( ppose.psi( i ), sfrd->torsion( 2 ), 1e-10 );
				TS_ASSERT_DELTA( ppose.omega( i ) , sfrd->torsion( 3 ), 1e-10 );
			}
		}

		core::pose::Pose end_pose;
		TS_ASSERT_THROWS_NOTHING( end_pose = env.end( ppose ) );
	}

	void test_xml() {

		protocols::filters::Filters_map f;
		core::pose::Pose tmp;
		basic::datacache::DataMap datamap;

		std::string const SELECTOR = "sele";

		datamap.add( "ResidueSelector",
			SELECTOR,
			core::select::residue_selector::ResidueSelectorOP( new core::select::residue_selector::ResidueIndexSelector( "1-10" ) ) );

		std::stringstream ss;
		ss << "<FragmentCM selector=\"" << SELECTOR << "\" fragments=\""
			<< FRAGFILE_LOCATION << "\"/>";
		utility::tag::TagPtr tag( new utility::tag::Tag );
		tag->read( ss );

		protocols::abinitio::abscript::FragmentCMOP fragmover( new protocols::abinitio::abscript::FragmentCM() );
		TS_ASSERT_THROWS_NOTHING( fragmover->parse_my_tag( tag , datamap, f, protocols::moves::Movers_map(), core::pose::Pose() ) );
		TS_ASSERT_THROWS_NOTHING( tag->die_for_unaccessed_options() );

		EnvironmentOP env_op( new Environment( "env" ) );
		Environment & env = *env_op;

		env.register_mover( fragmover );

		core::pose::Pose ppose;
		TS_ASSERT_THROWS_NOTHING( ppose = env.start( pose ) );

		for ( core::Size i = 1; i <= ppose.size(); ++i ) {
			TS_ASSERT_DIFFERS( ppose.phi( i ), pose.phi( i ) );
			TS_ASSERT_DIFFERS( ppose.psi( i ), pose.psi( i ) );
			TS_ASSERT_DIFFERS( ppose.omega( i ) , pose.omega( i ) );
		}

		core::pose::Pose end_pose;
		TS_ASSERT_THROWS_NOTHING( end_pose = env.end( ppose ) );
	}
};
