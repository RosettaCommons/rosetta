// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/denovo_design/ConsensusLoopDesignOperation.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::components::ConsensusLoopDesignOperation
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/task_operations/ConsensusLoopDesign.hh>

// Protocol headers
#include <protocols/denovo_design/util.hh>

// Core headers
#include <core/io/pdb/file_data.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/sequence/ABEGOManager.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// C++ headers

static THREAD_LOCAL basic::Tracer TR("protocols.denovo_design.ConsensusLoopDesignOperation.cxxtest");

// --------------- Test Class --------------- //
class ConsensusLoopDesignOperationTests : public CxxTest::TestSuite {
public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init();

		// set preserve header always
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_loop_info_and_task_manipulation() {
		using namespace protocols::denovo_design::task_operations;
		core::pose::Pose input_pose;
		// loop is 17-18
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/disulf_test.pdb" );

		ConsensusLoopDesignOperation loopdesign;
		std::string const loop_abego = "GB";
		SurroundingSS ss_around( 'E', 'H' );

		LoopAAs allowed_aas;
		allowed_aas.push_back( "G" );
		allowed_aas.push_back( "DNS" );
		loopdesign.set_allowed_aas( ss_around, loop_abego, allowed_aas );

		LoopAAs const & aalist = loopdesign.allowed_aas( ss_around, loop_abego );
		TS_ASSERT_EQUALS( aalist.size(), allowed_aas.size() );

		std::list< std::string >::const_iterator b = allowed_aas.begin();
		for ( LoopAAs::const_iterator a = aalist.begin(), enda = aalist.end(); a != enda; ++a, ++b ) {
			TS_ASSERT_EQUALS( *a, *b );
		}

		core::pack::task::TaskFactory tf;
		core::pack::task::PackerTaskOP task = tf.create_task_and_apply_taskoperations( input_pose );
		LoopInfo info;
		info.startres = 17;
		info.ss_around = ss_around;
		info.abego = loop_abego;
		loopdesign.disallow_aas( *task, info );

		TS_ASSERT_EQUALS( task->total_residue(), input_pose.total_residue() );
		for ( core::Size i=1, endi=input_pose.total_residue(); i<=endi; ++i ) {
			// all aas should be allowed except for residues 17, 18
			TS_ASSERT( task->being_packed( i ) );
			if ( i == 17 ) {
				TS_ASSERT_EQUALS( task->residue_task( i ).allowed_residue_types().size(), 1 );
				for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter t=task->residue_task( i ).allowed_residue_types_begin(), endt = task->residue_task( i ).allowed_residue_types_end(); t != endt; ++t ) {
					TS_ASSERT_EQUALS( (*t)->name1(), 'G' );
				}
			} else if ( i == 18 ) {
				std::string const & c = *(allowed_aas.rbegin());
				TS_ASSERT_EQUALS( task->residue_task( i ).allowed_residue_types().size(), c.size() );
				for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter t=task->residue_task( i ).allowed_residue_types_begin(), endt = task->residue_task( i ).allowed_residue_types_end(); t != endt; ++t ) {
					TS_ASSERT( c.find( (*t)->name1() ) != std::string::npos );
				}
			} else if ( input_pose.residue( i ).name() == "CYS:disulfide" ) {
				TS_ASSERT_EQUALS( task->residue_task( i ).allowed_residue_types().size(), 1 );
			} else {
				// the +1 is for HIS_D
				TS_ASSERT_EQUALS( task->residue_task( i ).allowed_residue_types().size(), core::Size( core::chemical::num_canonical_aas )+1 );
			}
		}

		// test getting loop info
		LoopInfoVec loopinfo = loopdesign.get_loop_info( input_pose );
		TS_ASSERT_EQUALS( loopinfo.size(), 4 );

		// get ss/abego for comparision
		core::scoring::dssp::Dssp dssp( input_pose );
		std::string const ss = dssp.get_dssp_secstruct();
		utility::vector1< std::string > const abego = core::sequence::get_abego( input_pose, 1 );

		for ( LoopInfoVec::const_iterator l=loopinfo.begin(), endl=loopinfo.end(); l != endl; ++l ) {
			TS_ASSERT_EQUALS( ss[ l->startres - 1 ], l->ss_around.before );
			TS_ASSERT_EQUALS( ss[ l->startres + l->abego.size() + 1 ], l->ss_around.after );
			for ( core::Size i=l->startres+1, endi=l->startres+l->abego.size()-1; i<endi; ++i ) {
				core::Size const abegoidx = i - l->startres;
				TS_ASSERT_EQUALS( ss[ i - 1 ], 'L' );
				TS_ASSERT_EQUALS( abego[ i ].size(), 1 );
				TS_ASSERT_EQUALS( abego[ i ][ 0 ], l->abego[ abegoidx ] );
			}
		}

	}

	void test_read_db() {
		using namespace protocols::denovo_design::task_operations;
		ConsensusLoopDesignOperation loopdesign;

		// get stats from a known motif -- GB connecting two helices
		LoopAAs aas = loopdesign.allowed_aas( SurroundingSS( 'H', 'H' ), "AGBA" );
		for ( LoopAAs::const_iterator r = aas.begin(), endr = aas.end(); r != endr; ++r ) {
			TR << *r << " ";
		}
		TR << std::endl;
	}


};
