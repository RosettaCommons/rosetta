// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file test/protocols/recon_design/FindConsensusSequenceTest.cxxtest.hh
/// @brief Unit test for FindConsensusSequence
/// @author Alex Sevy (alex.sevy@gmail.com)

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/recon_design/FindConsensusSequence.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <basic/Tracer.hh>
#include <util/pose_funcs.hh>

static basic::Tracer TR("protocols.recon_design.FindConsensusSequenceTest");

using namespace protocols::recon_design;
using namespace core::pose;
using namespace core::scoring::constraints;
using namespace protocols::minimization_packing;
using namespace protocols::simple_moves;

class FindConsensusSequenceTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		pose_ = create_trpcage_ideal_pose();
	}

	void tearDown() {
	}

	void test_resfile_parsing () {

		utility::vector1< PoseOP > poses;
		MutateResidueOP mutate;
		PoseOP pose1 = pose_.clone();
		PoseOP pose2 = pose_.clone();
		setPoseExtraScore( *pose1, "msd_job_dist_index", 1 );
		setPoseExtraScore( *pose2, "msd_job_dist_index", 2 );
		poses.push_back( pose1 );
		poses.push_back( pose2 );

		utility::vector1< std::string > resfiles;
		resfiles.push_back( "protocols/recon_design/design.resfile" );
		resfiles.push_back( "protocols/recon_design/design2.resfile" );
		utility::vector1< std::string > single_resfile;
		single_resfile.push_back( "protocols/recon_design/design.resfile" );

		FindConsensusSequenceOP fcs ( new FindConsensusSequence );
		fcs->set_poses( poses );

		// test resfile parsing when a single resfile is given
		fcs->resfiles( single_resfile );
		fcs->parse_resfiles();
		utility::vector1< utility::vector1< core::Size > > res_links = fcs->designable_residues();
		TS_ASSERT_EQUALS( res_links.size(), 2 );
		TS_ASSERT_EQUALS( res_links[ 1 ].size(), 1 );
		TS_ASSERT_EQUALS( res_links[ 1 ][ 1 ], 1 );
		TS_ASSERT_EQUALS( res_links[ 2 ].size(), 1 );
		TS_ASSERT_EQUALS( res_links[ 2 ][ 1 ], 1 );

		// test resfile parsing when multiple resfiles are given
		fcs->resfiles( resfiles );
		fcs->parse_resfiles();
		res_links = fcs->designable_residues();
		TS_ASSERT_EQUALS( res_links.size(), 2 );
		TS_ASSERT_EQUALS( res_links[ 1 ].size(), 1 );
		TS_ASSERT_EQUALS( res_links[ 1 ][ 1 ], 1 );
		TS_ASSERT_EQUALS( res_links[ 2 ].size(), 1 );
		TS_ASSERT_EQUALS( res_links[ 2 ][ 1 ], 2 );
	}

	void test_consensus_sequence () {
		FindConsensusSequenceOP fcs ( new FindConsensusSequence );
		utility::vector1< PoseOP > poses;
		PackRotamersMoverOP packer ( new PackRotamersMover );
		core::pack::task::TaskFactoryOP task_factory ( new core::pack::task::TaskFactory );
		core::pack::task::operation::RestrictToRepackingOP rtr ( new core::pack::task::operation::RestrictToRepacking );
		task_factory->push_back( rtr );
		packer->task_factory( task_factory );
		MutateResidueOP mutate;
		PoseOP pose1 = pose_.clone();
		mutate = MutateResidueOP ( new MutateResidue( 1, "ASN" ) );
		mutate->apply( *pose1 );
		setPoseExtraScore( *pose1, "msd_job_dist_index", 1 );

		poses.push_back( pose1 );
		packer->apply( *pose1 );
		// TR << "Asparagine: " << pose1->energies().total_energy() << std::endl;
		PoseOP pose2 = pose_.clone();
		mutate = MutateResidueOP ( new MutateResidue( 1, "ARG" ) );
		mutate->apply( *pose2 );
		setPoseExtraScore( *pose2, "msd_job_dist_index", 2 );
		poses.push_back( pose2 );
		packer->apply( *pose2 );
		// TR << "Arginine: " << pose2->energies().total_energy() << std::endl;
		PoseOP pose3 = pose_.clone();
		mutate = MutateResidueOP ( new MutateResidue( 1, "GLY" ) );
		mutate->apply( *pose3 );
		setPoseExtraScore( *pose3, "msd_job_dist_index", 3 );
		poses.push_back( pose3 );
		packer->apply( *pose3 );
		// TR << "Glycine: " << pose3->energies().total_energy() << std::endl;

		utility::vector1< PoseOP > poses2;
		poses2.push_back( pose2->clone() );
		poses2.push_back( pose3->clone() );
		poses2.push_back( pose1->clone() );

		fcs->set_poses( poses );
		utility::vector1< std::string > resfiles;
		resfiles.push_back( "protocols/recon_design/design.resfile" );
		fcs->resfiles( resfiles );
		fcs->apply( *poses[1] );
		// Make sure all poses end up with the same consensus sequence
		for ( core::Size i = 2; i <= poses.size(); ++i ) {
			TS_ASSERT_EQUALS( poses[ i ]->sequence(), poses[ 1 ]->sequence() );
		}

		// Make sure you see the same consensus sequence when the order is changed
		fcs->set_poses( poses2 );
		fcs->apply( *poses2[ 1 ] );
		for ( core::Size i = 2; i <= poses2.size(); ++i ) {
			TS_ASSERT_EQUALS( poses2[ i ]->sequence(), poses2[ 1 ]->sequence() );
			TS_ASSERT_EQUALS( poses2[ i ]->sequence(), poses[ 1 ]->sequence() );
		}

		// Since asparagine is lowest energy residue here, make sure it's picked every time
		// I think this should be fine since the energy difference btwn N and G or R is large
		// enough that it should be reproducible, but it's possible that the random nature of
		// rotamer packing could cause this to fail
		// hpark: in REF2015 the lowest energy AA has been changed to GLY (from ASN)
		TS_ASSERT_EQUALS( poses[ 1 ]->residue( 1 ).name3(), "GLY"/*"ASN"*/ );
		TS_ASSERT_EQUALS( poses2[ 1 ]->residue( 1 ).name3(), "GLY"/*"ASN"*/ );
	}

private:
	core::pose::Pose pose_;
};
