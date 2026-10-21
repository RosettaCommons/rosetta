// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file test/protocols/recon_design/MSDMoverTest.cxxtest.hh
/// @brief Unit test for MSDMover
/// @author Alex Sevy (alex.sevy@gmail.com)

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <protocols/recon_design/MSDMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <util/pose_funcs.hh>

static basic::Tracer TR("protocols.recon_design.MSDMoverTest");

// --------------- Test Class --------------- //

using namespace protocols::recon_design;
using namespace protocols::simple_moves;
using namespace core::pose;
using namespace core::scoring::constraints;
using namespace protocols::minimization_packing;

class MSDMoverTest : public CxxTest::TestSuite {

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
		MSDMoverOP msd_mover ( new MSDMover );
		utility::vector1< std::string > resfiles;
		resfiles.push_back( "protocols/recon_design/design.resfile" );
		resfiles.push_back( "protocols/recon_design/design2.resfile" );
		PackRotamersMoverOP packer ( new PackRotamersMover );
		utility::vector1< std::string > single_resfile;
		single_resfile.push_back( "protocols/recon_design/design.resfile" );

		// test resfile parsing with multiple resfiles for msd mover
		msd_mover->set_poses( poses );
		msd_mover->resfiles( resfiles );
		msd_mover->design_mover( packer );
		msd_mover->setup_mover( *pose1 );

		utility::vector1< utility::vector1< core::Size > > res_links = msd_mover->designable_residues();

		TS_ASSERT_EQUALS( res_links.size(), 2 );
		TS_ASSERT_EQUALS( res_links[ 1 ].size(), 1 );
		TS_ASSERT_EQUALS( res_links[ 1 ][ 1 ], 1 );
		TS_ASSERT_EQUALS( res_links[ 2 ].size(), 1 );
		TS_ASSERT_EQUALS( res_links[ 2 ][ 1 ], 2 );

		// test resfile parsing for single resfile
		msd_mover->resfiles( single_resfile );
		msd_mover->setup_mover( *pose1 );

		res_links = msd_mover->designable_residues();
		TS_ASSERT_EQUALS( res_links.size(), 2 );
		TS_ASSERT_EQUALS( res_links[ 1 ].size(), 1 );
		TS_ASSERT_EQUALS( res_links[ 1 ][ 1 ], 1 );
		TS_ASSERT_EQUALS( res_links[ 2 ].size(), 1 );
		TS_ASSERT_EQUALS( res_links[ 2 ][ 1 ], 1 );
	}

	void test_correct_constraints () {

		utility::vector1< PoseOP > poses;
		MutateResidueOP mutate;
		PoseOP pose1 = pose_.clone();
		mutate = MutateResidueOP ( new MutateResidue( 1, "ALA" ) );
		mutate->apply( *pose1 );
		setPoseExtraScore( *pose1, "msd_job_dist_index", 1 );
		poses.push_back( pose1 );

		PoseOP pose2 = pose_.clone();
		mutate = MutateResidueOP ( new MutateResidue( 1, "CYS" ) );
		mutate->apply( *pose2 );
		setPoseExtraScore( *pose2, "msd_job_dist_index", 2 );
		poses.push_back( pose2 );

		PoseOP pose3 = pose_.clone();
		mutate = MutateResidueOP ( new MutateResidue( 1, "ASP" ) );
		mutate->apply( *pose3 );
		setPoseExtraScore( *pose3, "msd_job_dist_index", 3 );
		poses.push_back( pose3 );

		MSDMoverOP msd_mover ( new MSDMover );
		msd_mover->set_poses( poses );
		PackRotamersMoverOP packer ( new PackRotamersMover );
		utility::vector1< std::string > resfiles;
		resfiles.push_back( "protocols/recon_design/design.resfile" );

		// set resfiles, design mover, current pose, all that shit
		msd_mover->resfiles( resfiles );
		msd_mover->design_mover( packer );

		msd_mover->setup_mover( *pose1 );

		// Find the sequence of the other poses
		utility::vector1< utility::vector1< std::string > > other_pose_sequences;
		//other_pose_sequences.push_back( "A" );
		//other_pose_sequences.push_back( "C" );
		//other_pose_sequences.push_back( "D" );
		utility::vector1< std::string > ala_seq;
		ala_seq.push_back("ALA");
		other_pose_sequences.push_back(ala_seq);
		utility::vector1< std::string > cys_seq;
		cys_seq.push_back("CYS");
		other_pose_sequences.push_back(cys_seq);
		utility::vector1< std::string > asp_seq;
		asp_seq.push_back("ASP");
		other_pose_sequences.push_back(asp_seq);

		utility::vector1< utility::vector1< core::Size > > designable_residues = msd_mover->designable_residues();

		utility::vector1< ConstraintCOP > constraints = msd_mover->apply_linked_constraints( *pose1, other_pose_sequences, designable_residues[ 1 ]);


		ResidueTypeConstraintOP res_type1 ( new ResidueTypeConstraint( *pose1, 1, "ALA", 0.5 ) );
		ResidueTypeConstraintOP res_type2 ( new ResidueTypeConstraint( *pose1, 1, "CYS", 0.5 ) );
		ResidueTypeConstraintOP res_type3 ( new ResidueTypeConstraint( *pose1, 1, "ASP", 0.5 ) );
		for ( core::Size i = 1; i <= constraints.size(); ++i ) {
			std::cout << "constraint " << i << std::endl;
			constraints[i]->show( std::cout );
			std::cout << std::endl;
			TS_ASSERT(
				*constraints[ i ] == *res_type1 ||
				*constraints[ i ] == *res_type2 ||
				*constraints[ i ] == *res_type3
			);
		}

		msd_mover->setup_mover( *pose2 );
		constraints = msd_mover->apply_linked_constraints( *pose2, other_pose_sequences,  designable_residues[ 1 ] );
		res_type1 = ResidueTypeConstraintOP ( new ResidueTypeConstraint( *pose2, 1, "ALA", 0.5 ) );
		res_type2 = ResidueTypeConstraintOP ( new ResidueTypeConstraint( *pose2, 1, "CYS", 0.5 ) );
		res_type3 = ResidueTypeConstraintOP ( new ResidueTypeConstraint( *pose2, 1, "ASP", 0.5 ) );
		for ( core::Size i = 1; i <= constraints.size(); ++i ) {
			TS_ASSERT(
				*constraints[ i ] == *res_type1 ||
				*constraints[ i ] == *res_type2 ||
				*constraints[ i ] == *res_type3
			);
		}

	}

	void test_keeps_constraints () {

		utility::vector1< PoseOP > poses;
		PoseOP pose1 = pose_.clone();
		setPoseExtraScore( *pose1, "msd_job_dist_index", 1 );
		PoseOP pose2 = pose_.clone();
		setPoseExtraScore( *pose2, "msd_job_dist_index", 2 );
		poses.push_back( pose1 );
		poses.push_back( pose2 );

		MSDMoverOP msd_mover ( new MSDMover );
		msd_mover->set_poses( poses );
		utility::vector1< std::string > resfiles;
		resfiles.push_back( "protocols/recon_design/design.resfile" );
		PackRotamersMoverOP packer ( new PackRotamersMover );
		msd_mover->resfiles( resfiles );
		msd_mover->design_mover( packer );
		msd_mover->setup_mover( *poses[1] );
		ResidueTypeConstraintOP cst1 ( new ResidueTypeConstraint( *poses[1], 1, "CYS", 5.0 ) );
		poses[1]->add_constraint( cst1 );
		msd_mover->apply( *poses[1] );
		utility::vector1< ConstraintCOP > csts = poses[1]->constraint_set()->get_all_constraints();
		TS_ASSERT( *cst1 == *csts[1] );
		TS_ASSERT_EQUALS( csts.size(), 1 );
	}


private:
	core::pose::Pose pose_;
};
