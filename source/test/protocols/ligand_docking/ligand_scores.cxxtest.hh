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
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/util.hh>

#include <protocols/ligand_docking/ligand_scores.hh>

#include <utility/vector1.hh>

static basic::Tracer TR("protocols.ligand_docking.ligand_scores.cxxtest");


std::string const phe1(
	"ATOM    214  N   PHE X   1       2.902  -4.187  20.617  1.00  1.00           N  \n"
	"ATOM    215  CA  PHE X   1       3.318  -2.920  21.197  1.00  1.00           C  \n"
	"ATOM    216  C   PHE X   1       2.488  -1.820  20.560  1.00  1.00           C  \n"
	"ATOM    217  O   PHE X   1       2.326  -1.787  19.346  1.00  1.00           O  \n"
	"ATOM    218  CB  PHE X   1       4.809  -2.665  20.969  1.00  1.00           C  \n"
	"ATOM    219  CG  PHE X   1       5.706  -3.632  21.689  1.00  1.00           C  \n"
	"ATOM    220  CD1 PHE X   1       6.164  -4.778  21.055  1.00  1.00           C  \n"
	"ATOM    221  CD2 PHE X   1       6.092  -3.399  23.000  1.00  1.00           C  \n"
	"ATOM    222  CE1 PHE X   1       6.989  -5.668  21.716  1.00  1.00           C  \n"
	"ATOM    223  CE2 PHE X   1       6.918  -4.287  23.662  1.00  1.00           C  \n"
	"ATOM    224  CZ  PHE X   1       7.366  -5.423  23.019  1.00  1.00           C  \n"
);

// 10 Ang from phe1
std::string const phe2(
	"ATOM    214  N   PHE X   2      12.902  -4.187  20.617  1.00  1.00           N  \n"
	"ATOM    215  CA  PHE X   2      13.318  -2.920  21.197  1.00  1.00           C  \n"
	"ATOM    216  C   PHE X   2      12.488  -1.820  20.560  1.00  1.00           C  \n"
	"ATOM    217  O   PHE X   2      12.326  -1.787  19.346  1.00  1.00           O  \n"
	"ATOM    218  CB  PHE X   2      14.809  -2.665  20.969  1.00  1.00           C  \n"
	"ATOM    219  CG  PHE X   2      15.706  -3.632  21.689  1.00  1.00           C  \n"
	"ATOM    220  CD1 PHE X   2      16.164  -4.778  21.055  1.00  1.00           C  \n"
	"ATOM    221  CD2 PHE X   2      16.092  -3.399  23.000  1.00  1.00           C  \n"
	"ATOM    222  CE1 PHE X   2      16.989  -5.668  21.716  1.00  1.00           C  \n"
	"ATOM    223  CE2 PHE X   2      16.918  -4.287  23.662  1.00  1.00           C  \n"
	"ATOM    224  CZ  PHE X   2      17.366  -5.423  23.019  1.00  1.00           C  \n"
);

// 10 Ang from phe2, 20 from phe1
std::string const phe3(
	"ATOM    214  N   PHE X   3      22.902  -4.187  20.617  1.00  1.00           N  \n"
	"ATOM    215  CA  PHE X   3      23.318  -2.920  21.197  1.00  1.00           C  \n"
	"ATOM    216  C   PHE X   3      22.488  -1.820  20.560  1.00  1.00           C  \n"
	"ATOM    217  O   PHE X   3      22.326  -1.787  19.346  1.00  1.00           O  \n"
	"ATOM    218  CB  PHE X   3      24.809  -2.665  20.969  1.00  1.00           C  \n"
	"ATOM    219  CG  PHE X   3      25.706  -3.632  21.689  1.00  1.00           C  \n"
	"ATOM    220  CD1 PHE X   3      26.164  -4.778  21.055  1.00  1.00           C  \n"
	"ATOM    221  CD2 PHE X   3      26.092  -3.399  23.000  1.00  1.00           C  \n"
	"ATOM    222  CE1 PHE X   3      26.989  -5.668  21.716  1.00  1.00           C  \n"
	"ATOM    223  CE2 PHE X   3      26.918  -4.287  23.662  1.00  1.00           C  \n"
	"ATOM    224  CZ  PHE X   3      27.366  -5.423  23.019  1.00  1.00           C  \n"
);

// Flipped ring from phe1, 0.1 in Z from phe1
std::string const phe4(
	"ATOM    214  N   PHE X   4       2.902  -4.187  20.717  1.00  1.00           N  \n"
	"ATOM    215  CA  PHE X   4       3.318  -2.920  21.297  1.00  1.00           C  \n"
	"ATOM    216  C   PHE X   4       2.488  -1.820  20.660  1.00  1.00           C  \n"
	"ATOM    217  O   PHE X   4       2.326  -1.787  19.446  1.00  1.00           O  \n"
	"ATOM    218  CB  PHE X   4       4.809  -2.665  21.069  1.00  1.00           C  \n"
	"ATOM    219  CG  PHE X   4       5.706  -3.632  21.789  1.00  1.00           C  \n"
	"ATOM    220  CD2 PHE X   4       6.164  -4.778  21.155  1.00  1.00           C  \n"
	"ATOM    221  CD1 PHE X   4       6.092  -3.399  23.100  1.00  1.00           C  \n"
	"ATOM    222  CE2 PHE X   4       6.989  -5.668  21.816  1.00  1.00           C  \n"
	"ATOM    223  CE1 PHE X   4       6.918  -4.287  23.762  1.00  1.00           C  \n"
	"ATOM    224  CZ  PHE X   4       7.366  -5.423  23.119  1.00  1.00           C  \n"
);

class LigandScoresUtilTests : public CxxTest::TestSuite {

	core::pose::PoseOP pose1_, pose2_, pose3_, pose4_;
	core::Real delta_ = 0.0001;

public:

	void setUp() {
		core_init_with_additional_options("");

		pose1_ = fullatom_poseop_from_string(phe1);
		pose2_ = fullatom_poseop_from_string(phe2);
		pose3_ = fullatom_poseop_from_string(phe3);
		pose4_ = fullatom_poseop_from_string(phe4);
	}

	void tearDown() {}

	void test_rmsds() {

		std::map< std::string, core::Real > rmsds;

		// Standard method.

		rmsds = protocols::ligand_docking::get_ligand_RMSDs('X', *pose1_, *pose2_ );
		TS_ASSERT_DELTA( rmsds["ligand_rms_no_super_X"], 10.0 , delta_);
		rmsds = protocols::ligand_docking::get_ligand_RMSDs('X', *pose1_, *pose3_ );
		TS_ASSERT_DELTA( rmsds["ligand_rms_no_super_X"], 20.0 , delta_);
		rmsds = protocols::ligand_docking::get_ligand_RMSDs('X', *pose1_, *pose4_ );
		TS_ASSERT_DELTA( rmsds["ligand_rms_no_super_X"], 0.1 , delta_);

		// Ensemble methods
		core::pose::Pose refpose( *pose3_ );

		rmsds = protocols::ligand_docking::get_ligand_RMSDs('X', *pose1_, refpose, "", true );
		TS_ASSERT_DELTA( rmsds["ligand_rms_no_super_X"], 20.0 , delta_);

		refpose.append_pose_by_jump( *pose2_, 1 );
		rmsds = protocols::ligand_docking::get_ligand_RMSDs('X', *pose1_, refpose, "", true );
		TS_ASSERT_DELTA( rmsds["ligand_rms_no_super_X"], 10.0 , delta_);

		refpose.append_pose_by_jump( *pose4_, 1 );
		rmsds = protocols::ligand_docking::get_ligand_RMSDs('X', *pose1_, refpose, "", true );
		TS_ASSERT_DELTA( rmsds["ligand_rms_no_super_X"], 0.1 , delta_);

		refpose.append_pose_by_jump( *pose1_, 1 );
		rmsds = protocols::ligand_docking::get_ligand_RMSDs('X', *pose1_, refpose, "", true );
		TS_ASSERT_DELTA( rmsds["ligand_rms_no_super_X"], 0.0 , delta_);

		// We still error out if we have mis-matched entries and aren't using the ensemble method
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING( protocols::ligand_docking::get_ligand_RMSDs('X', *pose1_, refpose, "", false ) );
	}

	void test_ligand_travel() {
		std::map< std::string, core::Real > travel;

		// Standard method.
		travel = protocols::ligand_docking::get_ligand_travel('X', *pose1_, *pose2_ );
		TS_ASSERT_DELTA( travel["ligand_centroid_travel_X"], 10.0 , delta_);
		travel = protocols::ligand_docking::get_ligand_travel('X', *pose1_, *pose3_ );
		TS_ASSERT_DELTA( travel["ligand_centroid_travel_X"], 20.0 , delta_);
		travel = protocols::ligand_docking::get_ligand_travel('X', *pose1_, *pose4_ );
		TS_ASSERT_DELTA( travel["ligand_centroid_travel_X"], 0.1 , delta_);

		// Ensemble methods
		core::pose::Pose refpose( *pose3_ );

		travel = protocols::ligand_docking::get_ligand_travel('X', *pose1_, refpose, "", true );
		TS_ASSERT_DELTA( travel["ligand_centroid_travel_X"], 20.0 , delta_);

		refpose.append_pose_by_jump( *pose2_, 1 );
		travel = protocols::ligand_docking::get_ligand_travel('X', *pose1_, refpose, "", true );
		TS_ASSERT_DELTA( travel["ligand_centroid_travel_X"], 10.0 , delta_);

		refpose.append_pose_by_jump( *pose4_, 1 );
		travel = protocols::ligand_docking::get_ligand_travel('X', *pose1_, refpose, "", true );
		TS_ASSERT_DELTA( travel["ligand_centroid_travel_X"], 0.1 , delta_);

		refpose.append_pose_by_jump( *pose1_, 1 );
		travel = protocols::ligand_docking::get_ligand_travel('X', *pose1_, refpose, "", true );
		TS_ASSERT_DELTA( travel["ligand_centroid_travel_X"], 0.0 , delta_);
	}


};

