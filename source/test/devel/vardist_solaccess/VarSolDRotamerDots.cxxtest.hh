// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/vardist_solaccess/VarSolDRotamerDots.cxxtest.hh
/// @brief  test suite for Variable-distance hydrogen-bonds to solvent calculator
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Core headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueKinWriter.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

// Basic headers
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>


static basic::Tracer TR("VarSolDRotamerDots.cxxtest");

using namespace core;
using namespace core::chemical;
using namespace core::conformation;

class VarSolDRotamerDotsTest : public CxxTest::TestSuite
{
public:
	typedef core::pose::PoseOP PoseOP;
	typedef devel::vardist_solaccess::VarSolDRotamerDots VarSolDRotamerDots; 
	//typedef core::pack::task::PackerTaskOP PackerTaskOP;
	//typedef core::pack::task::TaskFactory TaskFactory;
	//typedef protocols::pack_daemon::PackDaemon PackDaemon;

public:
	void setUp() {
		core_init();
	}

	void test_self_overlap() {

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueOP res10 = new core::conformation::Residue( trpcage.residue( 10 ) );
		VarSolDRotamerDots dots( res10 );
		dots.increment_self_overlap();
		//ResidueKinWriter writer;
		//std::ofstream fout( "test.kin" );
		//writer.write_kin_header( fout, *res10 );
		//writer.write_rsd_coords( fout, *res10 );
		//dots.write_dot_kinemage( fout );

	}

	void test_residue_pair_overlap() {

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueOP res10 = new core::conformation::Residue( trpcage.residue( 10 ) );
		ResidueOP res11 = new core::conformation::Residue( trpcage.residue( 11 ) );
		VarSolDRotamerDots dots10( res10 );
		VarSolDRotamerDots dots11( res11 );
		dots10.increment_self_overlap();
		dots11.increment_self_overlap();
		dots10.intersect_residues( dots11 );

		TS_ASSERT( true );
		//ResidueKinWriter writer;
		//std::ofstream fout( "test2.kin" );
		//writer.write_kin_header( fout, *res10 );
		//writer.write_rsd_coords( fout, *res10 );
		//writer.write_rsd_coords( fout, *res11 );
		//dots10.write_dot_kinemage( fout );
		//dots11.write_dot_kinemage( fout );

	}


};


