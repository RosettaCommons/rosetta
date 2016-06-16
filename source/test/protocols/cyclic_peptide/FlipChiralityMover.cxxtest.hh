// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/cyclic_peptide/FlipChiralityMover.cxxtest.hh
/// @brief  flips chirality of selected residues in a pmtion_se
/// @author Parisa Hosseinzadeh (parisah@uw.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/cyclic_peptide/FlipChiralityMover.hh>

// Core Headers
#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/util.hh>
#include <core/chemical/AA.hh>
// Protocol Headers
#include <basic/Tracer.hh>
#include <iostream>

static THREAD_LOCAL basic::Tracer TR("FlipChiralityMover");


class FlipChiralityMover : public CxxTest::TestSuite {
	//Define Variables

public:

	core::pose::PoseOP test_pose_;
	core::pose::PoseOP init_pose_;
	protocols::cyclic_peptide::FlipChiralityMoverOP flipper_;

	void setUp(){
		core_init();

		init_pose_ = core::pose::PoseOP( new core::pose::Pose() );
		test_pose_ = core::pose::PoseOP( new core::pose::Pose() );
		core::import_pose::pose_from_file( *init_pose_, "protocols/cyclic_peptide/dVdPN_A.pdb" , core::import_pose::PDB_file);
		core::import_pose::pose_from_file( *test_pose_, "protocols/cyclic_peptide/dVdPN_A.pdb" , core::import_pose::PDB_file);
		flipper_=protocols::cyclic_peptide::FlipChiralityMoverOP (new protocols::cyclic_peptide::FlipChiralityMover ());

		flipper_->apply(*test_pose_);


	}

	void tearDown(){

	}

	void test_size_range () {

		TS_ASSERT_EQUALS((*test_pose_).n_residue(),3);
	}

	void test_name_change (){

		TS_ASSERT((*test_pose_).residue(1).name3()=="VAL" && (*test_pose_).residue(2).name3()=="PRO" && (*test_pose_).residue(3).name3()=="DAN");
	}

	void test_mirrored () {


		core::Size x_value=0;
		core::Size y_value=0;
		core::Size z_value=0;


		for ( core::Size i=1; i<= init_pose_->n_residue(); ++i ) {

			core::conformation::Residue const & rsd_init((*init_pose_).residue(i) );
			core::conformation::Residue const & rsd_test((*test_pose_).residue(i) );
			numeric::xyzVector<core::Real> res1_coord_init=rsd_init.atom(1).xyz();
			numeric::xyzVector<core::Real> res1_coord_test=rsd_test.atom(1).xyz();
			core::Size z_value_init=(res1_coord_init.z()+res1_coord_test.z());

			for ( core::Size j=1; j<= rsd_init.nheavyatoms(); ++j ) {
				core::conformation::Atom const & atom_init( rsd_init.atom(j) );
				core::conformation::Atom const & atom_test( rsd_test.atom(j) );
				numeric::xyzVector< core::Real > init_coord( 0.0 );
				numeric::xyzVector< core::Real > test_coord( 0.0 );
				init_coord=atom_init.xyz();
				test_coord=atom_test.xyz();


				x_value+=(init_coord.x()-test_coord.x());
				y_value+=(init_coord.y()-test_coord.y());
				core::Size z_value_sec=init_coord.z()+test_coord.z();
				if ( z_value_init == z_value_sec ) {
					z_value+=0;
				} else {
					z_value+=1;
				}
			}
		}

		TS_ASSERT (x_value==0 && y_value==0 &&z_value==0);
	}

};



