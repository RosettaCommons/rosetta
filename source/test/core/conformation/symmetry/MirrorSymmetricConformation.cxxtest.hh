// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/Minimizer.cxxtest.hh
/// @brief  test suite for Minimizer
/// @author Phil Bradley
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <core/scoring/rms_util.hh>
#include <core/kinematics/Jump.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/Conformation.hh>
#include <src/core/id/AtomID.hh>
#include <src/core/id/TorsionID.hh>
#include <src/core/scoring/ScoreFunction.hh>
#include <src/core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <src/core/pack/pack_rotamers.hh>
#include <src/core/pack/make_symmetric_task.hh>
#include <src/core/pack/task/TaskFactory.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/scoring/aa_composition_energy/AACompositionConstraint.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.conformation.symmetry.SymmetricConformation.cxxtest");

using namespace core;

class MirrorSymmetricConformationTests : public CxxTest::TestSuite
{
	chemical::ResidueTypeSetCAP residue_set;

public:
	MirrorSymmetricConformationTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	///////////////////////////////////////////////////////////////////////////////
	// ------------------------------------------ //
	/// @brief simple test minimization
	void test_symmetric_conformation()
	{
		double rms_threshold = 1e-2;

		// 6 tests that should hit all the atomtree machinery
		// tests 1-5 should all have the same output (only test 6 perturbs the pose)

		pose::Pose pose1_ref, pose2_ref;
		import_pose::pose_from_file( pose1_ref, "core/conformation/symmetry/mtest1.pdb" , core::import_pose::PDB_file);
		import_pose::pose_from_file( pose2_ref, "core/conformation/symmetry/mtest6.pdb" , core::import_pose::PDB_file);

		// test 1: setup symmetric system
		pose::Pose start_pose;
		import_pose::pose_from_file( start_pose, "core/conformation/symmetry/Cs_INPUT.pdb" , core::import_pose::PDB_file);
		pose::Pose pose = start_pose;
		core::conformation::symmetry::SymmData symmdata1(  pose.size(),  pose.num_jump() );
		std::string symm_def1 = "core/conformation/symmetry/Cs.symm";
		symmdata1.read_symmetry_data_from_file(symm_def1);
		core::pose::symmetry::make_symmetric_pose( pose, symmdata1 );

		core::Real rms_to_restored = scoring::CA_rmsd( pose, pose1_ref );
		TR << "RMS difference after setup symmetric system: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < rms_threshold );

		// test 2: set_jump
		core::kinematics::Jump flexible_jump1 = pose.jump( 1 );
		pose.conformation().set_jump( 1, flexible_jump1 );
		rms_to_restored = scoring::CA_rmsd( pose, pose1_ref );
		TR << "RMS difference after set_jump: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < rms_threshold );


		// test 3: replace_residue (trigger jump atom replacement)
		core::conformation::Residue newres = pose.residue(1);
		pose.replace_residue(1, newres, false);
		rms_to_restored = scoring::CA_rmsd( pose, pose1_ref );
		TR << "RMS difference after replace_residue: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < rms_threshold );

		// test 4: set_xyz (trigger mirrored subunit refolding from inside subunit)
		pose.set_xyz(core::id::AtomID(4,1), numeric::xyzVector<core::Real> (-8.133 , -14.360 , 5.195) );
		rms_to_restored = scoring::CA_rmsd( pose, pose1_ref );
		TR << "RMS difference after set_xyz: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < rms_threshold );

		// test 5: set_jump (ensuring atom tree correct after test 3 & 4)
		core::kinematics::Jump flexible_jump2 = pose.jump( 1 );
		pose.conformation().set_jump( 1, flexible_jump2 );
		rms_to_restored = scoring::CA_rmsd( pose, pose1_ref );
		TR << "RMS difference after set_jump(2): " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < rms_threshold );

		// test 6: set_dof (trigger refold_from_rb_deltas)
		id::DOF_ID const & id_master1
			( pose.conformation().dof_id_from_torsion_id(id::TorsionID(1,id::JUMP,4)));
		pose.set_dof( id_master1, 60.0 );
		id::DOF_ID const & id_master2
			( pose.conformation().dof_id_from_torsion_id(id::TorsionID(1,id::JUMP,5)));
		pose.set_dof( id_master2, 60.0 );
		id::DOF_ID const & id_master3
			( pose.conformation().dof_id_from_torsion_id(id::TorsionID(1,id::JUMP,6)));
		pose.set_dof( id_master3, 60.0 );
		rms_to_restored = scoring::CA_rmsd( pose, pose2_ref );
		TR << "RMS difference after set_dof: " << rms_to_restored << std::endl;
		TS_ASSERT( rms_to_restored < rms_threshold );
	};

	/// @brief Tests design with the aa_composition score term.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void aacomp_maintest( bool const only_design_subset ) {
		core::pose::Pose pose;
		import_pose::pose_from_file( pose, "core/conformation/symmetry/Cs_INPUT.pdb" , core::import_pose::PDB_file);

		core::conformation::symmetry::SymmData symmdata1(  pose.size(),  pose.num_jump() );
		std::string symm_def1 = "core/conformation/symmetry/Cs.symm";
		symmdata1.read_symmetry_data_from_file(symm_def1);
		core::pose::symmetry::make_symmetric_pose( pose, symmdata1 );
		pose.update_residue_neighbors();
		//pose.dump_pdb("vmirrortemp1.pdb"); //DELETE ME

		//Little check: do we properly store whether glycines are mirrored?
		for ( core::Size ir=1, irmax=pose.total_residue(); ir<=irmax; ++ir ) {
			if ( pose.residue_type(ir).aa() == core::chemical::aa_gly ) {
				if ( pose.chain(ir) == 1 ) {
					TS_ASSERT( !pose.residue(ir).mirrored_relative_to_type() );
				} else {
					TS_ASSERT( pose.residue(ir).mirrored_relative_to_type() );
				}
			} else {
				TS_ASSERT( !pose.residue(ir).mirrored_relative_to_type() );
			}
		}

		// Set up score function
		core::scoring::symmetry::SymmetricScoreFunction scorefxn;
		scorefxn.set_weight( core::scoring::aa_composition, 1 );
		scorefxn.set_weight( core::scoring::fa_atr, 0.01 );
		scorefxn.set_weight( core::scoring::fa_rep, 0.01 );
		(scorefxn)(pose); //Score the pose once before calling the packer (necessary).

		// Set up the constraints in the pose
		core::scoring::aa_composition_energy::AACompositionConstraintOP comp_constraint( new core::scoring::aa_composition_energy::AACompositionConstraint );
		comp_constraint->initialize_from_file( "core/conformation/symmetry/mirror_symm.comp" ); //Requires exactly 5 DALA
		pose.add_constraint( comp_constraint );

		// Set up packer task and packer objects
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));
		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_ala ] = true;
		keep_aas[ core::chemical::aa_gly ] = true;
		for ( core::Size i = 1; i <= pose.size(); i++ ) {
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas, "PIKAA" );
			if ( only_design_subset && i < 20 ) {
				task->nonconst_residue_task( i ).prevent_repacking();
			}
		}
		task->request_symmetrize_by_intersection();
		task = core::pack::make_new_symmetric_PackerTask_by_requested_method( pose, task );


		core::pack::interaction_graph::AnnealableGraphBaseOP symmetric_ig;
		core::pack::rotamer_set::symmetry::SymmetricRotamerSetsOP sym_rotamer_sets( new core::pack::rotamer_set::symmetry::SymmetricRotamerSets );


		// Run the packer
		utility::vector0< int > rot_to_pack;
		core::pack::pack_rotamers_setup( pose, scorefxn, task, sym_rotamer_sets, symmetric_ig );
		for ( core::Size i=1; i<=5; ++i ) {
			core::pack::pack_rotamers_run( pose, task, sym_rotamer_sets, symmetric_ig, rot_to_pack );
			//Check that after packing, there are five ALA and five DALA residues in the pose.
			core::Size dalacount(0);
			core::Size alacount(0);
			TR << "Sequence:" << std::endl;
			for ( core::Size i=1, imax=pose.size(); i<=imax; ++i ) {
				TR << pose.residue(i).name3();
				if ( i<imax ) {
					if ( pose.residue(i).name3() == "DAL" || pose.residue(i).name3() == "ALA" ) TR << "* ";
					else TR << "  ";
				}
				if ( i % 5 == 0 ) TR << std::endl;
				if ( pose.residue(i).name3() == "DAL" ) ++dalacount;
				if ( pose.residue(i).name3() == "ALA" ) ++alacount;
			}
			TR << std::endl;
			TR.flush();
			TS_ASSERT_EQUALS( dalacount, 5);
			TS_ASSERT_EQUALS( alacount, 5);
		}
		//pose.dump_pdb("vmirrortemp2.pdb"); //DELETE ME
	}

	/// @brief Tests design with the aa_composition score term, with the whole pose being designed.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_aacomp() {
		TR << "Testing aa_composition design with mirror symmetry..." << std::endl;
		aacomp_maintest( false );
		TR.flush();
	}

	/// @brief Tests design with the aa_composition score term, with a subset of the pose being designed.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_aacomp_design_subset() {
		TR << "Testing aa_composition design with mirror symmetry, with only a subset of the pose being designed..." << std::endl;
		aacomp_maintest( true );
		TR.flush();
	}

};
