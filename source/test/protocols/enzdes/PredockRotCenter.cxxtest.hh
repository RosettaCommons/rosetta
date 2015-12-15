// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/DihedralConstraint.cxxtest.hh
/// @brief  test suite for constraints between protein and ligand
/// @author Florian Richter

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
//#include <core/conformation/Residue.hh>


//#include <core/kinematics/MoveMap.hh>

//#include <core/optimization/AtomTreeMinimizer.hh>
//#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh> //need for additional residue
#include <core/chemical/ResidueTypeSet.hh>
#include <basic/options/option.hh> //needed to set option
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //function for reading cstfiles
#include <protocols/enzdes/EnzdesMovers.hh> // for testing rot_center
#include <protocols/rigid/RigidBodyMover.hh>

//minimization stuff
#include <protocols/moves/Mover.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.LigInterfaceConstraints.cxxtest");

using namespace core;

class PredockRotCenterTest : public CxxTest::TestSuite
{

public:
	PredockRotCenterTest() {};
	protocols::toolbox::match_enzdes_util::EnzConstraintIOOP enz_io;
	protocols::enzdes::PredesignPerturbMoverOP predock;

	// Shared initialization goes here.
	void setUp() {
		core_init();
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCOP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if ( !residue_set.has_name("D2N") ) params_files.push_back("protocols/enzdes/D2N.params");
		residue_set.read_files_for_custom_residue_types(params_files);
		basic::options::option[basic::options::OptionKeys::run::preserve_header ].value(true);

		enz_io = protocols::toolbox::match_enzdes_util::EnzConstraintIOOP( new protocols::toolbox::match_enzdes_util::EnzConstraintIO(const_residue_set) );

		predock = protocols::enzdes::PredesignPerturbMoverOP( new protocols::enzdes::PredesignPerturbMover() );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_predock_rotcenter()
	{
		pose::Pose pose;
		core::import_pose::pose_from_pdb(pose, "protocols/enzdes/ligtest_it.pdb");
		scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );
		scorefxn->reset();
		enz_io->read_enzyme_cstfile("protocols/enzdes/ligtest_it.cst");
		enz_io->add_constraints_to_pose(pose, scorefxn, false);

		predock->set_ligand(107);
		predock->find_constraints_to_ligand(pose);

		core::Vector rot_center;

		rot_center = (predock->find_geometric_center_for_constrained_lig_atoms(pose));
		TS_ASSERT_DELTA(-4.987,rot_center[0],0.001);
		TS_ASSERT_DELTA(-4.415,rot_center[1],0.001);
		TS_ASSERT_DELTA(2.9765,rot_center[2],0.0001);
	}

};
