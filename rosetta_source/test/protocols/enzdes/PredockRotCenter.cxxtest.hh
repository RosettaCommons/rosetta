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

#include <core/io/pdb/pose_io.hh>

//#include <core/kinematics/MoveMap.hh>

//#include <core/optimization/AtomTreeMinimizer.hh>
//#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
//#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/chemical/ChemicalManager.hh> //need for additional residue
#include <core/chemical/ResidueTypeSet.hh>
#include <basic/options/option.hh> //needed to set option
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/BoundConstraint.hh> //need function in this file
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //function for reading cstfiles
#include <protocols/enzdes/EnzdesMovers.hh> // for testing rot_center
#include <protocols/moves/RigidBodyMover.hh>

//minimization stuff
#include <core/kinematics/MoveMap.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <math.h>  //need for sqrt taking

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>

//Auto Headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/io/pdb/file_data.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/scoring/types.hh>
#include <utility/stream_util.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
#include <utility/keys/Key2Tuple.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/format.hh>


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
    ResidueTypeSetCAP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
    ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
    if(!residue_set.has_name("D2N")) params_files.push_back("protocols/enzdes/D2N.params");
       residue_set.read_files(params_files,
       ChemicalManager::get_instance()->atom_type_set( FA_STANDARD ),
       ChemicalManager::get_instance()->element_set( FA_STANDARD ),
       ChemicalManager::get_instance()->mm_atom_type_set( FA_STANDARD ),
       ChemicalManager::get_instance()->orbital_type_set(FA_STANDARD));//,
       //ChemicalManager::get_instance()->csd_atom_type_set( FA_STANDARD ));
       basic::options::option[basic::options::OptionKeys::run::preserve_header ].value(true);

       enz_io = new protocols::toolbox::match_enzdes_util::EnzConstraintIO(& residue_set);

       predock = new protocols::enzdes::PredesignPerturbMover();
  }

  // Shared finalization goes here.
  void tearDown() {
  }

  void test_predock_rotcenter()
  {
    pose::Pose pose;
    core::import_pose::pose_from_pdb(pose, "protocols/enzdes/ligtest_it.pdb");
    scoring::ScoreFunctionOP scorefxn = new scoring::ScoreFunction;
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
