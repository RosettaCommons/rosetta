// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/relax/FastRelax.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/RandomTorsionMover.hh>
#include <protocols/simple_moves/hbs/HbsPatcher.hh>
#include <protocols/simple_moves/a3b_hbs/A3BHbsPatcher.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>


//Basic headers
#include <basic/resource_manager/ResourceManager.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.functions.hh>

// C++ headers
#include <string>
#include <sstream>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::a3b_hbs;
using namespace protocols::simple_moves::chiral;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer - used to replace cout
static THREAD_LOCAL basic::Tracer TR("BetaNonlocal");

int
main( int argc, char* argv[] )
{
	try {

		// init command line options
		//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
		devel::init(argc, argv);

		using namespace core::scoring::constraints;

		// create score function
		scoring::ScoreFunctionOP score_fxn( scoring::get_score_function() );
		score_fxn->set_weight( fa_rep, 0.3 );
		score_fxn->set_weight( atom_pair_constraint, 0.3 );
		score_fxn->set_weight( dihedral_constraint, 0.3 );

		chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

		Pose a3bpose;
		Pose b314pose;

		a3bpose.append_residue_by_jump( Residue( restype_set->name_map( "ALA:AcetylatedNtermProteinFull" ), true ), 1 );
		a3bpose.append_residue_by_bond( Residue( restype_set->name_map( "B3L" ), true ), true );
		a3bpose.append_residue_by_bond( Residue( restype_set->name_map( "ALA" ), true ), true );
		a3bpose.append_residue_by_bond( Residue( restype_set->name_map( "ALA" ), true ), true );
		a3bpose.append_residue_by_bond( Residue( restype_set->name_map( "ALA" ), true ), true );
		a3bpose.append_residue_by_bond( Residue( restype_set->name_map( "B3L" ), true ), true );
		a3bpose.append_residue_by_bond( Residue( restype_set->name_map( "ALA" ), true ), true );
		a3bpose.append_residue_by_bond( Residue( restype_set->name_map( "ALA" ), true ), true );
		a3bpose.append_residue_by_bond( Residue( restype_set->name_map( "ALA" ), true ), true );
		a3bpose.append_residue_by_bond( Residue( restype_set->name_map( "B3L" ), true ), true );
		a3bpose.append_residue_by_bond( Residue( restype_set->name_map( "ALA:MethylatedCtermProteinFull" ), true ), true );

		for ( Size ii = 1; ii <= a3bpose.size(); ++ii ) {
			if ( a3bpose.residue_type( ii ).is_beta_aa() ) {
				a3bpose.conformation().set_torsion( TorsionID( ii, id::BB, 1), -115);
				a3bpose.conformation().set_torsion( TorsionID( ii, id::BB, 2),   85);
				a3bpose.conformation().set_torsion( TorsionID( ii, id::BB, 3),  -95);
				a3bpose.conformation().set_torsion( TorsionID( ii, id::BB, 4),  180);
			} else {
				a3bpose.conformation().set_torsion( TorsionID( ii, id::BB, 1),  -65);
				a3bpose.conformation().set_torsion( TorsionID( ii, id::BB, 2),  -45);
				a3bpose.conformation().set_torsion( TorsionID( ii, id::BB, 3),  180);
			}
			if ( a3bpose.residue_type( ii ).name() == "B3L" ) {
				a3bpose.conformation().set_torsion( TorsionID( ii, id::CHI, 1),  180);
				a3bpose.conformation().set_torsion( TorsionID( ii, id::CHI, 2),   60);
			}
		}

		kinematics::MoveMapOP mm( new kinematics::MoveMap );
		mm->set_chi( false );
		mm->set_bb( true );
		protocols::simple_moves::MinMoverOP min( new protocols::simple_moves::MinMover( mm, score_fxn, "linmin_iterated", 0.01, true ) );
		min->apply( a3bpose );

		b314pose.append_residue_by_jump( Residue( restype_set->name_map( "B3A:AcetylatedNtermProteinFull" ), true ), 1 );
		b314pose.append_residue_by_bond( Residue( restype_set->name_map( "B3L" ), true ), true );
		b314pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		b314pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		b314pose.append_residue_by_bond( Residue( restype_set->name_map( "B3L" ), true ), true );
		b314pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		b314pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		b314pose.append_residue_by_bond( Residue( restype_set->name_map( "B3L" ), true ), true );
		b314pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		b314pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		b314pose.append_residue_by_bond( Residue( restype_set->name_map( "B3L" ), true ), true );
		b314pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A:MethylatedCtermProteinFull" ), true ), true );
		//-139.9 59.5 -138.7 180.0

		for ( Size ii = 1; ii <= b314pose.size(); ++ii ) {
			b314pose.conformation().set_torsion( TorsionID( ii, id::BB, 1), -139.9);
			b314pose.conformation().set_torsion( TorsionID( ii, id::BB, 2),   59.5);
			b314pose.conformation().set_torsion( TorsionID( ii, id::BB, 3), -138.7);
			b314pose.conformation().set_torsion( TorsionID( ii, id::BB, 4),  180);
			if ( b314pose.residue_type( ii ).name() == "B3L" ) {
				b314pose.conformation().set_torsion( TorsionID( ii, id::CHI, 1),  180);
				b314pose.conformation().set_torsion( TorsionID( ii, id::CHI, 2),   60);
			}
		}

		a3bpose.dump_pdb( "a3b_a_gp.pdb" );
		b314pose.dump_pdb( "b314_a_gp.pdb" );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main
