// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <core/scoring/constraints/AmbiguousConstraint.hh>
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
#include <protocols/moves/PyMolMover.hh>
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
static basic::Tracer TR("BetaNonlocal");

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
		score_fxn->set_weight( atom_pair_constraint, 1 );
		score_fxn->set_weight( dihedral_constraint, 0.3 );
		//score_fxn->set_weight( hbond_lr_bb, 50 );
		//score_fxn->set_weight( hbond_sr_bb, 0.2 );

		chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

		Pose pose;

		//pose.append_residue_by_jump( Residue( restype_set->name_map( "B3A:AcetylatedNtermProteinFull:MethylatedCtermProteinFull" ), true ), 1 );
		pose.append_residue_by_jump( Residue( restype_set->name_map( "B3A:AcetylatedNtermProteinFull" ), true ), 1 );
		pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		pose.append_residue_by_bond( Residue( restype_set->name_map( "NME" ), true ), true );

		pose.append_residue_by_jump( Residue( restype_set->name_map( "B3A:AcetylatedNtermProteinFull" ), true ), 1, "", "", true );
		pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A" ), true ), true );
		pose.append_residue_by_bond( Residue( restype_set->name_map( "B3A:MethylatedCtermProteinFull" ), true ), true );

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii )
			pose.set_torsion( TorsionID( ii, id::BB, 4 ), 180 );
		
		protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( pose, 1 ) );
		translate->step_size( 4 );
		translate->apply( pose );

		kinematics::MoveMapOP mm( new kinematics::MoveMap );
		mm->set_bb( false );
		mm->set_jump( true );
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii )
			mm->set( TorsionID( ii, id::BB, 4 ), false );
		
		protocols::simple_moves::MinMoverOP min( new protocols::simple_moves::MinMover( mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.001, true ) );
		
		min->apply( pose );
	
		utility::vector1< ConstraintOP > csts;
		
		utility::vector1< ConstraintOP > multi;
		for ( Size ii = 1; ii <= 5; ii += 2 ) {
			Size jj = 12 - ii;
			multi.push_back( ConstraintOP( new AtomPairConstraint(
				AtomID( pose.residue( ii ).atom_index( "O" ), ii ),
				AtomID( pose.residue( jj ).atom_index( "H" ), jj ),
				core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.1 ) )
			) ) );
			multi.push_back( ConstraintOP( new AtomPairConstraint(
				AtomID( pose.residue( ii ).atom_index( "H" ), ii ),
				AtomID( pose.residue( jj ).atom_index( "O" ), jj ),
				core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.1 ) )
			) ) );
		}
		csts.push_back( ConstraintOP( new MultiConstraint( multi ) ) );
		
		for ( Size ii = 1; ii <= 5; ii += 2 ) {
			Size jj = ii + 6;
			multi.push_back( ConstraintOP( new AtomPairConstraint(
				AtomID( pose.residue( ii ).atom_index( "O" ), ii ),
				AtomID( pose.residue( jj ).atom_index( "H" ), jj ),
				core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.1 ) )
			) ) );
			if ( ii == 5 ) continue;
			multi.push_back( ConstraintOP( new AtomPairConstraint(
				AtomID( pose.residue( ii ).atom_index( "H" ), ii ),
				AtomID( pose.residue( jj+1 ).atom_index( "O" ), jj ),
				core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.1 ) )
			) ) );
		}
		//pose.add_constraint( any_interchain_hbond );
		pose.add_constraint( ConstraintOP( new AmbiguousConstraint( csts ) ) );
		
		pose.dump_pdb( "out.pdb" );
		for ( Real phi = -170; phi <= 180; phi += 10 ) {
			for ( Real tht = -170; tht <= 180; tht += 10 ) {
				for ( Real psi = -170; psi <= 180; psi += 10 ) {
					for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
						if ( ii == 6 ) continue;
						pose.set_torsion( TorsionID( ii, id::BB, 1 ), phi );
						pose.set_torsion( TorsionID( ii, id::BB, 2 ), tht );
						pose.set_torsion( TorsionID( ii, id::BB, 3 ), psi );
					}
					
					min->apply( pose );
					std::stringstream fn;
					fn << "out_" << phi << "_" << tht << "_" << psi << ".pdb";
					if ( pose.energies().total_energy() < -5 ) {
						pose.dump_scored_pdb( fn.str(), *score_fxn );
					}
				}
			}
		}
		
		//AmbiguousConstraintOP any_interchain_hbond( new AmbiguousConstraint );
		/*
		 pose.add_constraint(
			ConstraintOP(
				new AtomPairConstraint(
					AtomID( pose.residue( 1 ).atom_index( "O" ), 1 ),
					AtomID( pose.residue( 3 ).atom_index( "H" ), 3 ),
					core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.04 ) )
		) ) );
		
		pose.add_constraint(
			ConstraintOP(
				new AtomPairConstraint(
					AtomID( pose.residue( 3 ).atom_index( "O" ), 3 ),
					AtomID( pose.residue( 1 ).atom_index( "H" ), 1 ),
					core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 1.9, 0.04 ) )
		) ) );
		
		
		pose.add_constraint(
			ConstraintOP(
				new DihedralConstraint(
					AtomID( pose.residue( 3 ).atom_index( "CA" ), 3 ),
					AtomID( pose.residue( 3 ).atom_index( "N" ), 3 ),
					AtomID( pose.residue( 3 ).atom_index( "CO" ), 3 ),
					AtomID( pose.residue( 3 ).atom_index( "CP2" ), 3 ),
					core::scoring::func::CircularHarmonicFuncOP( new core::scoring::func::CircularHarmonicFunc( 3.14159, 0.04 ) )
		) ) );
		pose.add_constraint(
			ConstraintOP(
				new DihedralConstraint(
					AtomID( pose.residue( 3 ).atom_index( "H" ), 3 ),
					AtomID( pose.residue( 3 ).atom_index( "N" ), 3 ),
					AtomID( pose.residue( 3 ).atom_index( "CO" ), 3 ),
					AtomID( pose.residue( 3 ).atom_index( "OP1" ), 3 ),
					core::scoring::func::CircularHarmonicFuncOP( new core::scoring::func::CircularHarmonicFunc( 3.14159, 0.04 ) )
		) ) );
		pose.add_constraint(
							ConstraintOP(
										 new DihedralConstraint(
																AtomID( pose.residue( 3 ).atom_index( "CA" ), 3 ),
																AtomID( pose.residue( 3 ).atom_index( "N" ), 3 ),
																AtomID( pose.residue( 3 ).atom_index( "CO" ), 3 ),
																AtomID( pose.residue( 3 ).atom_index( "CP2" ), 3 ),
																core::scoring::func::CircularHarmonicFuncOP( new core::scoring::func::CircularHarmonicFunc( 3.14159, 0.04 ) )
																) ) );
		*/
		
		
		

		
		// Long MC trajectory



	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main
