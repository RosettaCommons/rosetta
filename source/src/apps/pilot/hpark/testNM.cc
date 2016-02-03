// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <core/id/AtomID.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/normalmode/NormalMode.hh>
#include <protocols/normalmode/NormalModeRelaxMover.hh>
#include <protocols/normalmode/NormalModeMinimizer.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>

#include <core/types.hh>
#include <devel/init.hh>
#include <sys/time.h>

OPT_1GRP_KEY(IntegerVector, fpd, rsd_changing)
OPT_1GRP_KEY(IntegerVector, fpd, modes)

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

void test_NM( pose::Pose pose ){
	//NM.torsion( true );
	//NM.eckart_correction( false );

	protocols::normalmode::NormalMode NM( "CA", 9.0 );
	
	NM.solve( pose );

	utility::vector1< Real > eigv;
	eigv = NM.get_eigvec_tor( 1 );

	std::cout << "NDOF: " << eigv.size() << std::endl;
	std::cout << "Mode1!" << std::endl;

	for( Size itor = 1; itor <= eigv.size(); ++itor ){
		if( std::abs(eigv[itor]) > 0.05 ) 
			std::cout << itor << " " << eigv[itor] << std::endl;
	}

}

void test_NMrelaxer( pose::Pose pose ){

	core::scoring::ScoreFunctionCOP sfxn 
    = core::scoring::ScoreFunctionFactory::create_score_function( option[score::weights] );

	// Movemap setup by reading fpd options
	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	// Changing residue setup
	if( option[ fpd::rsd_changing ].user() ){
		Size const nres = option[ fpd::rsd_changing ].size();
		for( int i = 1; i <= nres; ++i ){
			int resno( option[ fpd::rsd_changing ]()[i] );
			mm->set_bb( (Size)(resno), true );
		}
	} else {
		mm->set_bb( true );
	}


	// Modes setup
	utility::vector1< Size > modes;
	if( option[ fpd::modes ].user() ){
		Size const nmode = option[ fpd::modes ].size();
		for( int i = 1; i <= nmode; ++i ){
			int modeno( option[ fpd::modes ]()[i] );
			modes.push_back( (Size)(modeno) );
		}
	} else {
		modes.push_back( 1 );
	}

	// Run CartNM & TorsionMin

	// Way 1
	protocols::normalmode::TorsionNormalModeMover 
		NMMover( pose, sfxn, mm, "CA", 10.0, "extrapolate" );

	// Way 2
	/*
	protocols::normalmode::CartesianNormalModeMover
		NMMover( pose, sfxn, mm, "CA", 10.0, "min" );
	NMMover.set_cartesian_minimize( false );
	*/

	optimization::MinimizerOptionsOP minoption
		= new optimization::MinimizerOptions( "lbfgs_armijo_nonmonotone",
																					0.001, true, false, false );

	minoption->max_iter( 50 );
	NMMover.set_minoption( minoption );

	// Iter over different scales
	for(Size i_mode = 1; i_mode <= modes.size(); ++i_mode ){
		Size const modeno( modes[i_mode] );
		NMMover.set_mode( modeno );

		// Positive direction
		for( Size i = 1; i <= 25; ++i ){
			// Refresh normal mode for current pose
			pose::Pose pose_tmp( pose );

			NMMover.set_extrapolate_scale( i*0.1 );
			std::cout << "positive, " << i << std::endl;
			NMMover.apply( pose_tmp );

			std::stringstream pdbname("");
			pdbname << "relax.p" << modeno << "." << i << ".pdb";
			pose_tmp.dump_pdb( pdbname.str() );
		}

		// Negative direction
		NMMover.invert_direction();
		for( Size i = 1; i <= 25; ++i ){
			pose::Pose pose_tmp( pose );

			NMMover.set_extrapolate_scale( i*0.1 );
			std::cout << "negative, " << i << std::endl;
			NMMover.apply( pose_tmp );

			std::stringstream pdbname("");
			pdbname << "relax.n" << modeno << "." << i << ".pdb";
			pose_tmp.dump_pdb( pdbname.str() );
		}
	}
}

void test_NMmin( pose::Pose pose ){

	protocols::normalmode::NormalModeMinimizer NMmin;
	core::optimization::AtomTreeMinimizer ATmin;
	core::kinematics::MoveMap mm;
	mm.set_bb ( true ); // Set only backbone torsions
	mm.set_chi ( true );
	mm.set_jump( true );
	//mm.set( core::id::THETA, option[ OptionKeys::relax::minimize_mainchain_bond_angles ]() );
	//mm.set( core::id::D, option[ OptionKeys::relax::minimize_mainchain_bond_lengths ]() );

	core::scoring::ScoreFunctionCOP sfxn 
    = core::scoring::ScoreFunctionFactory::create_score_function( option[score::weights] );

	core::optimization::MinimizerOptions options( "dfpmin_armijo", 0.01, 
																								true, 
																								false, false );

	// Setup for density scoring
	protocols::electron_density::SetupForDensityScoringMover mover;
	mover.apply( pose );

	/*
	utility::vector1< Size > modes_using;
	modes_using.push_back( 1 );
	modes_using.push_back( 2 );
	modes_using.push_back( 3 );
	modes_using.push_back( 4 );
	modes_using.push_back( 5 );
	NMmin.set_modes( modes_using );
	*/

	pose::Pose pose_work;

	// Try same in AtomTree
	pose_work = pose;

	std::cout << "Starting score!" << std::endl;
	sfxn->show(pose_work);

	ATmin.run( pose_work, mm, *sfxn, options );
	pose_work.dump_pdb( "ATmin.pdb" );

	sfxn->show(pose_work);

	// Same for NMmin
	pose_work = pose;
	NMmin.run( pose_work, mm, *sfxn, options );
	pose_work.dump_pdb( "NMmin.pdb" );

	sfxn->show(pose_work);

}

int main( int argc, char *argv [] ){
  using namespace protocols::normalmode;
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	NEW_OPT( fpd::rsd_changing, "", 0 );
	NEW_OPT( fpd::modes, "", 0 );

	devel::init(argc, argv);

	protocols::moves::MoverOP tocen 
		= new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );

  core::pose::Pose pose;

  core::chemical::ResidueTypeSetCAP rsd_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

  core::import_pose::pose_from_file( pose, *rsd_set, option[ in::file::s ](1) , core::import_pose::PDB_file); 

	//test_NMmin( pose );
	test_NMrelaxer( pose );

	return 0;
}

