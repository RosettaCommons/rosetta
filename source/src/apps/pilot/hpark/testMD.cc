// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/md/CartesianMD.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/kinematics/MoveMap.hh>

#include <devel/init.hh>

int main( int argc, char *argv [] ){
  using namespace protocols::md;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	devel::init(argc, argv);

  core::pose::Pose pose;

	core::scoring::ScoreFunctionCOP sfxnOP 
		= core::scoring::ScoreFunctionFactory::create_score_function( option[score::weights] );

  core::chemical::ResidueTypeSetCAP rsd_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

  core::import_pose::pose_from_file( pose, *rsd_set, option[ in::file::s ](1) , core::import_pose::PDB_file); 

	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_bb( true ); movemap->set_chi( true ); 	movemap->set_jump( true );
	movemap->set( core::id::THETA, true ); movemap->set( core::id::D, true );

	CartesianMD MD( pose, sfxnOP, movemap );
	MD.use_rattle( true );
	MD.set_reportstep( 10 );

	core::scoring::ScoreFunctionOP obj_sfxnOP 
		= core::scoring::ScoreFunctionFactory::create_score_function( "score0" );
	obj_sfxnOP->set_weight( core::scoring::vdw, 0.0 );
	obj_sfxnOP->set_weight( core::scoring::coordinate_constraint, 1.0 );
	obj_sfxnOP->set_weight( core::scoring::atom_pair_constraint, 1.0 );

	MD.set_scorefxn_obj( obj_sfxnOP );
	MD.set_report_scorecomp( true );

	if( !option[ in::file::md_schfile ].user() ){
		MD.set_nstep( 500 );
		MD.set_temperature( 300.0 );
	}

	Size nstruct( option[ out::nstruct ]() );

	for( Size i = 1; i <= nstruct; ++i ){
		pose::Pose pose_work( pose );
		MD.apply( pose_work );

		std::stringstream pdbname;
		pdbname << option[ out::prefix ]() << "." << i << ".final.pdb";
		pose_work.dump_pdb( pdbname.str() );
	}

	return 0;
}

