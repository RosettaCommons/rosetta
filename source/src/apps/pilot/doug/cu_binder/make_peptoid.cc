// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/doug/peptoid_rotlibs/peptoid_rotlib_test1.cc
/// @brief Test the cyclic peptide/peptoid patches
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// core headers
#include <core/init.hh>
#include <core/types.hh>

#include <core/pose/Pose.hh>

#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResiduePeptoidLibrary.hh>
#include <core/pack/rotamers/RotamericSingleResiduePeptoidLibrary.hh>
#include <core/pack/rotamers/RotamericSingleResiduePeptoidLibrary.tmpl.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/util.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <core/graph/Graph.hh>

#include <core/kinematics/MoveMap.hh>

// protocols headers
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <protocols/simple_moves/MinMover.hh>

// basic headers
#include <basic/database/open.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// utility headers
#include <utility/io/izstream.hh>

// c++ headers
#include <iostream>
#include <string>

// namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

int
main( int argc, char * argv [] )
{

	// init options, rng, etc.
	core::init(argc, argv);

	std::cout << "SCORE FUNCTION" << std::endl;

// create score function
	core::scoring::ScoreFunctionOP score_fxn( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::MM_STD_WTS ) );
	score_fxn->set_weight( unfolded, 0.0 );

	// get a ResidueTypeSet
	std::cout << "RTS" << std::endl;
	ResidueTypeSetCAP rsd_type_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	std::cout << "NAME LIST" << std::endl;

	// load array with names
	utility::vector1<std::string> aa_names;
	aa_names.push_back( "SAR_p:NtermPeptoidFull" );
	aa_names.push_back( "P28" );
	aa_names.push_back( "SAR" );
	aa_names.push_back( "P28" );
	aa_names.push_back( "SAR" );
	aa_names.push_back( "P28_p:CtermPeptoidFull" );

	// turn array in to pose
	std::cout << "BUILD POSE" << std::endl;
	pose::Pose pose;

	for ( utility::vector1<std::string>::const_iterator i( aa_names.begin() ), end( aa_names.end() ); i != end; ++i ) {

		ResidueType const & rsd_type( rsd_type_set->name_map( *i ) );

		Residue rsd( rsd_type, true );

		if ( i == aa_names.begin() ) {
			pose.append_residue_by_jump( rsd, 1 );
		} else {
			pose.append_residue_by_bond( rsd, true );
		}
	}

	pose.conformation().detect_bonds();
	pose.conformation().show_residue_connections();

	std::string before_filename( "before.pdb" );
	pose.dump_scored_pdb( before_filename, *score_fxn );

	std::cout << "MINIMIZING" << std::endl;

	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;

	movemap->set_bb( true );
	movemap->set_chi( true );

	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( movemap, score_fxn, basic::options::option[ basic::options::OptionKeys::run::min_type ].value(), 0.01, true );

	min_mover->apply( pose );

	std::string after_min_filename( "after_min.pdb" );
	pose.dump_scored_pdb( after_min_filename, *score_fxn );

	std::cout << "DONE" << std::endl;

	return 0;
}
