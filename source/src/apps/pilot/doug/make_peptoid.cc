// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/doug/peptoid_rotlibs/peptoid_rotlib_test1.cc
/// @brief Test the cyclic peptide/peptoid patches
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// core headers
#include <devel/init.hh>
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

#include <utility/graph/Graph.hh>

#include <core/kinematics/MoveMap.hh>

// protocols headers
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/PyMOLMover.hh>

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

//local options
basic::options::StringOptionKey const peptoid_tlc( "peptoid_tlc" );

int
main( int argc, char * argv [] )
{
	try {
		// add local options
		option.add( peptoid_tlc, "peptoid three letter code" ).def("01");

		// init options, rng, etc.
		devel::init(argc, argv);

		std::cout << "SCORE FUNCTION" << std::endl;

		// create score function
		core::scoring::ScoreFunctionOP score_fxn( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::MM_STD_WTS ) );
		score_fxn->set_weight( unfolded, 0.0 );

		// get a ResidueTypeSet
		std::cout << "RTS" << std::endl;
		ResidueTypeSetCOP rsd_type_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

		// make a list of names
		std::cout << "NAME LIST" << std::endl;

		std::string name( option[ peptoid_tlc ].value() );
		//std::stringstream nterm_name; nterm_name << name << "_p:NtermPeptoidFull";
		std::stringstream nterm_name; nterm_name << name << "_p:AcetylatedPeptoidNterm";
		std::stringstream cterm_name; cterm_name << name << "_p:CtermPeptoidFull";

		// load array with names
		utility::vector1<std::string> aa_names;
		aa_names.push_back( nterm_name.str() );
		aa_names.push_back( name );
		aa_names.push_back( cterm_name.str() );

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

		// pymol mover
		protocols::moves::PyMOLMoverOP pmm( new protocols::moves::PyMOLMover() );
		pmm->keep_history( true );
		pmm->apply( pose );

		// rot precedding omg
		Real orig_pomg( pose.omega(1) );
		for ( Real i(0); i <= 360; i += 10 ) {
			pose.set_omega( 1, orig_pomg + i );
			pmm->apply( pose );
		}

		// rot phi
		Real orig_phi( pose.phi(2) );
		for ( Real i(0); i <= 360; i += 10 ) {
			pose.set_phi( 2, orig_phi + i );
			pmm->apply( pose );
		}

		// rot psi
		Real orig_psi( pose.psi(2) );
		for ( Real i(0); i <= 360; i += 10 ) {
			pose.set_psi( 2, orig_psi + i );
			pmm->apply( pose );
		}

		// rot omg
		Real orig_omg( pose.omega(2) );
		for ( Real i(0); i <= 360; i += 10 ) {
			pose.set_omega( 2, orig_omg + i );
			pmm->apply( pose );
		}

		// rot chis
		if ( pose.residue(2).type().nchi() >= 1 ) {
			for ( Size j(1); j <=  pose.residue(2).type().nchi(); ++j ) {
				Real orig_chi( pose.chi( j, 2 ) );
				for ( Real i(0); i <= 360; i += 10 ) {
					pose.set_chi( j, 2, orig_chi + i );
					pmm->apply( pose );
				}
			}
		}

		// dump pose to file
		std::stringstream filename;
		filename << name << ".pdb";
		pose.dump_scored_pdb( filename.str(), *score_fxn );

		// std::cout << "MINIMIZING" << std::endl;

	// core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;

	// movemap->set_bb( true );
	// movemap->set_chi( true );

	// protocols::minimization_packing::MinMoverOP min_mover = new protocols::minimization_packing::MinMover( movemap, score_fxn, basic::options::option[ basic::options::OptionKeys::run::min_type ].value(), 0.01, true );

	// min_mover->apply( pose );

	// std::string after_min_filename( "after_min.pdb" );
	// pose.dump_scored_pdb( after_min_filename, *score_fxn );

		std::cout << "************************************d**o**n**e***********************************" << std::endl;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
