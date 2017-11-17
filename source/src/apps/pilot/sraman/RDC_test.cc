// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

//core
#include <core/types.hh>
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/sequence/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

//devel
#include <devel/residual_dipolar_coupling/RDC_main.hh>
#include <devel/residual_dipolar_coupling/RDC_energies.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>
 #include <map>
#include <vector>


// option key includes

#include <basic/options/keys/RDC.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/looprelax.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::sequence;
	using namespace core::scoring;

	// Parses command line options and inits RNG.
	devel::init(argc, argv);

	std::string filename( option[ OptionKeys::RDC::RDC_file ]().name() );
	std::string sequence( core::sequence::read_fasta_file( option[ in::file::fasta ]() )[1]);
	std::cout << "Input file name " << filename << std::endl;
	//	devel::residual_dipolar_coupling::read_RDC_file( filename );
	//	devel::residual_dipolar_coupling::RDC_data my_RDC_data( filename );

	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	core::pose::Pose pose, extended_pose;
	core::import_pose::pose_from_file( pose, *rsd_set, option [ OptionKeys::looprelax::input_pdb ]().name() , core::import_pose::PDB_file);
	core::pose::make_pose_from_sequence( 	extended_pose, sequence,
		*( core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" ))
	);

	// make extended chain
	for ( core::Size pos = 1; pos <= extended_pose.size(); pos++ ) {
		extended_pose.set_phi( pos, -45 );
		extended_pose.set_psi( pos, -45 );
		extended_pose.set_omega( pos, 180 );
	}


	core::scoring::ScoreFunctionOP scorefxn( new ScoreFunction() );
	scorefxn->set_weight( vdw, 1.0 );
	scorefxn->set_weight( env, 1.0 );
	scorefxn->set_weight( cbeta, 1.0 );
	scorefxn->set_weight( pair, 1.0 );


	//	std::map< core::Size, utility::vector1<devel::residual_dipolar_coupling::RDC> > RDC_data_lines( my_RDC_data.get_RDC_data() );


	//v	utility::vector1< devel::residual_dipolar_coupling::RDC > All_RDC_lines( devel::residual_dipolar_coupling::read_RDC_file( filename ) );
	//v	devel::residual_dipolar_coupling::eval_dipolar( pose, All_RDC_lines );
	//v	devel::residual_dipolar_coupling::eval_dipolar( extended_pose, All_RDC_lines );


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


