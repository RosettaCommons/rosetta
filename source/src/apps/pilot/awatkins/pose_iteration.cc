// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   cov_hbs.cc
/// @brief  Sidechain conjugation to acryl amides
/// @author Andy Watkins (amw579@nyu.edu)

// includes
#include <iostream>
#include <fstream>
#include <string>

#include <devel/init.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>

// Mover headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

int main ( int argc, char* argv[] )
{
	try {
		devel::init(argc, argv);

		using namespace core;
		using namespace pose;
		using namespace chemical;
		using namespace conformation;

		core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueType const & ala = residue_set_cap->name_map( "ALA" );
		ResidueOP res_ala( new Residue( ala, true ) ); 
		Pose pose;
		pose.append_residue_by_jump( *res_ala, 1 );
		pose.append_residue_by_bond( *res_ala, true );
		pose.append_residue_by_bond( *res_ala, true );
		pose.append_residue_by_bond( *res_ala, true );
		pose.append_residue_by_bond( *res_ala, true );
		pose.append_residue_by_bond( *res_ala, true );
		pose.append_residue_by_bond( *res_ala, true );
		pose.append_residue_by_bond( *res_ala, true );
		pose.append_residue_by_bond( *res_ala, true );
		
		for ( Residue const & res : pose ) {
			std::cout << "I am happy about " << res.seqpos() << " and here is why: " << res << std::endl;
		}
		
		// Demonstrate non ranged for syntax
		for ( auto res_iter = pose.begin(), pose_end = pose.end(); res_iter != pose_end; ++res_iter ) {
			// res_iter is a special vector1< ResidueOP >::iterator such that operator* applies	
			// an extra dereference, so it gives us our const & style access.
			std::cout << "I am happy about " << res_iter->seqpos() << " and here is why: " << *res_iter << std::endl;
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
