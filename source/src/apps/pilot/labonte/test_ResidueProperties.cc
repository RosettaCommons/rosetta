// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file   test_ResidueProperties.cc
/// @brief  This is simple pilot app for testing ResidueProperties.
/// @note   I intend to convert this into unit tests once everything is sound.
/// @author Labonte <JWLabonte@jhu.edu>

// Package headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <devel/init.hh>

// Utility header
#include <utility/vector1.hh>

// C++ headers
#include <iostream>


int
main( int argc, char *argv[] )
{
	using namespace std;
	using namespace utility;
	using namespace core;
	using namespace chemical;
	using namespace orbitals;
	using namespace pose;


	try {
		// Initialize core.
		devel::init( argc, argv );

		ChemicalManager * cm( ChemicalManager::get_instance() );
		string const tag( FA_STANDARD );
		AtomTypeSetCOP atom_types = cm->atom_type_set( tag );
		ElementSetCOP element_types = cm->element_set( "default" );
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set( tag );
		OrbitalTypeSetCOP orbital_types = cm->orbital_type_set( tag );

		ResidueType rsd( atom_types, element_types, mm_atom_types, orbital_types );

		rsd.add_property( "PROTEIN" );

		if ( rsd.is_protein() ) {
			cout << "I'm a protein!" << endl;
		}

		if ( rsd.has_property( "PROTEIN" ) ) {
			cout << "I'm a protein!" << endl;
		}

		if ( rsd.properties().has_property( PROTEIN ) ) {
			cout << "I'm a protein!" << endl;
		}

		rsd.add_property( "POLAR" );
		rsd.add_property( "CHARGED" );

		//cout << rsd.properties().get_list_of_properties() << endl;
		cout << rsd.properties() << endl;

		rsd.add_variant_type( UPPER_TERMINUS_VARIANT );
		if ( rsd.has_variant_type( UPPER_TERMINUS_VARIANT ) ) {
			cout << "I'm an upper terminus!" << endl;
		}

		Pose pose;
		make_pose_from_sequence( pose, "AAA", "fa_standard" );

		cout << pose.residue( 1 ) << endl;

	} catch ( excn::EXCN_Base const & e ) {
		cerr << "Caught exception: " << e.msg() << endl;
		return -1;
	}
	return 0;
}
