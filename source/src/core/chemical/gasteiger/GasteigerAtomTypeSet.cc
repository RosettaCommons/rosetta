// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/gasteiger/GasteigerAtomTypeSet.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>

#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/Element.hh>

// Project headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

namespace core {
namespace chemical {
namespace gasteiger {

static thread_local basic::Tracer tr( "core.chemical.gasteiger.GasteigerAtomTypeSet" );

GasteigerAtomTypeSet::GasteigerAtomTypeSet() {}

GasteigerAtomTypeSet::GasteigerAtomTypeSet( ElementSetCAP element_set ):
	element_set_( element_set )
{}

GasteigerAtomTypeSet::GasteigerAtomTypeSet( GasteigerAtomTypeSet const & other ):
	utility::pointer::ReferenceCount(),
	element_set_( other.element_set_ ),
	atom_type_index_( other.atom_type_index_ ),
	atom_types_( other.atom_types_ )
{}

GasteigerAtomTypeSet::~GasteigerAtomTypeSet() {}

/// @details Initialize a GasteigerAtomTypeSet from an external file "filename",
/// and set parameters and properties for each Atom.
/// Refer to rosetta_database/chemical/gasteiger/atom_types/atom_properties.txt
/// for file format
void
GasteigerAtomTypeSet::read_file( std::string const & filename )
{
	utility::io::izstream data( filename.c_str() );

	if ( !data.good() ) utility_exit_with_message( "Unable to open gasteiger atom type file: "+filename );

	// now parse the rest of the file
	{
		using namespace basic;

		char next;
		std::string line;
		while ( data.good() ) {
			while ( next=data.peek(), next == ' ' || next == '\n' || next == '\t' ) { data.get(); } // Discard leading whitespace
			if ( ! data.good() ) {
				break;
			}
			if ( data.peek() == '#' ) {
				getline( data,line ); // Discard the comment line
				continue;
			}
			GasteigerAtomTypeDataOP atom( new GasteigerAtomTypeData );
			atom->read( data, element_set_);
			if ( data.good() ) {
				atom_types_.push_back( atom );
				std::string symbol( atom->get_name() );
				if ( atom_type_index_.count( symbol ) ) {
					utility_exit_with_message("GasteigerAtomTypeSet:: duplicate atom type symbol "+symbol);
				}
				atom_type_index_[ symbol ] = atom_types_.size();
				tr.Debug << "New Bcl atom type: " << symbol << std::endl;
			}
		}
	}
}

void
GasteigerAtomTypeSet::read_bond_file( std::string const & filename ){

	//setup map needed to parse information
	// make a mapping from element type, # bonds, # electrons in bonds to a listing of all atom types with those
	// properties
	std::map< std::string,  utility::vector1<gasteiger::GasteigerAtomTypeDataOP> > phenotype_to_atom_type;
	for ( utility::vector1< GasteigerAtomTypeDataOP >::iterator itr(atom_types_.begin()); itr != atom_types_.end(); ++itr ) {
		gasteiger::GasteigerAtomTypeDataOP & atom_type( *itr);

		// create the triplet of phenotypic info about this atom type
		const std::string phenotype
			(
			atom_type->get_element_type()->get_chemical_symbol() + utility::to_string(atom_type->get_number_bonds()) +
			utility::to_string(atom_type->get_number_electrons_in_bonds())
		);

		phenotype_to_atom_type[ phenotype].push_back( atom_type);
	}


	utility::io::izstream data( filename.c_str() );
	if ( !data.good() ) utility_exit_with_message( "Unable to open bond file: "+filename );

	// now parse the rest of the file
	{
		using namespace basic;

		char next;
		std::string line;
		while ( data.good() ) {
			while ( next=data.peek(), next == ' ' || next == '\n' || next == '\t' ) { data.get(); } // Discard leading whitespace
			if ( ! data.good() ) {
				break;
			}
			if ( data.peek() == '#' ) {
				getline( data,line ); // Discard the comment line
				continue;
			}


			std::string key; //phenotype key
			data >> key;
			std::string bond_type; //the bond type
			data >> bond_type;
			Real bond_dist; //the distance between bonds
			data >> bond_dist;


			//utility::vector1<gasteiger::GasteigerAtomTypeData> & vector_atype = phenotype_to_atom_type[key];
			gasteiger::GasteigerAtomTypeData::Properties property;
			if ( bond_type == "VdWaalsRadiusCSD" ) {
				property = gasteiger::GasteigerAtomTypeData::VdWaalsRadiusCSD;
			} else if ( bond_type == "CovalentRadiusAromaticBond" ) {
				property = gasteiger::GasteigerAtomTypeData::CovalentRadiusAromaticBond;
			} else if ( bond_type == "CovalentRadiusTripleBond" ) {
				property = gasteiger::GasteigerAtomTypeData::CovalentRadiusTripleBond;
			} else if ( bond_type == "CovalentRadiusDoubleBond" ) {
				property = gasteiger::GasteigerAtomTypeData::CovalentRadiusDoubleBond;
			} else if ( bond_type == "CovalentRadiusSingleBond" ) {
				property =  gasteiger::GasteigerAtomTypeData::CovalentRadiusSingleBond;
			} else {
				utility_exit_with_message("GasteigerAtomTypeSet:: unknown bond type "+bond_type);
			}

			//set the property
			for ( Size i=1; i<=phenotype_to_atom_type[key].size(); ++i ) {
				phenotype_to_atom_type[key][i]->set_property(property, bond_dist);
			}
		}
	}
}


/// @brief Check if there is an element_type associated with an element_symbol string
bool
GasteigerAtomTypeSet::contains_atom_type( std::string const & atom_type_name ) const
{
	std::map< std::string, core::Size >::const_iterator
		iter( atom_type_index_.find( atom_type_name ) );
	return iter != atom_type_index_.end();
}


/// @brief Lookup the element index by the element_symbol string
Size
GasteigerAtomTypeSet::atom_type_index( std::string const & atom_type_name ) const
{
	std::map< std::string, core::Size >::const_iterator
		iter( atom_type_index_.find( atom_type_name ) );
	if ( iter == atom_type_index_.end() ) {
		utility_exit_with_message("unrecognized atom_type_name "+atom_type_name);
	}
	return iter->second;
}

/// @brief Lookup the element index by the element_symbol string
GasteigerAtomTypeDataCOP
GasteigerAtomTypeSet::atom_type( std::string const & atom_type_name ) const
{
	return atom_types_[ atom_type_index( atom_type_name ) ];
}


/// @brief Lookup an Atom by 1-based indexing
GasteigerAtomTypeDataCOP
GasteigerAtomTypeSet::operator[] ( Size const index ) const
{
	return atom_types_[ index ];
}

/// @brief Return the associated element type set.
ElementSetCAP
GasteigerAtomTypeSet::element_set() const {
	return element_set_;
}

/// @brief Return the type that's assigned for fake atoms. (Virtual atoms and the like.)
GasteigerAtomTypeDataCOP
GasteigerAtomTypeSet::type_for_fake_atoms() const {
	return atom_types_[ atom_type_index( "FAKE" ) ];
}


} // gasteiger
} // chemical
} // core
