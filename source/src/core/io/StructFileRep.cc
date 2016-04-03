// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/StructFileRep.cc
/// @brief  Class/structure method definitions for StructFileRep
/// @author Sergey Lyskov
/// @author Andy Watkins
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @author Labonte <JWLabonte@jhu.edu>


// Unit headers
#include <core/io/HeaderInformation.hh>
#include <core/io/StructFileRep.hh>

// Package headers
#include <core/io/Remarks.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>


namespace core {
namespace io {

// Standard Methods ///////////////////////////////////////////////////////////
StructFileRep::StructFileRep() : utility::pointer::ReferenceCount(),
	filename_( "" ),
	modeltag_( "" ),
	header_( HeaderInformationOP( new HeaderInformation() ) ),
	remarks_( RemarksOP( new Remarks ) ),
	heterogen_names_(),
	residue_type_base_names_(),
	ssbond_map_(),
	link_map_(),
	crystinfo_(),
	chains_(),
	foldtree_string_( "" ),
	pdb_comments_(),
	additional_string_output_( "" )
{}

StructFileRep::~StructFileRep()
{}

/// @details Uses the compiler-default copy constructor.
StructFileRepOP
StructFileRep::clone() const {
	return StructFileRepOP( new StructFileRep( *this ) );
}


// Append more string data to the additional_string_output_ string in the SFR.
/// @details Retains all string data already added.
void
StructFileRep::append_to_additional_string_output( std::string const &input_string )
{
	additional_string_output_ = additional_string_output_ + input_string;
}


// Helper Functions ///////////////////////////////////////////////////////////
// Debug printing, serializing to Tracer-like object
std::ostream &
operator<<( std::ostream &os, AtomInformation const & ai )
{
	os << "<AtomInformation>{" << "isHet=" << ai.isHet << " serial=" << ai.serial << " name=" << ai.name << " resName=" << ai.resName
		<< " chainID=" << ai.chainID << " resSeq=" << ai.resSeq
		<< " x=" << ai.x << " y=" << ai.y << " z=" << ai.z
		<< " temperature=" << ai.temperature
		<< " occupancy=" << ai.occupancy
		<< " segmentID=" << ai.segmentID
		<< " element=" << ai.element
		<< "}";
	return os;
}

// Output StructFileRep object to TR-like stream in human-readable format.
std::ostream &
operator<<( std::ostream & os, StructFileRep const & sfr )
{
	os << "<StructFileRep>{";
	for ( uint i( 0 ); i < sfr.chains().size(); ++i ) {
		os << "Chain<" << i << ">";
		for ( uint j( 0 ); j < sfr.chains()[ i ].size(); ++j ) {
			os << "[" << j << ":" << sfr.chains()[ i ][ j ] << "]" << std::endl;
		}
	}
	os << "}";
	return os;
}

} // namespace io
} // namespace core
