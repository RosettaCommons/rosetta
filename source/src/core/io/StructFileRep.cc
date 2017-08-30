// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

/// @brief Do these link information correspond to the same connection (ignoring order)
bool same_link( LinkInformation const & l1, LinkInformation const & l2 ) {
	// Order quick compares (char/int) before string compares
	if ( l1.resSeq1 == l2.resSeq1 &&
			l1.chainID1 == l2.chainID1 &&
			l1.iCode1 == l2.iCode1 &&
			l1.resSeq2 == l2.resSeq2 &&
			l1.chainID2 == l2.chainID2 &&
			l1.iCode2 == l2.iCode2 &&
			l1.name1 == l2.name1 &&
			l1.resName1 == l2.resName1 &&
			l1.name2 == l2.name2 &&
			l1.resName2 == l2.resName2 ) {
		return true;
	}
	// Test flipped.
	if ( l1.resSeq1 == l2.resSeq2 &&
			l1.chainID1 == l2.chainID2 &&
			l1.iCode1 == l2.iCode2 &&
			l1.resSeq2 == l2.resSeq1 &&
			l1.chainID2 == l2.chainID1 &&
			l1.iCode2 == l2.iCode1 &&
			l1.name1 == l2.name2 &&
			l1.resName1 == l2.resName2 &&
			l1.name2 == l2.name1 &&
			l1.resName2 == l2.resName1 ) {
		return true;
	}
	return false;
}

bool link_in_vector( utility::vector1< LinkInformation > const & link_vector, LinkInformation const & link ) {
	for ( LinkInformation const & vect_link: link_vector ) {
		if ( same_link( vect_link, link ) ) {
			return true;
		}
	}
	return false;
}

// Standard Methods ///////////////////////////////////////////////////////////
StructFileRep::StructFileRep() : utility::pointer::ReferenceCount(),
	filename_( "" ),
	modeltag_( "" ),
	header_( HeaderInformationOP( new HeaderInformation() ) ),
	remarks_( RemarksOP( new Remarks ) ),
	heterogen_names_(),
	residue_type_base_names_(),
	HELIXInformations_(),
	SHEETInformations_(),
	ssbond_map_(),
	link_map_(),
	crystinfo_(),
	chains_(),
	foldtree_string_( "" ),
	pdb_comments_(),
	additional_string_output_( "" )
{}

StructFileRep::~StructFileRep() = default;

/// @details Uses the compiler-default copy constructor.
StructFileRepOP
StructFileRep::clone() const {
	return StructFileRepOP( new StructFileRep( *this ) );
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

/// @brief Append more string data to the additional_string_output_ string in the SFR.
/// @details Retains all string data already added.
void
StructFileRep::append_to_additional_string_output(
	std::string const &input_string
) {
	additional_string_output_ = additional_string_output_ + input_string;
}

/// @brief Debugging output for LinkInformation
std::ostream & operator<<( std::ostream & os, LinkInformation const & li ) {
	// Debugging-style output - we're not overly concerned with column alignment
	os << "LINK ";
	os << li.name1 << " " << li.resName1 << " " << li.chainID1 << " " << li.resSeq1 << li.iCode1 << "   ";
	os << li.name2 << " " << li.resName2 << " " << li.chainID2 << " " << li.resSeq2 << li.iCode2 << "   ";
	os << li.length;
	return os;
}

} // namespace io
} // namespace core
