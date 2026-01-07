// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/DockingPartners.cc
/// @brief  Simple utility class
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/pose/DockingPartners.hh>

#include <basic/Tracer.hh>

#include <utility/string_util.hh>


namespace core {
namespace pose {

static basic::Tracer TR( "core.pose.DockingPartners" );

DockingPartners
DockingPartners::docking_partners_from_string(std::string const & partner_string) {
	DockingPartners partners;

	utility::vector1< std::string > const split = utility::string_split( partner_string, '_' );

	if ( split.size() != 2 ) {
		TR.Error << "DockingPartners only works with two-body specifications. Designation `" << partner_string << "` implies other than two-body docking." << std::endl;
		utility_exit_with_message( "Can only work with two-body specification. Can't understand " + partner_string );
	}

	// Right now only understand single letter format for chains, when parsed from the partner_string format.
	for ( char c: split[1] ) {
		partners.partner1.push_back( std::string{c} );
	}
	for ( char c: split[2] ) {
		partners.partner2.push_back( std::string{c} );
	}

	return partners;
}

bool
DockingPartners::is_empty() const {
	return partner1.empty() && partner2.empty();
}

bool
DockingPartners::has_both() const {
	return (!partner1.empty()) && (!partner2.empty());
}

std::string
DockingPartners::str() const {
	return partner1_str() + "_" + partner2_str();
}

std::string
DockingPartners::partner1_str() const {
	std::string out = "";
	for ( std::string const & c: partner1 ) {
		if ( c.size() > 1 ) {
			TR.Warning << "Chain `" << c << "` has multiple letters. DockingPartners::partner_str() will be incorrect." << std::endl;
		}
		out += c;
	}
	return out;
}

std::string
DockingPartners::partner2_str() const {
	std::string out = "";
	for ( std::string const & c: partner2 ) {
		if ( c.size() > 1 ) {
			TR.Warning << "Chain `" << c << "` has multiple letters. DockingPartners::partner_str() will be incorrect." << std::endl;
		}
		out += c;
	}
	return out;
}

bool
DockingPartners::operator<( DockingPartners const & other ) const {
	return partner1 < other.partner1 || (partner1 == other.partner1 && partner2 < other.partner2);
}

bool
DockingPartners::operator==( DockingPartners const & other ) const {
	return partner1 == other.partner1 && partner2 == other.partner2;
}

std::ostream &
operator<<( std::ostream & output, DockingPartners const & object_to_output ) {
	output << object_to_output.str();
	return output;
}

} // pose
} // core
