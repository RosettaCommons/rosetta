// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   ResidueSelector.cc
/// @author Phil Bradley


// Unit headers
#include <core/chemical/ResidueSelector.hh>

// Package headers
#include <core/chemical/ResidueTypeSet.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>

// Utility headers
#include <utility/vector1.hh>


namespace core {
namespace chemical {

/// @details Auto-generated virtual destructor
ResidueSelector::~ResidueSelector() {}

/// @details Auto-generated virtual destructor
ResidueSelectorSingle::~ResidueSelectorSingle() {}


////////////////////////////////////////////////////////////////////////////////////////

ResidueSelectorSingleOP
residue_selector_single_from_line( std::string const & line )
{
	std::istringstream l( line );
	std::string tag;
	l >> tag;
	if ( l.fail() ) return 0;

	bool desired_result( true );
	if ( tag == "NOT" ) {
		desired_result = false;
		l >> tag;
	}

	if ( tag == "AA" ) {
		AA aa;
		utility::vector1< AA > aas;
		l >> aa;
		while ( !l.fail() ) {
			aas.push_back( aa );
			l >> aa;
		}
		//std::cout << "AAline: " << aas.size() << ' ' << line << std::endl;
		if ( !aas.empty() ) return ResidueSelectorSingleOP( new Selector_AA( aas, desired_result ) );

	} else if ( tag == "NAME3" ) {
		std::string name3;
		utility::vector1< std::string > name3s;
		l >> name3;
		while ( !l.fail() ) {
			name3s.push_back( name3 );
			l >> name3;
		}
		if( !name3s.empty() )  return ResidueSelectorSingleOP( new Selector_NAME3( name3s, desired_result ) );

	} else if ( tag == "PROPERTY" ) {
		std::string property;
		utility::vector1< std::string > properties;
		l >> property;
		while ( !l.fail() && property[0] != '#' ) {
			properties.push_back( property );
			l >> property;
		}
		if ( !properties.empty() ) return ResidueSelectorSingleOP( new Selector_PROPERTY( properties, desired_result ) );

	} else if ( tag == "VARIANT_TYPE" ) {
		std::string variant_type;
		utility::vector1< std::string > variant_types;
		l >> variant_type;
		while ( !l.fail() && variant_type[0] != '#' ) {
			variant_types.push_back( variant_type );
			l >> variant_type;
		}
		if ( !variant_types.empty() ) return ResidueSelectorSingleOP( new Selector_VARIANT_TYPE( variant_types, desired_result ) );

	} else if (tag == "UPPER_POSITION") {  // This is the position label at which the upper connection is attached.
		uint position;
		l >> position;
		if (position) {
			return ResidueSelectorSingleOP( new Selector_UPPER_POSITION(position, desired_result) );
		}

	} else if ( tag == "CMDLINE_SELECTOR" ) {
		std::string selector_string;
		l >> selector_string; // if one wants AND logical operation make this a vector of selector_strings...
		if ( !selector_string.empty() ) return ResidueSelectorSingleOP( new Selector_CMDFLAG( selector_string, desired_result ) );
	}
	std::cout << "residue_selector_single: unrecognized line: " << line << std::endl;
	return 0;
}

///
ResidueTypeCOPs
ResidueSelector::select( ResidueTypeSet const & rsd_set )
{
	ResidueTypeCOPs rsd_list;
	for ( ResidueTypeCOPs::const_iterator it= rsd_set.residue_types().begin(), ite= rsd_set.residue_types().end();
				it != ite; ++it ) {
		if ( operator[]( **it ) ) rsd_list.push_back( *it );
	}
	return rsd_list;
}

Selector_CMDFLAG::Selector_CMDFLAG(std::string  const & flag_in, bool const result) : ResidueSelectorSingle( result )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	b_flag_is_present_ = false;
	if ( !option[ OptionKeys::chemical::patch_selectors ].user() ) return;
	for ( StringVectorOption::const_iterator it = option[ OptionKeys::chemical::patch_selectors ]().begin(),
				eit = option[ OptionKeys::chemical::patch_selectors ]().end(); it != eit ; ++ it ) {
		if ( *it == flag_in ) {
			b_flag_is_present_ = true;
			break;
		}
	}
}

} // chemical
} // core
