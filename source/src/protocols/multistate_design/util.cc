// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistate_design/util.cc
/// @brief
/// @author Justin Ashworth wrote this, but Andrew Leaver-Fay copy-and-pasted it into this file.

// Unit headers
#include <protocols/multistate_design/util.hh>

// Package headers
#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/multistate_design/MultiStatePacker.hh>


// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/pack/task/PackerTask.hh>

// C++ headers
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace multistate_design {

protocols::genetic_algorithm::EntityElements
list_amino_acid_options(
	core::Size i, // aka, residue index
	core::pack::task::ResidueLevelTask const & rtask
)
{
	using namespace core;
	using namespace core::chemical;
	using namespace protocols::genetic_algorithm;

	EntityElements choices;
	// to avoid duplicate AA's (such as for multiple histidine ResidueTypes)
	std::set< core::chemical::AA > aaset;
	std::list< ResidueTypeCOP > const & allowed( rtask.allowed_residue_types() );
	 for ( auto const & t : allowed ) {
		core::chemical::AA aa( t->aa() );
		// avoid duplicate AA's (such as for multiple histidine ResidueTypes)
		if ( aaset.find( aa ) != aaset.end() ) continue;
		aaset.insert(aa);
		//TR(t_debug) << "adding choice " << aa << std::endl;
		choices.push_back( protocols::genetic_algorithm::EntityElementOP( new PosType( i, aa ) ) );
	}
	return choices;
}

protocols::genetic_algorithm::EntityElements
entity_elements_from_1letterstring(
	std::string const & input
)
{
	protocols::genetic_algorithm::EntityElements elements( input.size() );
	for ( core::Size ii = 0, count = 1; ii < input.size(); ++ii, ++count ) {
		std::ostringstream output;
		output << "AA:" << count << ":" << input[ ii ];
		elements[ count ] = protocols::genetic_algorithm::EntityElementFactory::get_instance()->element_from_string( output.str() );
	}
	return elements;
}


}
}

