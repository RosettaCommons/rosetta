// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/multistate_design/util.cc
/// @brief
/// @author Justin Ashworth wrote this, but Andrew Leaver-Fay copy-and-pasted it into this file.

// Unit headers
#include <protocols/multistate_design/util.hh>

// Package headers
#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/multistate_design/MultiStatePacker.hh>
#include <protocols/multistate_design/SingleState.hh> // REQUIRED FOR WINDOWS

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh> // REQUIRED FOR WINDOWS 

// C++ headers
#include <set>

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
	std::list< ResidueTypeCAP > const & allowed( rtask.allowed_residue_types() );
	for ( std::list< ResidueTypeCAP >::const_iterator t( allowed.begin() ), end( allowed.end() );
				t != end; ++t ) {
		core::chemical::AA aa( (*t)->aa() );
		// avoid duplicate AA's (such as for multiple histidine ResidueTypes)
		if ( aaset.find( aa ) != aaset.end() ) continue;
		aaset.insert(aa);
		//TR(t_debug) << "adding choice " << aa << std::endl;
		choices.push_back( new PosType( i, aa ) );
	}
	return choices;
}

}
}

