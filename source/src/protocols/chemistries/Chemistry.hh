// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/Chemistry.hh
/// @brief  Base class for chemical modifications such as add atoms, delete atoms, replace atom, hydrogen manipulationsds
/// @author Steven Combs
/// @author Rocco Moretti

#ifndef INCLUDED_protocols_chemistries_Chemistry_hh
#define INCLUDED_protocols_chemistries_Chemistry_hh

#include <protocols/chemistries/Chemistry.fwd.hh>
#include <core/chemical/modifications/ChemistryBase.hh>

#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/AtomRefMapping.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <string>
#include <map>

namespace protocols {
namespace chemistries {

class Chemistry : public core::chemical::modifications::ChemistryBase
{
public:

	Chemistry(std::string const & name):
		ChemistryBase( name )
	{}

	/// @brief Modify the passed ResidueType
	/// @details if you have a Pose, call the one which takes the Pose context
	void apply(core::chemical::MutableResidueType & ) override = 0;

	/// @brief Modify the passed ResidueType, context sensitive
	/// @details By default, not context sensitive - will just redirect to the other function.
	virtual void apply(core::chemical::MutableResidueType & restype, core::pose::Pose const & ) {
		apply( restype );
	}

	/// @brief Initialize any data members of this instance from an input tag
	/// and a DataMap object
	virtual void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &) = 0;

private:

};


}
}

#endif
