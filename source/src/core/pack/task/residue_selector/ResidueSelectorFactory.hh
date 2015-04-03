// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ResidueSelectorFactory.hh
/// @brief  Class for instantiating arbitrary ResidueSelectors from a string --> ResidueSelectorCreator map
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_task_residue_selector_ResidueSelectorFactory_HH
#define INCLUDED_core_pack_task_residue_selector_ResidueSelectorFactory_HH

// Package headers
#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>
#include <core/pack/task/residue_selector/ResidueSelectorCreator.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>

// C++ headers
#include <map>
#include <string>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

class ResidueSelectorFactory {
public:
	static ResidueSelectorFactory * get_instance();

	void factory_register( ResidueSelectorCreatorOP creator );
	bool has_type( std::string const & ) const;

	ResidueSelectorOP new_residue_selector(
		std::string const & selector_name,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) const;

	void set_throw_on_double_registration();

private:
	static ResidueSelectorFactory * instance_;
	ResidueSelectorFactory();

private:
	std::map< std::string, ResidueSelectorCreatorOP > creator_map_;
	bool throw_on_double_registration_;
};


} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif
