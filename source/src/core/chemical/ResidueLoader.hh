// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/ResidueLoader.fwd.hh
/// @brief headers for the Residue Loader
/// @author Sam DeLuca

#ifndef INCLUDED_core_chemical_ResidueLoader_hh
#define INCLUDED_core_chemical_ResidueLoader_hh

#include <core/chemical/ResidueLoader.fwd.hh>
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <istream>

namespace core {
namespace chemical {

class ResidueLoader : public basic::resource_manager::ResourceLoader
{
public:
	ResidueLoader() {}
	~ResidueLoader() override = default;

	/// @brief returns a pointer to a ResidueType object originated from the data_source specified to the ResourceManager
	
	utility::pointer::ReferenceCountOP
	create_resource(
		basic::resource_manager::ResourceOptions const & options,
		basic::resource_manager::LocatorID const & locator_id,
		std::istream & istream
	) const override;

	
	basic::resource_manager::ResourceOptionsOP
	default_options() const override;

};

}
}

#endif

