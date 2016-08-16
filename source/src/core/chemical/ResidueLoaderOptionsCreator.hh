// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/ResidueLoaderOptionsCreator.hh
/// @brief
/// @author Sam DeLuca

#ifndef INCLUDED_core_chemical_ResidueLoaderOptionsCreator_hh
#define INCLUDED_core_chemical_ResidueLoaderOptionsCreator_hh

//package headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

namespace core {
namespace chemical {

class ResidueLoaderOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
{
public:

	virtual std::string options_type() const;
	virtual basic::resource_manager::ResourceOptionsOP create_options() const;

};

}
}

#endif
