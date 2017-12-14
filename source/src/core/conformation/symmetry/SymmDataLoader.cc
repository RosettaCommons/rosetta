// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/SymmDataLoader.cc
/// @brief  load the SymmData data-structure, which is used to configure symmetric poses.
/// @author Matthew O'Meara (mattjomeara@gmail.com)

//unit headers
#include <core/conformation/symmetry/SymmDataLoader.hh>
#include <core/conformation/symmetry/SymmDataLoaderCreator.hh>
#include <core/conformation/symmetry/SymmDataOptionsCreator.hh>
#include <core/conformation/symmetry/SymmData.hh>

//package headers
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>

//C++ headers
#include <istream>

namespace core {
namespace conformation {
namespace symmetry {


basic::resource_manager::ResourceOP
SymmDataLoader::create_resource(
	basic::resource_manager::ResourceOptions const &,
	basic::resource_manager::LocatorID const &,
	std::istream & istream
) const {
	SymmDataOP symm_data( new SymmData() );
	symm_data->read_symmetry_data_from_stream(istream);
	return symm_data;
}


///// SymmDataOptionsCreator /////
SymmDataOptionsCreator::SymmDataOptionsCreator() = default;

SymmDataOptionsCreator::~SymmDataOptionsCreator() = default;

basic::resource_manager::ResourceOptionsOP
SymmDataOptionsCreator::create_options() const {
	return basic::resource_manager::ResourceOptionsOP( new SymmDataOptions );
}

std::string
SymmDataOptionsCreator::options_type() const {
	return "SymmDataOptions";
}

//// SymmDataLoaderCreator
basic::resource_manager::ResourceLoaderOP
SymmDataLoaderCreator::create_resource_loader() const
{
	return basic::resource_manager::ResourceLoaderOP( new SymmDataLoader() );
}

std::string SymmDataLoaderCreator::loader_type() const
{
	return "SymmData";
}


} // namespace
} // namespace
} // namespace
