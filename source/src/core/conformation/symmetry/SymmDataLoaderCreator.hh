// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/SymmDataLoaderCreator.hh
/// @brief  load the SymmData data-structure, which is used to configure symmetric poses.
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_core_conformation_symmetry_SymmDataLoaderCreator_hh
#define INCLUDED_core_conformation_symmetry_SymmDataLoaderCreator_hh

//unit headers
#include <basic/resource_manager/ResourceLoaderCreator.hh>


namespace core {
namespace conformation {
namespace symmetry {

class SymmDataLoaderCreator : public basic::resource_manager::ResourceLoaderCreator
{
public:

	basic::resource_manager::ResourceLoaderOP
	create_resource_loader() const override;

	std::string loader_type() const override;

	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} // namespace
} // namespace
} // namespace

#endif // include guard
