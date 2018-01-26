// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/SymmDataLoader.hh
/// @brief  load the SymmData data-structure, which is used to configure symmetric poses.
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_core_conformation_symmetry_SymmDataLoader_hh
#define INCLUDED_core_conformation_symmetry_SymmDataLoader_hh

//unit headers
#include <core/conformation/symmetry/SymmDataLoader.fwd.hh>
#include <basic/resource_manager/ResourceLoader.hh>

//package headers
#include <basic/resource_manager/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

//C++ headers
#include <istream>

namespace core {
namespace conformation {
namespace symmetry {

class SymmDataLoader : public basic::resource_manager::ResourceLoader
{
public:
	SymmDataLoader();
	virtual ~SymmDataLoader();

	basic::resource_manager::ResourceCOP
	create_resource(
		basic::resource_manager::ResourceManager & resource_manager,
		utility::tag::TagCOP resource_tag,
		std::string const & input_id,
		std::istream & input_stream
	) const override;

	static
	std::string
	classname();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};

} // namespace
} // namespace
} // namespace

#endif // include guard
