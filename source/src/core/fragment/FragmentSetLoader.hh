// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/FragmentSetLoader.hh
/// @brief
/// @author

#ifndef INCLUDED_core_fragment_FragmentSetLoader_HH
#define INCLUDED_core_fragment_FragmentSetLoader_HH

//unit headers
#include <core/fragment/FragmentSetLoader.fwd.hh>

//project headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

//C++ headers
#include <istream>

namespace core {
namespace fragment {

/// @brief %FragmentSetLoader constructs a FragSet instance from data provided by the %ResourceManager.
class FragmentSetLoader : public basic::resource_manager::ResourceLoader
{
public:
	/// @brief Construct the %FragmentSetLoader.
	FragmentSetLoader();

	/// @brief Destructor.
	~FragmentSetLoader() override;

	/// @brief Return a FragmentFileDataOP constructed from the given input stream (istream).

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

} // namespace fragment
} // namespace core

#endif //INCLUDED_core_fragment_FragmentSetLoader_HH
