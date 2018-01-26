// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/FragmentSetLoaderCreator.hh
/// @brief
/// @author

#ifndef INCLUDED_core_fragment_FragmentSetLoaderCreator_HH
#define INCLUDED_core_fragment_FragmentSetLoaderCreator_HH

//unit headers
#include <basic/resource_manager/ResourceLoaderCreator.hh>

//project headers

//utility headers


//C++ headers

namespace core {
namespace fragment {

/// @brief %FragmentSetLoaderCreator allows the ResourceLoaderFactory to create a FragmentSetLoader instance,
/// which in turn can load a FragmentSet.
/// @details The FragmentSetLoader class can be constructed from the string "FragmentSet", which enables a user to specify
/// that this type of %resource is required for a particular %job in their XML input file.
class FragmentSetLoaderCreator : public basic::resource_manager::ResourceLoaderCreator
{
public:

	/// @brief Return a up-casted owning pointer (ResourceLoaderOP) to the resource loader.
	basic::resource_manager::ResourceLoaderOP
	create_resource_loader() const override;

	/// @brief Return the string identifier for the associated ResourceLoader (FragmentSet).
	std::string loader_type() const override;

	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;


};

} // namespace fragment
} // namespace core


#endif //INCLUDED_core_fragment_FragmentSetLoaderCreator_HH
