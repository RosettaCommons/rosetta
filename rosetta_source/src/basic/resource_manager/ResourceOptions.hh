// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceOptions.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceOptions_hh
#define INCLUDED_basic_resource_manager_ResourceOptions_hh

//unit headers

//project headers

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//C++ headers
#include <string>

namespace basic {
namespace resource_manager {

class ResourceOptions : public utility::pointer::ReferenceCount
{
public:
	ResourceOptions();

	ResourceOptions(
		std::string const & name
	);

	virtual ~ResourceOptions();

	virtual
	void
	parse_my_tag(
		utility::tag::TagPtr tag
	) = 0;

	/// @brief The class name for a particular ResourceOptions instance.
	/// This function allows for better error message delivery
	virtual
	std::string
	type() const = 0;

	/// @brief A name given to a particular ResourceOptions instance.
	/// This function allows for better error message delivery
	std::string
	name() const;

	void name( std::string const & setting );

private:
	std::string name_;

};

} // namespace resource_manager
} // namespace basic



#endif //INCLUDED_basic_resource_manager_ResourceOptions_hh
