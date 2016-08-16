// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

/// @brief The %ResourceOptions class is responsible for describing all
/// the data requried for instanting a particular resource, except for the
/// data stream (i.e. the file) that contains the data for the resource.
/// For example, when reading in a PDB file, there are 30 different options
/// for how that PDB file should be turned into a Pose.  That data is now
/// held in an ImportPoseOptions object.  The purpose of the %ResourceOptions
/// class is to allow different jobs to load resources in different ways, or
/// for one job to load two different resources of the same type in different
/// ways. For example, a protocol may need both a centroid pose and a
/// fullatom pose to be loaded in from disk; however, if the logic for
/// loading a pose in from disk is controlled by the options system alone,
/// this becomes impossible.
class ResourceOptions : public utility::pointer::ReferenceCount
{
public:
	ResourceOptions();

	/// @brief Assign a name to an instance of the resource options object.
	/// Usefull for identifying flaws in input files defining ResourceOptions.
	ResourceOptions(
		std::string const & name
	);

	virtual ~ResourceOptions();

	/// @brief Describe this instance to a given output stream
	virtual
	void
	show(
		std::ostream & out) const;

	/// @brief Friend output-operator function that invokes the show() function
	friend
	std::ostream &
	operator<< (
		std::ostream & out,
		const ResourceOptions & resource_manager);


	/// @brief Initialize from the recursive "tag" structure.
	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag
	) = 0;

	/// @brief The class name (its type) for a particular ResourceOptions instance.
	/// This function allows for better error message delivery.
	virtual
	std::string
	type() const = 0;

	/// @brief A name given to a particular ResourceOptions instance.
	/// This function allows for better error message delivery.
	std::string
	name() const;

	/// @brief Set the name for this %ResoureOptions instance.
	void name( std::string const & setting );

private:
	std::string name_;

};

} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_manager_ResourceOptions_hh
