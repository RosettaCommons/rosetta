// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/silent/SilentFileOptions.hh
/// @brief  Options for constructing a pose from a silent file
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_core_io_silent_SilentFileOptions_hh
#define INCLUDED_core_io_silent_SilentFileOptions_hh


// Unit Headers
#include <core/io/silent/SilentFileOptions.fwd.hh>
#include <basic/resource_manager/ResourceOptions.hh>

// Platform Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <string>

namespace core {
namespace io {
namespace silent {

class SilentFileOptions : public basic::resource_manager::ResourceOptions
{
public:
	SilentFileOptions();

	SilentFileOptions(
		std::string const & name);

	/* Undefined commention out to fix PyRosetta build SilentFileOptions(
	std::string const & name,
	Real mapreso,
	Real grid_spacing); */

	~SilentFileOptions();

	SilentFileOptions(
		SilentFileOptions const & src);

	std::string
	get_silent_struct_type() const;

	void
	set_silent_struct_type(
		std::string const & silent_struct_type
	);

public: // The ResourceOptions public interface
	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag
	);

	/// @brief The class name for a particular ResourceOptions instance.
	/// This function allows for better error message delivery
	virtual
	std::string
	type() const { return "SilentFileOptions"; }

private:
	std::string silent_struct_type_;

};


} // namespace
} // namespace
} // namespace


#endif
