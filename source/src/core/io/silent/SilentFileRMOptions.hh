// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/silent/SilentFileRMOptions.hh
/// @brief  Options initialized by the ResourceManager for constructing a pose from a silent file
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_core_io_silent_SilentFileRMOptions_hh
#define INCLUDED_core_io_silent_SilentFileRMOptions_hh


// Unit Headers
#include <core/io/silent/SilentFileRMOptions.fwd.hh>
#include <basic/resource_manager/ResourceOptions.hh>

// Project Headers
#include <core/types.hh>
#include <core/io/silent/SilentFileOptions.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <string>

namespace core {
namespace io {
namespace silent {

class SilentFileRMOptions : public basic::resource_manager::ResourceOptions
{
public:
	SilentFileRMOptions();

	SilentFileRMOptions( std::string const & name );

	~SilentFileRMOptions();

	SilentFileRMOptions( SilentFileRMOptions const & src );

	SilentFileOptions const &
	opts() const;

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
	type() const;

	static
	std::string
	class_name();

private:
	SilentFileOptionsOP options_;
};


} // namespace
} // namespace
} // namespace


#endif
