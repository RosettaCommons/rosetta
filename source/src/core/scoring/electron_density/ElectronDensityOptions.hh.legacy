// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/electron_density/ElectronDensityOptions.hh
/// @brief  Options for constructing an electron density map
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_core_scoring_electron_density_ElectronDensityOptions_hh
#define INCLUDED_core_scoring_electron_density_ElectronDensityOptions_hh


// Unit Headers
#include <core/scoring/electron_density/ElectronDensityOptions.fwd.hh>
#include <basic/resource_manager/ResourceOptions.hh>

// Platform Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <string>

namespace core {
namespace scoring {
namespace electron_density {

class ElectronDensityOptions : public basic::resource_manager::ResourceOptions
{
public:
	ElectronDensityOptions();

	ElectronDensityOptions(
		std::string const & name);

	ElectronDensityOptions(
		std::string const & name,
		Real mapreso,
		Real grid_spacing);

	~ElectronDensityOptions();

	ElectronDensityOptions(
		ElectronDensityOptions const & src);

	Real
	get_mapreso() const;

	void
	set_mapreso( Real mapreso );

	Real
	get_grid_spacing() const;

	void
	set_grid_spacing( Real grid_spacing);

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
	type() const { return "ElectronDensityOptions"; }

private:
	Real mapreso_;
	Real grid_spacing_;

};


} // namespace
} // namespace
} // namespace


#endif
