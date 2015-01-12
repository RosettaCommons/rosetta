// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffractionOptions.hh
/// @brief  Options for fiber diffraction data
/// @author Wojciech Potrzebowski and Ingemar Andre

#ifndef INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionOptions_hh
#define INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionOptions_hh


// Unit Headers
#include <core/scoring/fiber_diffraction/FiberDiffractionOptions.fwd.hh>
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
namespace fiber_diffraction {

class FiberDiffractionOptions : public basic::resource_manager::ResourceOptions
{
public:
	FiberDiffractionOptions();

	FiberDiffractionOptions(
		std::string const & name);

	FiberDiffractionOptions(
		std::string const & name,
		Real c,
		Real res_cutoff_high,
  	Real res_cutoff_low );

	FiberDiffractionOptions(
		FiberDiffractionOptions const & src);

	~FiberDiffractionOptions();

	Real get_res_high() const;

  void set_res_high( Real res_cutoff_high_ );

  Real get_res_low() const;

  void set_res_low( Real res_cutoff_low_ );

	Real get_c_repeat() const;

  void set_c_repeat( Real c_ );


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
	type() const { return "FiberDiffractionOptions"; }

private:
	Real c;
	Real res_cutoff_high;
	Real res_cutoff_low;

};


} // namespace
} // namespace
} // namespace


#endif
