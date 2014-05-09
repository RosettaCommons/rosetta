// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange ( oliver.lange@tum.de )

#ifndef INCLUDED_protocols_canonical_sampling_WTEBiasEnergy_hh
#define INCLUDED_protocols_canonical_sampling_WTEBiasEnergy_hh

// Unit Headers
//#include <protocols/canonical_sampling/WTEBiasEnergy.fwd.hh>
#include <protocols/canonical_sampling/BiasEnergy.hh>

#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <core/types.hh>

namespace protocols {
namespace canonical_sampling {

///@details
class WTEBiasEnergy : public BiasEnergy {
	typedef BiasEnergy Parent;
public:
	WTEBiasEnergy( core::Size stride, core::Real omega, core::Real gamma );
	WTEBiasEnergy();

	std::string
	get_name() const { return "WTEBiasEnergy"; };

protected:
	virtual core::Real extract_collective_var( core::pose::Pose const& ) const;
};

} //namespace canonical_sampling
} //namespace protocols

#endif //INCLUDED_protocols_canonical_sampling_WTEBiasEnergy_HH
