// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/RsdEnergiesBreakdown.fwd.hh
///
///	@brief		Residue Energies Breakdown
///	@details	Returns individual one body, two-body and whole strucutre energies after a pose is
///				sored. Useful for unit testing. Rocco wrote an awesome app for this in
///				public/analysis/residue_energy_breakdown.cc. I just made it into a class.
///				Last Modified: 4/24/14
///
///	@author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_RsdEnergiesBreakdown_fwd_hh
#define INCLUDED_core_scoring_RsdEnergiesBreakdown_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh> 

namespace core {
namespace scoring {

class RsdEnergiesBreakdown;
typedef utility::pointer::owning_ptr< RsdEnergiesBreakdown > RsdEnergiesBreakdownOP;
typedef utility::pointer::owning_ptr< RsdEnergiesBreakdown const > RsdEnergiesBreakdownCOP;

} // scoring
} // core

#endif // ICNLUDED_core_scoring_RsdEnergiesBreakdown_fwd_hh
