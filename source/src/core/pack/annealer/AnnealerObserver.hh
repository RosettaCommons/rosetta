// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/annealer/AnnealerObserver.hh
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_pack_annealer_AnnealerObserver_hh
#define INCLUDED_core_pack_annealer_AnnealerObserver_hh

/// Unit headers
#include <core/pack/annealer/AnnealerObserver.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>

/// Package headers

/// C++ headers


#include <utility/VirtualBase.hh>

#include <core/conformation/Residue.fwd.hh> // AUTO IWYU For Residue

namespace core {
namespace pack {
namespace annealer {

class AnnealerObserver : public utility::VirtualBase
{
public:

	///@brief Do whatever you need to do before monte carlo starts. This will be called exactly once before observe_mc
	virtual
	void setup_for_mc(
		rotamer_set::FixbbRotamerSetsCOP const & rotamer_sets
	) = 0;

	///@brief This will be called every time a Monte Carlo move is considered
	virtual
	void observe_mc(
		int resid, //Residue ID
		int mres,  //Molten Residue ID
		core::conformation::Residue const & existing_rotamer,
		core::conformation::Residue const & candidate_rotamer,
		float temperature, //Monte Carlo Temperature
		bool passed //True if the Monte Carlo move was accepted
	) = 0;

};

}// namespace annealer
}// namespace pack
}// namespace core

#endif
