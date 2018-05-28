// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/FixbbPwatSimAnnealer.cc
/// @brief  fixed temperature annealer for packing with statistical PWAT water model
/// @author Ryan Pavlovicz (rpavlov@uw.edu)

#ifndef INCLUDED_core_pack_annealer_FixbbPwatSimAnnealer_hh
#define INCLUDED_core_pack_annealer_FixbbPwatSimAnnealer_hh

// Unit Headers
#include <core/pack/annealer/FixbbPwatSimAnnealer.fwd.hh>

// Package Headers
#include <core/pack/annealer/RotamerAssigningAnnealer.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace annealer {

// store dwelltimes of water during packing in this struct
struct PointDwell {
	core::Vector xyz;
	core::Real dwell;
};


class FixbbPwatSimAnnealer : public RotamerAssigningAnnealer
{
public:
	typedef interaction_graph::AnnealableGraphBaseOP AnnealableGraphBaseOP;

public:
	FixbbPwatSimAnnealer(
		utility::vector0< int > & rot_to_pack,
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		AnnealableGraphBaseOP ig,
		FixbbRotamerSetsCOP rotamer_sets,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
	);

	virtual ~FixbbPwatSimAnnealer();
	void run();

	void
	set_min_dwell(core::Real mindwell) {
		min_dwell_ = mindwell;
	}

	utility::vector1< PointDwell > &
	get_dwell_times() {
		return all_rot_;
	}

private:
	AnnealableGraphBaseOP ig_;
	FixbbPwatSimAnnealer(const FixbbPwatSimAnnealer& rhs);

	utility::vector1< PointDwell > all_rot_;
	core::Real min_dwell_;
};

}//end namespace annealer
}//end namespace pack
}//end namespace core

#endif
