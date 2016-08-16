// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/annealer/FixbbLinkingRotamerSimAnnealer.hh
/// @brief  Packer's standard annealer class declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_annealer_FixbbLinkingRotamerSimAnnealer_hh
#define INCLUDED_core_pack_annealer_FixbbLinkingRotamerSimAnnealer_hh

// Unit Headers
#include <core/pack/annealer/FixbbLinkingRotamerSimAnnealer.fwd.hh>

// Package Headers
#include <core/pack/annealer/RotamerAssigningAnnealer.hh>

#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>

#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerLinks.fwd.hh>

namespace core {
namespace pack {
namespace annealer {

class FixbbLinkingRotamerSimAnnealer;

class FixbbLinkingRotamerSimAnnealer : public RotamerAssigningAnnealer
{
public:
	typedef interaction_graph::AnnealableGraphBaseOP AnnealableGraphBaseOP;
	typedef rotamer_set::RotamerSetsCOP RotamerSetsCOP;
	typedef rotamer_set::RotamerSetCOP RotamerSetCOP;
	typedef rotamer_set::RotamerLinksOP   RotamerLinksOP;
	typedef rotamer_set::RotamerLinksCOP  RotamerLinksCOP;

public:
	FixbbLinkingRotamerSimAnnealer(
		utility::vector0<int> & rot_to_pack,
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		AnnealableGraphBaseOP ig,
		FixbbRotamerSetsCOP rotamer_sets,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq,
		RotamerLinksCOP rotamer_links
	);

	FixbbLinkingRotamerSimAnnealer(
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		AnnealableGraphBaseOP ig,
		FixbbRotamerSetsCOP rotamer_sets,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq,
		RotamerLinksCOP rotamer_links
	);


	virtual ~FixbbLinkingRotamerSimAnnealer();

	/// @brief sim_annealing for fixed backbone design mode
	void run();

private:
	void
	setup_rotamer_links( RotamerLinksCOP rotamer_links );

	AnnealableGraphBaseOP ig_;
	FixbbLinkingRotamerSimAnnealer(const FixbbLinkingRotamerSimAnnealer& rhs);
	RotamerLinksOP rotamer_links_;
};

}//end namespace annealer
}//end namespace pack
}//end namespace core

#endif
