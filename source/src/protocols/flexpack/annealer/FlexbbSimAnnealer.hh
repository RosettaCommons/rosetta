// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/annealer/FlexbbSimAnnealer
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_flexpack_annealer_FlexbbSimAnnealer_hh
#define INCLUDED_protocols_flexpack_annealer_FlexbbSimAnnealer_hh

/// Unit headers
#include <protocols/flexpack/annealer/FlexbbSimAnnealer.fwd.hh>

/// Package headers
#include <protocols/flexpack/interaction_graph/FlexbbInteractionGraph.fwd.hh>
#include <protocols/flexpack/rotamer_set/FlexbbRotamerSets.fwd.hh>

/// Project headers
#include <core/pack/annealer/SimAnnealerBase.hh>

#include <utility/vector1.fwd.hh>


namespace protocols {
namespace flexpack {
namespace annealer {

class FlexbbSimAnnealer : public core::pack::annealer::SimAnnealerBase
{
public:
	typedef core::pack::annealer::SimAnnealerBase parent;
	typedef core::Size Size;
	typedef core::PackerEnergy PackerEnergy;

public:
	FlexbbSimAnnealer(
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		interaction_graph::FlexbbInteractionGraphOP ig,
		rotamer_set::FlexbbRotamerSetsCOP rotsets,
		ObjexxFCL::FArray1D_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< PackerEnergy > & rot_freq
	);

	virtual ~FlexbbSimAnnealer();

	void run();

protected:

	Size
	pick_a_rotamer(
		Size outercycle,
		Size innercycle,
		Size cycle_number,
		utility::vector1< Size > & accessible_state_list
	) const;

	bool
	pass_metropolis_multiple_nodes_changing(
		PackerEnergy previous_fragmentE,
		PackerEnergy deltaE,
		Size num_changing_nodes
	) const;

private:

	interaction_graph::FlexbbInteractionGraphOP ig_;
	rotamer_set::FlexbbRotamerSetsCOP rotsets_;

	FlexbbSimAnnealer(const FlexbbSimAnnealer & rhs);

};


}
}
}

#endif
