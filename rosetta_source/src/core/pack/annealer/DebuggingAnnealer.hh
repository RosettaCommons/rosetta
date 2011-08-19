// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/DebuggingAnnealer.hh
/// @brief  Debugging annealer class declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_pack_annealer_DebuggingAnnealer_hh
#define INCLUDED_core_pack_annealer_DebuggingAnnealer_hh

// Unit Headers
#include <core/pack/annealer/DebuggingAnnealer.fwd.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/annealer/RotamerAssigningAnnealer.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>

// ObjexxFCL Headers

namespace core {
namespace pack {
namespace annealer {

class DebuggingAnnealer : public RotamerAssigningAnnealer
{
 public:
	DebuggingAnnealer(
		utility::vector0< int > & rot_to_pack,
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		float & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		interaction_graph::InteractionGraphBaseOP ig,
		const rotamer_set::RotamerSet * p_rotamer_set,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D_float & rot_freq
	);

	DebuggingAnnealer(
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		float & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		interaction_graph::InteractionGraphBaseOP ig,
		rotamer_set::RotamerSetCOP p_rotamer_set,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D_float & rot_freq
	);

	virtual ~DebuggingAnnealer();
	void run();

 private:
	interaction_graph::InteractionGraphBaseOP ig_;
	DebuggingAnnealer(const DebuggingAnnealer& rhs);
};

}//end namespace annealer
}//end namespace pack
}//end namespace core

#endif
