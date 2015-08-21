// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/FixbbSimAnnealer.hh
/// @brief  Packer's standard annealer class declaration, originally written by Brian Kuhlman
/// and factored into base and derived classes by Jenny Hu.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_annealer_FixbbSimAnnealer_hh
#define INCLUDED_core_pack_annealer_FixbbSimAnnealer_hh

// Unit Headers
#include <core/pack/annealer/FixbbSimAnnealer.fwd.hh>

// Package Headers
#include <core/pack/annealer/RotamerAssigningAnnealer.hh>

#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>

#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace annealer {

class FixbbSimAnnealer;

class FixbbSimAnnealer : public RotamerAssigningAnnealer
{
public:
	typedef interaction_graph::InteractionGraphBaseOP InteractionGraphBaseOP;

public:
	FixbbSimAnnealer(
		utility::vector0< int > & rot_to_pack,
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		InteractionGraphBaseOP ig,
		FixbbRotamerSetsCOP rotamer_sets,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
	);

	FixbbSimAnnealer(
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		InteractionGraphBaseOP ig,
		FixbbRotamerSetsCOP rotamer_sets,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
	);

	virtual ~FixbbSimAnnealer();
	void run();

	void record_annealer_trajectory( bool setting );
	void trajectory_file_name( std::string const & setting );

private:
	InteractionGraphBaseOP ig_;
	bool record_annealer_trajectory_;
	std::string trajectory_file_name_;
	FixbbSimAnnealer(const FixbbSimAnnealer& rhs);
};

}//end namespace annealer
}//end namespace pack
}//end namespace core

#endif
