// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/annealer/SmartFixbbSimAnnealer.hh
/// @brief  Derivation of the FixbbSimAnnealer with some fancy Tensorflow logic to adaptively decrease sample space during the run
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_pack_annealer_SmartFixbbSimAnnealer_hh
#define INCLUDED_core_pack_annealer_SmartFixbbSimAnnealer_hh

// Unit Headers
#include <core/pack/annealer/SmartFixbbSimAnnealer.fwd.hh>

// Package Headers
#include <core/pack/annealer/RotamerAssigningAnnealer.hh>

#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace annealer {

class SmartFixbbSimAnnealer;

class SmartFixbbSimAnnealer : public RotamerAssigningAnnealer
{
public:
	typedef interaction_graph::AnnealableGraphBaseOP AnnealableGraphBaseOP;

public:
	SmartFixbbSimAnnealer(
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

	SmartFixbbSimAnnealer(
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		AnnealableGraphBaseOP ig,
		FixbbRotamerSetsCOP rotamer_sets,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
	);

	~SmartFixbbSimAnnealer() override;
	void run() override;

	void record_annealer_trajectory( bool setting );
	void trajectory_file_name( std::string const & setting );

	void perform_validation_test() const;

	void initialize_from_task( core::pack::task::PackerTask const & task );

private:
	AnnealableGraphBaseOP ig_;
	bool record_annealer_trajectory_;
	std::string trajectory_file_name_;
	SmartFixbbSimAnnealer(const SmartFixbbSimAnnealer& rhs);

	bool has_been_initialized_from_task_ = false;
	std::string model_ = "generation2";
	core::Real cutoff_ = 0.25;
	bool pick_again_ = true;
	bool disable_during_quench_ = true;

};

}//end namespace annealer
}//end namespace pack
}//end namespace core

#endif
