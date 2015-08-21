// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/MultiCoolAnnealer.hh
/// @brief  Multiple low-temperature cooling cycles annealer class declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_annealer_MultiCoolAnnealer_hh
#define INCLUDED_core_pack_annealer_MultiCoolAnnealer_hh

/// Unit headers
#include <core/pack/annealer/MultiCoolAnnealer.fwd.hh>

/// Package headers
#include <core/pack/annealer/RotamerAssigningAnnealer.hh>

#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace annealer {


class MultiCoolAnnealer : public RotamerAssigningAnnealer
{
public:
	MultiCoolAnnealer(
		task::PackerTaskCOP task,
		utility::vector0<int> & rot_to_pack,
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		interaction_graph::InteractionGraphBaseOP ig,
		FixbbRotamerSetsCOP rotamer_sets,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
	);

	MultiCoolAnnealer(
		task::PackerTaskCOP task,
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		interaction_graph::InteractionGraphBaseOP ig,
		FixbbRotamerSetsCOP rotamer_sets,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
	);

	virtual ~MultiCoolAnnealer();
	void run();


private:

	void cool();
	void run_quench(
		ObjexxFCL::FArray1D_int & state_on_node,
		ObjexxFCL::FArray1D_int & best_state_on_node,
		core::PackerEnergy & best_energy,
		int num_cycles );

	void run_constant_temp_rotamer_substitutions(
		ObjexxFCL::FArray1D_int & state_on_node,
		ObjexxFCL::FArray1D_int & best_state_on_node,
		core::PackerEnergy & best_energy,
		int num_cycles
	);

	void store_top_energy(
		ObjexxFCL::FArray1D_int const & state_on_node,
		core::PackerEnergy energy );

	/// @brief unimplemented, private copy ctor -- uncopyable
	MultiCoolAnnealer( MultiCoolAnnealer const & rhs);
	/// @brief unimplemented, private assignment operator -- uncopyable
	MultiCoolAnnealer const & operator = ( MultiCoolAnnealer const & rhs );

private:

	static core::PackerEnergy const uninitialized_energy;

	interaction_graph::InteractionGraphBaseOP ig_;
	//std::vector<int> rot_to_pack_;

	//static int top_to_keep_static;

	ObjexxFCL::FArray1D_int nsteps_for_rot_;
	int nsteps_;
	Size top_to_keep;
	ObjexxFCL::FArray2D_int top_netstates_;
	ObjexxFCL::FArray1D< core::PackerEnergy > energy_top_;
	core::PackerEnergy worst_top_energy_;
	int which_netstate_worst_top_;
	int num_top_kept_;

};

}// namespace annealer
}// namespace pack
}// namespace core

#endif
