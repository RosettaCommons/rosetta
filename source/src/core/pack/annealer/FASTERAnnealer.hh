// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/annealer/FASTERAnnealer.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_annealer_FASTERAnnealer_hh
#define INCLUDED_core_pack_annealer_FASTERAnnealer_hh

// Unit headers
#include <core/pack/annealer/FASTERAnnealer.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/pack/annealer/RotamerAssigningAnnealer.hh>
#include <core/pack/interaction_graph/FASTERInteractionGraph.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>

// Utility headers
#include <utility/vector0.hh>

// Objexx Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace annealer {

class FASTERAnnealer : public RotamerAssigningAnnealer
{
public:
	typedef RotamerAssigningAnnealer parent;
	typedef rotamer_set::FixbbRotamerSetsCOP FixbbRotamerSetsCOP;

public:
	FASTERAnnealer(
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		interaction_graph::FASTERInteractionGraphOP ig,
		FixbbRotamerSetsCOP rotamer_sets,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
	);

	virtual ~FASTERAnnealer();

	void run();

private:
	void iBR();
	void trySeveral_ciBRs();
	void ciBR();
	void sBR();
	void dBR();
	void finalize_output();

	void run_quench_cycles();

	void reset_recent_network_state_history();
	void note_current_network_state( ObjexxFCL::FArray1_int const & netstate );
	bool stuck_in_network_state_loop();
	int pick_rotamer_for_node( int node );
	int pick_a_rotamer_for_sBR();
	void shuffle_sBR_rotamers();

	int hash_recent_history( int history_index );

public:
	void set_ciBR_only( bool setting );
	void set_num_sa_trajectories( int setting );
	void set_sa_length_scale( Real setting );
	//void set_sBR_limit( int setting ); // -1 for "run until completion"
private:

	interaction_graph::FASTERInteractionGraphOP ig_;
	FixbbRotamerSetsCOP rotamer_sets_;
	int const num_nodes_;

	ObjexxFCL::FArray2D_int recent_network_state_history_;
	ObjexxFCL::FArray1D_int recent_history_hash_values_;
	ObjexxFCL::FArray1D_int recent_history_hash_count_;
	int recent_history_head_;
	int curr_in_recent_history_;
	bool netstate_duplicated_;

	static int const recent_history_size_ = 100;
	static int const hash_size_ = 2017; //prime

	int progress_through_sBR_;
	utility::vector0< int > sBR_rotamers_;

	bool ciBR_only_;
	int num_sa_trajectories_; // default of 5
	Real sa_inner_iterations_length_scale_; // default of 0.05
	int sBR_limit_;

};

}
}
}

#endif
