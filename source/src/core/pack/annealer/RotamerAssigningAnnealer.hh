// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/RotamerAssigningAnnealer.hh
/// @brief  Residue assigning annealer class declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_annealer_RotamerAssigningAnnealer_hh
#define INCLUDED_core_pack_annealer_RotamerAssigningAnnealer_hh

// Unit Headers
#include <core/pack/annealer/RotamerAssigningAnnealer.fwd.hh>

// Package Headers
#include <core/pack/annealer/SimAnnealerBase.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>


#include <utility/vector0.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace annealer {

class RotamerAssigningAnnealer;

class RotamerAssigningAnnealer : public SimAnnealerBase
{
public:
	typedef rotamer_set::FixbbRotamerSetsCOP FixbbRotamerSetsCOP;

public:
	RotamerAssigningAnnealer(
		int num_of_rot_to_pack,
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		FixbbRotamerSetsCOP p_rotamer_set,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
	);

	RotamerAssigningAnnealer(
		utility::vector0< int > & rot_to_pack,
		int num_of_rot_to_pack,
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		FixbbRotamerSetsCOP p_rotamer_set,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
	);

	virtual ~RotamerAssigningAnnealer();

	int pick_a_rotamer( int cycle );
	int pick_a_rotamer_for_node( int node ) const;

	void set_assign_state_to_all_nodes_immediately( bool setting );

protected:
	FixbbRotamerSetsCOP rotamer_sets() const;
	utility::vector0< int > const & rot_to_pack() const;

private:

	void setup_rots_for_node(
		FixbbRotamerSetsCOP rotamer_sets
	);

private:
	FixbbRotamerSetsCOP rotamer_sets_;
	utility::vector0< int > rot_to_pack_;
	utility::vector1< utility::vector1< int > > rots_for_nodes_;
	int current_to_pick_;
	int  n_assigned_at_start_;
	bool assign_state_to_all_nodes_immediately_;

};

} // ene namespace annealer
} // end namespace pack
} // end namespace core

#endif
