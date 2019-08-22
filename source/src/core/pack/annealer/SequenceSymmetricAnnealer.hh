// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/annealer/SequenceSymmetricAnnealer.hh
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_pack_annealer_SequenceSymmetricAnnealer_hh
#define INCLUDED_core_pack_annealer_SequenceSymmetricAnnealer_hh

// Unit Headers
#include <core/pack/annealer/SequenceSymmetricAnnealer.fwd.hh>

// Package Headers
#include <core/pack/annealer/RotamerAssigningAnnealer.hh>

#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>

#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace annealer {

/*
This annealer was created to perform the same purpose of link residues,
but hopefully in a more user-friendly way.

When used, the packer will enforce that all chains end up with the same sequence.
It uses pdb info to link residues together, so all residues with the same pdb number
will be the same amino acid in the end.

If a residue does not have a partner on every chain, it will not be allowed to mutate.

Like traditional symmetry, this assumes that all chains are part of the same symmetric system.
It is impossible to have, say, chains A+B+C where A+B are symmetric and C is separate.
*/

class SequenceSymmetricAnnealer : public RotamerAssigningAnnealer
{
public:
	typedef interaction_graph::AnnealableGraphBaseOP AnnealableGraphBaseOP;

public:
	SequenceSymmetricAnnealer(
		core::pose::Pose const & pose,
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

	SequenceSymmetricAnnealer(
		core::pose::Pose const & pose,
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		AnnealableGraphBaseOP ig,
		FixbbRotamerSetsCOP rotamer_sets,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
	);

	virtual ~SequenceSymmetricAnnealer();
	void run();

	void record_annealer_trajectory( bool setting );
	void trajectory_file_name( std::string const & setting );

private:
	std::string starting_sequence_;
	core::Size num_chains_;
	pose::PDBInfoCOP pdb_info_;
	AnnealableGraphBaseOP ig_;
	bool record_annealer_trajectory_;
	std::string trajectory_file_name_;
	SequenceSymmetricAnnealer( SequenceSymmetricAnnealer const & rhs );
};

}//end namespace annealer
}//end namespace pack
}//end namespace core

#endif
