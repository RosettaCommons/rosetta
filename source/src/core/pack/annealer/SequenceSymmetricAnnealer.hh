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
/// @author Updated by Tim Neary, timdot10@gmail.com


#ifndef INCLUDED_core_pack_annealer_SequenceSymmetricAnnealer_hh
#define INCLUDED_core_pack_annealer_SequenceSymmetricAnnealer_hh

// Unit Headers
#include <core/pack/annealer/SequenceSymmetricAnnealer.fwd.hh>

// Package Headers
#include <core/pack/annealer/RotamerAssigningAnnealer.hh>

#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/cxx_versioning_macros.hh>

// External Headers
#include <boost/function.hpp>
#include <boost/bind/bind.hpp>

// C++ headers
#include <string>
#include <unordered_set>
#include <unordered_map>
// #include <map> // not needed unless using the callback map

namespace core {
namespace pack {
namespace annealer {

typedef utility::pointer::shared_ptr< std::unordered_set< char > > CharSetOP;

struct SeqSymmAnnealerSetup {
	utility::vector1< utility::vector1< core::Size > > const corresponding_mress_for_mres;
	utility::vector1< CharSetOP > const common_res_types;
};

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
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq,
		std::string const & rs_prefix
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
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq,
		std::string const & rs_prefix
	);

	~SequenceSymmetricAnnealer() override;
	void run() override;

	void record_annealer_trajectory( bool setting );
	void trajectory_file_name( std::string const & setting );

	NODISCARD_ATTR
	utility::vector1< utility::vector1< Size > >
	create_corresponding_mress_for_mres() const;

private: // Member methods

	SequenceSymmetricAnnealer( SequenceSymmetricAnnealer const & rhs );

	/// @brief Determines whether the selection logic passed to the annealer is valid, must be equal to "type_only", "identity", or "try_identity".
	// void sele_logic_valid( std::string const & rot_sele_logic ) const;

	/// @brief Uses the rs_prefix_ to search the pose for stored residue subsets.
	/// It is expected that all relevant residue subsets will be prefixed by the rs_prefix_ and be in the form:
	/// <rs_prefix_>_<linked_residue_subset_id>_<residue_subset_number_in_group>
	void search_pose_for_residue_subsets( core::pose::Pose const & pose );

	/// @brief Prints the linked residues in a useful, human readable format.
	/// Uses the Rosetta number for linked residues.
	void print_linked_residues( utility::vector1< utility::vector1< core::Size > > const & corresponding_mress_for_mres ) const;

	/// @brief Updates map with all residue types for moltenresidue, id.
	/// Map contains a count for each type of residue previously seen
	void
	update_shared_residue_map(
		core::Size const id,
		std::unordered_map< char, core::Size > & map
	) const;

	/// @brief Returns the common set of residues based on a vector of linked moltenred ids.
	std::unordered_set< char >
	get_shared_residue_types(
		core::Size const curr_mres,
		utility::vector1< core::Size > const & linked_res
	) const;

	/// @brief Sets up the necessary objects required for running the annealer.
	/// Returns vector of vector of mres for linked residues and common residue types for each linked res set.
	SeqSymmAnnealerSetup
	setup_for_linked_residues();

private: // Member variables

	std::string starting_sequence_;
	pose::PDBInfoCOP pdb_info_;
	AnnealableGraphBaseOP ig_;
	bool record_annealer_trajectory_;
	std::string trajectory_file_name_;

	//Documented on wiki site for SetupForSequenceSymmetryMover
	bool power_mode_ = false;
	utility::vector1< utility::vector1< utility::vector1< bool > > > residue_subsets_;
	std::string rs_prefix_ = "SequenceSymmetricAnnealer__";

};

}//end namespace annealer
}//end namespace pack
}//end namespace core

#endif
