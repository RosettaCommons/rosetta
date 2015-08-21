// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


#ifndef INCLUDED_devel_dna_util_hh
#define INCLUDED_devel_dna_util_hh

#include <protocols/dna/typedefs.hh>


#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.fwd.hh>


#include <iosfwd>

#include <utility/vector1.hh>


namespace devel {
namespace dna {

void
make_base_pair_mutation(
	core::pose::Pose & pose,
	core::Size const seqpos,
	core::chemical::AA const & na
);

void
randomize_motif_sequence(
	core::pose::Pose & pose,
	core::Size const motif_begin,
	core::Size const motif_size
);

void
detect_interface_residues(
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & pos_list,
	core::Real const contact_threshold,
	utility::vector1< bool > & interface
);


void
detect_interface_by_nbrs(
	core::pose::Pose const & pose,
	utility::vector1< bool > & interface
);

void
detect_allatom_interface_residues(
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & pos_list,
	core::Real const contact_threshold,
	utility::vector1< bool > & interface
);

void
calc_protein_DNA_rmsd(
	core::pose::Pose const & pose,
	core::pose::Pose const & reference_pose,
	core::Real & ca_rmsd,
	core::Real & interface_ca_rmsd,
	core::Real & interface_allatom_rmsd
);


void
calc_DNA_bb_rmsd(
	core::pose::Pose const & pose,
	core::pose::Pose const & reference_pose,
	core::Real & dna_bb_rmsd
);


void
analyze_interface_sasa(
	core::pose::Pose const & complex_pose, // since we need to update nbr info
	int const split_jump,
	core::Real & bsasa14,
	core::Real & bsasa5,
	utility::vector1< core::id::AtomID > & buried_unsatisfied_donors,
	utility::vector1< core::id::AtomID > & buried_unsatisfied_acceptors
	//  int & buried_unsatisfied_donors,
	//  int & buried_unsatisfied_acceptors
);

void
analyze_interface_sasa(
	core::pose::Pose const & complex_pose,
	core::pose::Pose const & partner1,
	core::pose::Pose const & partner2,
	core::Real & bsasa14,
	core::Real & bsasa5,
	utility::vector1< core::id::AtomID > & buried_unsatisfied_donors,
	utility::vector1< core::id::AtomID > & buried_unsatisfied_acceptors
	//  int & buried_unsatisfied_donors,
	//  int & buried_unsatisfied_acceptors
);


/// uses default cutoff values, see .cc file
void
check_residue_proximity_to_dna(
	core::Size const ppos,
	protocols::dna::Positions const & dna_design_positions,
	core::pose::Pose const & pose,
	bool & close,  // output
	bool & contact // output
);


void
check_residue_proximity_to_dna(
	core::Size const ppos,
	protocols::dna::Positions const & dna_design_positions,
	core::pose::Pose const & pose,
	core::Real const close_threshold2,    // A**2
	core::Real const contact_threshold2,  // A**2
	core::Real const z_cutoff,            // A
	bool & close,  // output
	bool & contact // output
);

} // namespace dna
} // namespace devel

#endif
