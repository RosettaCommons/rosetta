// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Align a random jump to template
/// @detailed
/// @author Yifan Song

#ifndef INCLUDED_protocols_moves_InsertChunkMover_hh
#define INCLUDED_protocols_moves_InsertChunkMover_hh

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <numeric/xyz.functions.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace comparative_modeling {
namespace hybridize {

using namespace core;

class InsertChunkMover: public protocols::moves::Mover
{

public:
	
InsertChunkMover();
~InsertChunkMover();
	
// atom_map: from mod_pose to ref_pose
void
get_superposition_transformation(
								 pose::Pose const & mod_pose,
								 pose::Pose const & ref_pose,
								 id::AtomID_Map< id::AtomID > const & atom_map,
								 numeric::xyzMatrix< core::Real > &R, numeric::xyzVector< core::Real > &preT, numeric::xyzVector< core::Real > &postT );

void
apply_transform(
				pose::Pose & mod_pose,
				std::list <Size> const & residue_list,
				numeric::xyzMatrix< core::Real > const & R, numeric::xyzVector< core::Real > const & preT, numeric::xyzVector< core::Real > const & postT
				);
	
void align_chunk(core::pose::Pose & pose);

void set_bb_xyz_aligned(core::pose::Pose & pose);

bool success() {
    return success_;
}
    
void set_template(core::pose::PoseCOP template_pose,
				  std::map <core::Size, core::Size> const & sequence_alignment );

void set_aligned_chunk(core::pose::Pose const & pose, Size const jump_number);

void set_reset_torsion_unaligned(bool reset_torsion_unaligned);
	
void steal_torsion_from_template(core::pose::Pose & pose);
	
void steal_torsion_and_bonds_from_template(core::pose::Pose & pose);
	
bool get_local_sequence_mapping(core::pose::Pose const & pose,
								int registry_shift = 0,
								Size MAX_TRIAL = 10 );
	
void set_registry_shift(int registry_shift);
	
Size trial_counter(Size ires);
	
void apply(core::pose::Pose & pose);
	
std::string get_name() const;
	
private:
	core::pose::PoseCOP template_pose_;
	std::map <core::Size, core::Size> sequence_alignment_;

	Size jump_number_; // the jump to be realigned
	Size seqpos_start_; // start and end seqpose of the chunk
	Size seqpos_stop_;
    bool success_;
	int registry_shift_;
	char secstruct_;

	// parameters of the protocol
	bool reset_torsion_unaligned_; // reset torsion of unaligned region to the default value for the given secondary structure
	bool align_to_ss_only_; // only use the secondary structure portion to align to the template
	bool copy_ss_torsion_only_; // only copy the secondary structure information from the template
	
	std::map <core::Size, core::Size> sequence_alignment_local_; // with registry shift of the aligned chunk
	core::id::AtomID_Map< core::id::AtomID > atom_map_; // atom map for superposition
	utility::vector1 <Size> align_trial_counter_;
};
	
} // hybridize 
} // comparative_modeling 
} // protocols

#endif
