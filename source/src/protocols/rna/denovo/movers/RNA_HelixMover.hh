// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RNA_HelixMover.hh
/// @brief
/// @details
///
/// @author Kalli Kappel


#ifndef INCLUDED_protocols_rna_RNA_HelixMover_hh
#define INCLUDED_protocols_rna_RNA_HelixMover_hh

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>

//// C++ headers
#include <string>
#include <vector>

#include <protocols/rna/denovo/movers/RNA_HelixMover.fwd.hh>
#include <protocols/rna/denovo/base_pairs/RNA_BasePairHandler.fwd.hh>
#include <core/kinematics/FoldTree.hh>

#include <utility/vector1.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>


namespace protocols {
namespace rna {
namespace denovo {
namespace movers {

/// @brief rotations and translations of helix along helical axis
class RNA_HelixMover: public protocols::moves::Mover {
public:
	/// @brief Construct the mover object
	RNA_HelixMover( utility::vector1< utility::vector1< Size > > helix_regions,
		protocols::rna::denovo::base_pairs::RNA_BasePairHandlerCOP rna_base_pair_handler,
		bool const & move_first_rigid_body );

	~RNA_HelixMover();

	void apply( core::pose::Pose & pose );

	void set_pose( core::pose::Pose const & pose );

	void get_helix_ends();

	virtual std::string get_name() const;

	void set_rot_magnitude( core::Real const & rot_mag ) { rot_mag_ = rot_mag; }

	void set_trans_magnitude( core::Real const & trans_mag ) { trans_mag_ = trans_mag; }

private:
	//private methods
	std::pair< core::Vector, core::Vector >
	get_helical_axis_and_center( core::pose::Pose const & pose,
		core::Size const & region ) const;

	core::Vector
	get_bp_center( core::pose::Pose const & pose,
		std::pair< core::Size, core::Size > const & bp_res ) const;

	core::Vector
	get_bb_pos( core::pose::Pose const & pose,
		core::Size const & res ) const;

	core::Vector
	get_backbone_centroid( core::pose::Pose const & pose,
		core::Size const & res ) const;

	//data
	utility::vector1< utility::vector1< Size > > helix_regions_;
	utility::vector1< utility::vector1< Size > > helix_regions_with_jumps_and_ends_;
	utility::vector1< Size > helix_regions_jumps_;
	utility::vector1< Size > helix_regions_jumps_final_;
	protocols::rna::denovo::base_pairs::RNA_BasePairHandlerCOP rna_base_pair_handler_;
	utility::vector1< std::pair< std::pair< Size, Size>, std::pair< Size, Size > > > helix_ends_;
	utility::vector1< std::pair< std::pair< Size, Size>, std::pair< Size, Size > > > helix_ends_final_;
	core::kinematics::FoldTree pose_fold_tree_;
	bool pose_is_set_;
	bool move_first_rigid_body_;
	utility::vector1< protocols::rigid::RigidBodySpinMoverOP > spin_movers_;
	utility::vector1< protocols::rigid::RigidBodyTransMoverOP > trans_movers_;
	core::Real rot_mag_;
	core::Real trans_mag_;

}; // class RNA_HelixMover


} //movers
} //denovo
} //rna
} //protocols

#endif
