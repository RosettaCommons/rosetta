// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/align/StepWisePoseAligner.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_align_StepWisePoseAligner_HH
#define INCLUDED_protocols_stepwise_modeler_align_StepWisePoseAligner_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/modeler/align/StepWisePoseAligner.fwd.hh>
#include <core/pose/Pose.fwd.hh>
// Need full header for map hash
#include <core/id/AtomID.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <map>
#include <string>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace align {

class StepWisePoseAligner: public utility::pointer::ReferenceCount {

public:

	//constructor
	StepWisePoseAligner( core::pose::Pose const & reference_pose );

	//destructor
	~StepWisePoseAligner();

public:

	// could probably make this a mover?
	void
	apply( core::pose::Pose & pose );

	core::Real
	get_rmsd_over_all_poses( core::pose::Pose & pose );

	void
	initialize( core::pose::Pose const & pose );

	core::Real
	get_rmsd_no_superimpose( core::pose::Pose const & pose,
		bool const check_align = true );

	core::Real
	get_rmsd_no_superimpose( core::pose::Pose const & pose,
		core::pose::Pose const & pose_reference,
		bool const check_align = true );

	void
	create_coordinate_constraints( core::pose::Pose & pose,
		core::Real const rmsd_screen );

	std::map< core::id::AtomID, core::id::AtomID>
	create_coordinate_constraint_atom_id_map( core::pose::Pose const & pose );

	void set_superimpose_over_all_instantiated( bool const & setting ){ superimpose_over_all_instantiated_ = setting; }
	bool superimpose_over_all_instantiated() const{ return superimpose_over_all_instantiated_; }

	void set_root_partition_res( utility::vector1< core::Size > const & setting ){ root_partition_res_ = setting; }
	void set_user_defined_calc_rms_res( utility::vector1< core::Size > const & setting ){ user_defined_calc_rms_res_ = setting; }

	core::Real rmsd() const { return rmsd_;}
	core::Size natoms_rmsd() const { return calc_rms_atom_id_map_.size();}

	core::Real superimpose_rmsd() const { return superimpose_rmsd_;}
	core::Size natoms_superimpose_rmsd() const { return superimpose_atom_id_map_.size();}

	std::map< core::id::AtomID, core::id::AtomID > const & superimpose_atom_id_map(){ return superimpose_atom_id_map_; }

	core::Real rmsd_over_alignment_atoms() const { return superimpose_rmsd_;}

	bool check_matching_atom_names( core::pose::Pose const & pose1, core::pose::Pose const & pose2,
		bool const verbose = true );

	bool check_matching_atom_names( core::pose::Pose const & pose1, core::pose::Pose const & pose2,
		std::map< core::id::AtomID, core::id::AtomID > const & atom_id_map, bool const verbose = true );

	utility::vector1< core::Size > const & rmsd_res_in_pose(){ return rmsd_res_in_pose_;}
	utility::vector1< core::Size > const & superimpose_res_in_pose(){ return superimpose_res_in_pose_;}

private:

	void
	superimpose_recursively( core::pose::Pose & pose,
		core::Real & rmsd_over_all, core::Size & natoms_over_all );

	core::Real
	do_superimposition( core::pose::Pose & pose );

	void
	update_reference_pose_local( core::pose::Pose const & pose );

	utility::vector1< core::Size >
	get_res_list_in_reference( core::pose::Pose const & pose ) const;

	void
	update_calc_rms_atom_id_map( core::pose::Pose const & pose );

	void
	get_calc_rms_atom_id_map( std::map< core::id::AtomID, core::id::AtomID > & calc_rms_atom_id_map,
		core::pose::Pose const & pose,
		utility::vector1 < core::Size > const & calc_rms_res ) const;

	void
	update_superimpose_atom_id_map( core::pose::Pose const & pose );

	core::Size
	get_rmsd_res_and_superimpose_res_in_pose( core::pose::Pose const & pose );

	void
	add_coordinate_constraints_from_map( core::pose::Pose & pose, core::pose::Pose const & reference_pose,
		std::map< core::id::AtomID, core::id::AtomID > const & atom_id_map,
		core::Real const & constraint_x0, core::Real const & constraint_tol ) const;

	bool
	do_checks( std::string const & atom_name, core::Size const n, core::pose::Pose const & pose ) const;

	bool
	add_to_atom_id_map_after_checks( std::map< core::id::AtomID, core::id::AtomID> & atom_id_map,
		std::string const & atom_name,
		core::Size const n1, core::Size const n2,
		core::pose::Pose const & pose1, core::pose::Pose const & pose2,
		bool const do_the_checks = true ) const;

	std::map< core::id::AtomID, core::id::AtomID >
	get_root_triad_atom_id_map( core::pose::Pose const & pose ) const;

	void
	output_atom_id_map( std::map< core::id::AtomID, core::id::AtomID > const & atom_id_map ) const;

	void
	output_atom_id_map( std::map< core::id::AtomID, core::id::AtomID > const & atom_id_map,
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2 ) const;

	std::map< core::id::AtomID, core::id::AtomID > filter_virtual_atoms( std::map< core::id::AtomID, core::id::AtomID > const & atom_id_map,
		core::pose::Pose const & pose );

private:

	core::pose::Pose const & reference_pose_;
	core::pose::PoseCOP reference_pose_local_;
	core::pose::PoseOP mod_reference_pose_local_ = nullptr;
	utility::vector1< core::Size > root_partition_res_;
	utility::vector1< core::Size > user_defined_calc_rms_res_;

	utility::vector1< core::Size > rmsd_res_in_pose_;
	utility::vector1< core::Size > superimpose_res_in_pose_;
	std::map< core::id::AtomID, core::id::AtomID > calc_rms_atom_id_map_;
	std::map< core::id::AtomID, core::id::AtomID > complete_moving_atom_id_map_;
	std::map< core::id::AtomID, core::id::AtomID > superimpose_atom_id_map_;
	std::string annotated_sequence_used_for_atom_id_maps_;

	core::Real rmsd_;
	core::Real superimpose_rmsd_;

	core::Real const check_alignment_tolerance_;
	bool superimpose_over_all_instantiated_;

};

} //align
} //modeler
} //stepwise
} //protocols

#endif
