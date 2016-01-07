// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/align/StepWisePoseAligner.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_align_StepWisePoseAligner_HH
#define INCLUDED_protocols_stepwise_modeler_align_StepWisePoseAligner_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/modeler/align/StepWisePoseAligner.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <map>
#include <string>

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace align {

class StepWisePoseAligner: public utility::pointer::ReferenceCount {

public:

	//constructor
	StepWisePoseAligner( pose::Pose const & reference_pose );

	//destructor
	~StepWisePoseAligner();

public:

	// could probably make this a mover?
	void
	apply( pose::Pose & pose );

	Real
	get_rmsd_over_all_poses( pose::Pose & pose );

	void
	initialize( pose::Pose const & pose );

	Real
	get_rmsd_no_superimpose( pose::Pose const & pose,
		bool const check_align  = true );

	Real
	get_rmsd_no_superimpose( pose::Pose const & pose,
		pose::Pose const & pose_reference,
		bool const check_align  = true );

	void
	create_coordinate_constraints( pose::Pose & pose,
		Real const rmsd_screen );

	void set_superimpose_over_all_instantiated( bool const & setting ){ superimpose_over_all_instantiated_ = setting; }
	bool superimpose_over_all_instantiated() const{ return superimpose_over_all_instantiated_; }

	void set_root_partition_res( utility::vector1< Size > const & setting ){ root_partition_res_ = setting; }
	void set_user_defined_calc_rms_res( utility::vector1< Size > const & setting ){ user_defined_calc_rms_res_ = setting; }

	Real rmsd() const { return rmsd_;}
	Size natoms_rmsd() const { return calc_rms_atom_id_map_.size();}

	Real superimpose_rmsd() const { return superimpose_rmsd_;}
	Size natoms_superimpose_rmsd() const { return superimpose_atom_id_map_.size();}

	std::map< id::AtomID, id::AtomID > const & superimpose_atom_id_map(){ return superimpose_atom_id_map_; }

	Real rmsd_over_alignment_atoms() const { return superimpose_rmsd_;}

	bool check_matching_atom_names( pose::Pose const & pose1, pose::Pose const & pose2,
		bool const verbose = true );

	bool check_matching_atom_names( pose::Pose const & pose1, pose::Pose const & pose2,
		std::map< id::AtomID, id::AtomID > const & atom_id_map, bool const verbose = true );

	utility::vector1< Size > const & rmsd_res_in_pose(){ return rmsd_res_in_pose_;}
	utility::vector1< Size > const & superimpose_res_in_pose(){ return superimpose_res_in_pose_;}

private:

	void
	superimpose_recursively( pose::Pose & pose,
		Real & rmsd_over_all, Size & natoms_over_all );

	Real
	do_superimposition( pose::Pose & pose );

	void
	update_reference_pose_local( pose::Pose const & pose );

	utility::vector1< Size >
	get_res_list_in_reference( pose::Pose const & pose ) const;

	void
	update_calc_rms_atom_id_map( pose::Pose const & pose );

	void
	get_calc_rms_atom_id_map( std::map< id::AtomID, id::AtomID > & calc_rms_atom_id_map,
		core::pose::Pose const & pose,
		utility::vector1 < Size > const & calc_rms_res ) const;

	void
	update_superimpose_atom_id_map( pose::Pose const & pose );

	Size
	get_rmsd_res_and_superimpose_res_in_pose( pose::Pose const & pose );

	void
	add_coordinate_constraints_from_map( pose::Pose & pose, pose::Pose const & reference_pose,
		std::map< id::AtomID, id::AtomID > const & atom_id_map,
		core::Real const & constraint_x0, core::Real const & constraint_tol ) const;

	bool
	do_checks( std::string const & atom_name, Size const & n, pose::Pose const & pose ) const;

	bool
	add_to_atom_id_map_after_checks( std::map< id::AtomID, id::AtomID> & atom_id_map,
		std::string const & atom_name,
		Size const & n1, Size const & n2,
		pose::Pose const & pose1, pose::Pose const & pose2,
		bool const do_the_checks = true ) const;

	std::map< id::AtomID, id::AtomID >
	get_root_triad_atom_id_map( pose::Pose const & pose ) const;

	void
	output_atom_id_map( std::map< id::AtomID, id::AtomID > const & atom_id_map ) const;

	void
	output_atom_id_map( std::map< id::AtomID, id::AtomID > const & atom_id_map,
		pose::Pose const & pose1,
		pose::Pose const & pose2 ) const;

	std::map< id::AtomID, id::AtomID > filter_virtual_atoms( std::map< id::AtomID, id::AtomID > const & atom_id_map,
		pose::Pose const & pose );

private:

	pose::Pose const & reference_pose_;
	pose::PoseOP reference_pose_local_;
	utility::vector1< Size > root_partition_res_;
	utility::vector1< Size > user_defined_calc_rms_res_;

	utility::vector1< Size > rmsd_res_in_pose_;
	utility::vector1< Size > superimpose_res_in_pose_;
	std::map< id::AtomID, id::AtomID > calc_rms_atom_id_map_;
	std::map< id::AtomID, id::AtomID > complete_moving_atom_id_map_;
	std::map< id::AtomID, id::AtomID > superimpose_atom_id_map_;
	std::string annotated_sequence_used_for_atom_id_maps_;

	Real rmsd_;
	Real superimpose_rmsd_;

	Real const check_alignment_tolerance_;
	bool superimpose_over_all_instantiated_;

};

} //align
} //modeler
} //stepwise
} //protocols

#endif
