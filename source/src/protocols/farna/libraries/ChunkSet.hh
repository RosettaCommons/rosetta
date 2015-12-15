// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/libraries/ChunkSet.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_libraries_ChunkSet_HH
#define INCLUDED_protocols_farna_libraries_ChunkSet_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/farna/libraries/ChunkSet.fwd.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>
#ifdef WIN32
#include <core/pose/MiniPose.hh>
#endif
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>
#include <protocols/toolbox/AtomID_Mapper.fwd.hh>


namespace protocols {
namespace farna {
namespace libraries {

class ChunkSet : public utility::pointer::ReferenceCount  {
public:

	//constructor!
	ChunkSet( utility::vector1< core::pose::PoseOP > const & pose_list,
		core::pose::ResMap const & res_map );

	ChunkSet( utility::vector1< core::pose::MiniPoseOP > const & mini_pose_list,
		core::pose::Pose const & example_pose,
		core::pose::ResMap const & res_map );


	// Need a clone();

	//destructor -- necessary?
	~ChunkSet();

	void
	insert_chunk_into_pose( core::pose::Pose & pose, Size const & chunk_pose_index, protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
													bool do_rosetta_library_domain_check = true ) const;

	Size
	num_chunks() const{ return mini_pose_list_.size(); };

	std::map< core::id::AtomID, core::id::AtomID >
	get_atom_id_map(  core::pose::Pose & pose, protocols::toolbox::AtomID_Mapper const & atom_id_mapper_to_target_vanilla_pose) const;

	core::pose::MiniPoseOP const mini_pose( Size const idx ) const;

	bool
	check_fold_tree_OK( core::pose::Pose const & pose ) const;

	void set_user_input( bool const & setting ){ user_input_ = setting; }
	bool user_input() const { return user_input_; }

private:

	void
	setup_atom_id_mask_and_mapper( core::pose::Pose const & pose );

	void
	setup_atom_id_mask( core::pose::Pose const & pose );

	std::map< core::id::AtomID, core::Size >
	get_atom_id_domain_map_for_rosetta_library_chunk(
												 std::map< core::id::AtomID, core::id::AtomID > atom_id_map,
												 core::pose::Pose const & pose, toolbox::AtomLevelDomainMap const & atom_level_domain_map,
												 bool do_rosetta_library_domain_check = true ) const;

	void filter_atom_id_map_with_mask( std::map< core::id::AtomID, core::id::AtomID > & atom_id_map ) const;

	void setup_atom_id_mapper_to_vanilla_chunk_pose( core::pose::Pose const & pose );

	void filter_poses_have_same_sequence_and_variants();

private:

	utility::vector1< core::pose::MiniPoseOP > mini_pose_list_;
	core::pose::ResMap res_map_; // goes from big pose into chunk (mini) pose.
	std::map< core::id::AtomID, bool > atom_id_mask_;
	bool user_input_;
	toolbox::AtomID_MapperCOP atom_id_mapper_to_vanilla_chunk_pose_;

};


} //libraries
} //farna
} //protocols

#endif
