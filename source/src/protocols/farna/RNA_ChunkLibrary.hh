// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Rhiju Das

#ifndef INCLUDED_protocols_rna_RNA_ChunkLibrary_HH
#define INCLUDED_protocols_rna_RNA_ChunkLibrary_HH

#include <protocols/farna/RNA_ChunkLibrary.fwd.hh>
#include <protocols/farna/BasePairStepLibrary.fwd.hh>
#include <protocols/farna/BasePairStep.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>
#ifdef WIN32
#include <core/pose/MiniPose.hh>
#endif
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>
#include <protocols/toolbox/AllowInsert.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// C++ Headers
#include <string>
#include <map>

#include <utility/vector1.hh>


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//
// This is a kind of generalized fragment library, where mini poses that
//  have putative solutions for junctions, loops, internal loops etc. can be switched into a
//  larger pose. Does not require ideal geometry
// Note that the class is actually not that RNA specific, so perhaps may be useful for others,
//  although I'm currently only developing it for pieces connected to RNA stems.
// Maybe has overlap with CartesianFragment class used by DNA people, but I find this somewhat
//  less confusing.
// The only thing to be wary of ... the fold tree must be set up with chainbreaks in the "right places"...
//  basically one at every junction.
//
//
namespace protocols {
namespace farna {

class ChunkSet : public utility::pointer::ReferenceCount  {
public:

	//constructor!
	ChunkSet( utility::vector1< core::pose::MiniPoseOP > const & mini_pose_list,
		core::pose::ResMap const & res_map );

	ChunkSet( utility::vector1< core::pose::PoseOP > const & pose_list,
		core::pose::ResMap const & res_map );

	// Need a clone();

	//destructor -- necessary?
	~ChunkSet();

	void
	insert_chunk_into_pose( core::pose::Pose & pose, Size const & chunk_pose_index, protocols::toolbox::AllowInsertOP const & allow_insert ) const;

	Size
	num_chunks() const{ return mini_pose_list_.size(); };

	std::map< core::id::AtomID, core::id::AtomID >
	get_atom_id_map(  core::pose::Pose & pose, protocols::toolbox::AllowInsertOP const & allow_insert ) const;

	core::pose::MiniPoseOP const mini_pose( Size const idx ) const;

	bool
	check_fold_tree_OK( core::pose::Pose const & pose );

private:

	void filter_atom_id_map_with_mask( std::map< core::id::AtomID, core::id::AtomID > & atom_id_map ) const;

	utility::vector1< core::pose::MiniPoseOP > mini_pose_list_;
	core::pose::ResMap res_map_;
	std::map< core::id::AtomID, bool > atom_id_mask_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class RNA_ChunkLibrary : public utility::pointer::ReferenceCount  {
public:

	RNA_ChunkLibrary();

	// default constructor.
	RNA_ChunkLibrary(
		utility::vector1 < std::string > const & pdb_files,
		utility::vector1 < std::string > const & silent_files,
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & input_res );

	// should not be in use in the future...
	RNA_ChunkLibrary(
		utility::vector1 < std::string > const & silent_files,
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & input_res );

	//destructor
	~RNA_ChunkLibrary();

	// default constructor.
	void
	initialize_rna_chunk_library(
		utility::vector1 < std::string > const & pdb_files,
		utility::vector1 < std::string > const & silent_files,
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & input_res );

	Size num_chunk_sets() const { return chunk_sets_.size(); };

	Size num_chunks( Size const n ) const { return chunk_sets_[ n ]->num_chunks(); };

	ChunkSetOP chunk_set( Size const n ) const { return chunk_sets_[ n ];  };

	void add_chunk_set( std::string const & silent_file,
		core::pose::ResMap const & res_map,
		core::pose::Pose const & big_pose );

	void insert_chunk_into_pose( core::pose::Pose & pose,
		Size const & chunk_list_index,
		Size const & chunk_pose_index ) const;

	utility::vector1< Size >
	get_indices_of_moving_chunks() const;

	Size num_moving_chunks() const ;

	bool
	random_chunk_insertion( core::pose::Pose & pose ) const;

	void
	initialize_random_chunks( core::pose::Pose & pose, bool const dump_pdb = false ) const;

	toolbox::AllowInsertOP allow_insert(){ return allow_insert_; };

	void set_allow_insert(toolbox::AllowInsertOP allow_insert );

	core::Real const & chunk_coverage() const{ return chunk_coverage_; };

	void
	superimpose_to_first_chunk( core::pose::Pose & pose ) const;

	bool
	check_fold_tree_OK( core::pose::Pose const & pose );

	void
	setup_base_pair_step_chunks( core::pose::Pose const & pose, utility::vector1< BasePairStep > base_pair_steps );

private:

	void
	update_allow_insert( core::pose::ResMap const & res_map,
		core::pose::Pose const & pose,
		core::pose::Pose const & scratch_pose,
		Size const domain_num );

	bool
	check_fold_tree_OK( core::pose::ResMap const & res_map,
		core::pose::Pose const & pose,
		core::pose::Pose const & scratch_pose );
	void
	figure_out_chunk_coverage();

	bool
	check_res_map( core::pose::ResMap const & res_map, core::pose::Pose const & scratch_pose, std::string const & sequence ) const;

	void
	align_to_chunk( core::pose::Pose & pose, ChunkSet const & chunk_set, core::Size const chunk_index ) const;

private:

	utility::vector1< ChunkSetOP > chunk_sets_;
	toolbox::AllowInsertOP allow_insert_;
	ObjexxFCL::FArray1D <bool> covered_by_chunk_;
	core::Real chunk_coverage_;
	bool coarse_rna_;
	BasePairStepLibraryOP base_pair_step_library_;

};


} //farna
} //protocols

#endif
