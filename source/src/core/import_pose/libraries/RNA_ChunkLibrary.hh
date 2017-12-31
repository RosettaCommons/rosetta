// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Rhiju Das

#ifndef INCLUDED_core_import_pose_RNA_ChunkLibrary_HH
#define INCLUDED_core_import_pose_RNA_ChunkLibrary_HH

#include <core/import_pose/libraries/RNA_ChunkLibrary.fwd.hh>
#include <core/import_pose/libraries/ChunkSet.fwd.hh>
#include <core/import_pose/libraries/BasePairStepLibrary.fwd.hh>
#include <core/pose/rna/BasePairStep.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/pose/toolbox/AtomID_Mapper.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// C++ Headers
#include <string>
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace import_pose {
namespace libraries {

extern core::Size const ROSETTA_LIBRARY_DOMAIN;

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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class RNA_ChunkLibrary : public utility::pointer::ReferenceCount  {
public:

	RNA_ChunkLibrary();

	RNA_ChunkLibrary( utility::vector1 < std::string > const & pdb_files,
		utility::vector1 < std::string > const & silent_files,
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & input_res,
		utility::vector1< core::Size > const & allow_insert_res = utility::vector1< core::Size >() );

	// should not be in use in the future...
	RNA_ChunkLibrary(
		utility::vector1 < std::string > const & silent_files,
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & input_res );

	RNA_ChunkLibrary( core::pose::Pose const & pose );

	//destructoro
	~RNA_ChunkLibrary();

	/// @brief clone the ChunkLibrary
	RNA_ChunkLibraryOP
	clone() const;

	// default constructor.
	void
	initialize_rna_chunk_library(
		utility::vector1 < std::string > const & pdb_files,
		utility::vector1 < std::string > const & silent_files,
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & input_res,
		utility::vector1< core::Size > const & allow_insert_res = utility::vector1< core::Size >() /* blank */ );

	Size num_chunk_sets() const { return chunk_sets_.size(); };

	Size num_chunks( Size const n ) const;

	ChunkSetOP chunk_set( Size const n ) const { return chunk_sets_[ n ];  };

	utility::vector1< ChunkSetOP > chunk_sets() const { return chunk_sets_;  };

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

	core::Size
	get_alignment_domain( core::pose::Pose const & pose ) const;

	void
	initialize_random_chunks( core::pose::Pose & pose, bool const dump_pdb = false ) const;

	void
	insert_random_protein_chunks( core::pose::Pose & pose ) const;

	core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map() const { return atom_level_domain_map_; };


	void set_atom_level_domain_map(core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map );

	core::Real const & chunk_coverage() const { return chunk_coverage_; };

	Size
	single_user_input_chunk() const;

	bool
	superimpose_to_single_user_input_chunk( core::pose::Pose & pose ) const;

	bool
	check_fold_tree_OK( core::pose::Pose const & pose ) const;

	void
	setup_base_pair_step_chunks( core::pose::Pose const & pose,
		utility::vector1< core::pose::rna::BasePairStep > const & base_pair_steps,
		BasePairStepLibrary const & base_pair_step_library );

	void
	update_to_move_rosetta_library_chunks();

private:

	void
	update_atom_level_domain_map( core::pose::ResMap const & res_map,
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
	core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map_;
	ObjexxFCL::FArray1D <bool> covered_by_chunk_;
	core::Real chunk_coverage_;
	bool coarse_rna_;
	bool do_rosetta_library_domain_check_;

};


} //libraries
} //denovo
} //protocols

#endif
