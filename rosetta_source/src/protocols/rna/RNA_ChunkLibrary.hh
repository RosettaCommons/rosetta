// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//  CVS information:
//  $Revision: 1.1.2.1 $
//  $Date: 2005/11/07 21:05:35 $
//  $Author: rhiju $
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_rna_RNA_ChunkLibrary_HH
#define INCLUDED_protocols_rna_RNA_ChunkLibrary_HH

#include <protocols/rna/RNA_ChunkLibrary.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>
#ifdef WIN32
#include <core/pose/MiniPose.hh>
#endif
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>
#include <numeric/xyzVector.hh>
#include <protocols/rna/AllowInsert.fwd.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/vector1.fwd.hh>
#include <core/id/AtomID.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// C++ Headers
#include <string>
#include <map>
// AUTO-REMOVED #include <vector>

//Auto Headers
#include <core/id/AtomID.fwd.hh>
#include <utility/vector1_bool.hh>
#include <iostream>



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
//  basically one at eery junction.
//
//

namespace protocols{
namespace rna{

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
		insert_chunk_into_pose( core::pose::Pose & pose, Size const & chunk_pose_index, protocols::rna::AllowInsertOP const & allow_insert ) const;

		Size
		num_chunks() const{ return mini_pose_list_.size(); };

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

		// constructor -- needs a list of silent files. Each silent file
		//  has solutions for a particular piece of the desired pose.
		RNA_ChunkLibrary( utility::vector1 < std::string > const & silent_files,
											core::pose::Pose const & pose,
											std::map< Size, Size > const & connections_in_big_pose /* to figure out mapping to big pose*/
											);

		RNA_ChunkLibrary(
								utility::vector1 < std::string > const & silent_files,
								core::pose::Pose const & pose,
								utility::vector1< core::Size > const & input_res );

		//destructor -- necessary?
		// ~RNA_ChunkLibrary();

 		Size num_chunk_sets() const { return chunk_sets_.size(); };

 		Size num_chunks( Size const n ) const { return chunk_sets_[ n ]->num_chunks(); };

		ChunkSetOP chunk_set( Size const n ) const { return chunk_sets_[ n ];  };

		void add_chunk_set( std::string const & silent_file,
												core::pose::ResMap const & res_map,
												core::pose::Pose const & big_pose );

		void insert_chunk_into_pose( core::pose::Pose & pose,
																 Size const & chunk_list_index,
																 Size const & chunk_pose_index ) const;

		void
		random_chunk_insertion( core::pose::Pose & pose ) const;

		void
		initialize_random_chunks( core::pose::Pose & pose, bool const dump_pdb = false ) const;

		AllowInsertOP allow_insert(){ return allow_insert_; };

		void set_allow_insert( AllowInsertOP allow_insert );

		core::Real const & chunk_coverage() const{ return chunk_coverage_; };

	private:

		void
		zero_out_allow_insert( core::pose::ResMap const & res_map,
													 core::pose::Pose const & pose,
													 core::pose::Pose const & scratch_pose,
													 Size const domain_num );


		void
		figure_out_chunk_coverage();

		void
		get_component_sequences(
							 core::pose::Pose const & pose,
							 utility::vector1< std::string > & sequences,
							 utility::vector1< core::Size > & chain_id,
							 utility::vector1< core::Size > & sequence_start ) const;

		void
		process_silent_file( std::string const & silent_file,
												 utility::vector1< core::pose::PoseOP > & pose_list ) const;

		void
		figure_out_possible_res_maps(
																 utility::vector1< core::pose::ResMap > & res_maps,
																 core::pose::Pose const & scratch_pose,
																 std::string const & sequence_of_big_pose,
																 std::map< Size, Size > const & connections_in_big_pose ) const;

		void
		find_res_maps(
									utility::vector1< Size > const & chain_id,
									utility::vector1< Size > const & scratch_sequence_start,
									utility::vector1< std::string > const & scratch_sequences,
									utility::vector1< utility::vector1< Size > > const & matches_to_each_scratch_sequence,
								core::pose::Pose const & scratch_pose,
									std::map< Size, Size > const & connections_in_big_pose,
									utility::vector1< core::pose::ResMap > & res_maps ) const;

			void
			get_sequence_matches( 	utility::vector1< utility::vector1< Size > > & matches_to_each_scratch_sequence,
															utility::vector1< std::string > const &	scratch_sequences,
															std::string const & sequence_of_big_pose ) const;

		void
		check_connections( Size const & num_chain, core::pose::ResMap & res_map,
											 utility::vector1< Size > const & chain_id,
											 utility::vector1< Size > const & scratch_sequence_start,
											 utility::vector1< std::string > const & scratch_sequences,
											 utility::vector1< utility::vector1< Size > > const & matches_to_each_scratch_sequence,
											 core::pose::Pose const & scratch_pose,
											 std::map< Size, Size > const & connections_in_big_pose,
											 utility::vector1< core::pose::ResMap > & res_maps ) const;

		bool
		fill_res_map( core::pose::ResMap & res_map, Size const & match_pos, Size const & scratch_start_pos, Size const & scratch_sequence_length ) const;

		void
		test_matches( Size const & res1, Size const & res2, core::pose::ResMap & res_map,
									utility::vector1< Size > const & chain_id,
									utility::vector1< Size > const & scratch_sequence_start,
									utility::vector1< std::string > const & scratch_sequences,
									utility::vector1< utility::vector1< Size > > const & matches_to_each_scratch_sequence,
									core::pose::Pose const & scratch_pose,
									std::map< Size, Size > const & connections_in_big_pose,
									utility::vector1< core::pose::ResMap > & res_maps ) const;

		bool
		check_res_map( core::pose::ResMap const & res_map, core::pose::Pose const & scratch_pose, std::string const & sequence ) const;

		void
		check_res_map_recursively( core::pose::ResMap const & res_map_old,
															 utility::vector1< std::string > const & scratch_sequences,
															 utility::vector1< utility::vector1< Size > > const & matches_to_each_scratch_sequence,
															 core::pose::Pose const & scratch_pose,
															 std::map< Size, Size > const & connections_in_big_pose,
															 utility::vector1< core::Size > const & chain_id,
															 core::Size const & num_sequence,
															 core::Size const & num_match,
															 utility::vector1< core::pose::ResMap > & res_maps ) const;
		bool
		check_jump_match(
										 core::pose::Pose const & scratch_pose,
										 std::map< core::Size, core::Size > const & connections_in_big_pose,
										 core::pose::ResMap const & res_map,
										 utility::vector1< Size > const & chain_id ) const;

	private:

		utility::vector1< ChunkSetOP > chunk_sets_;
		AllowInsertOP allow_insert_;
		ObjexxFCL::FArray1D <bool> covered_by_chunk_;
		core::Real chunk_coverage_;
		bool coarse_rna_;

	};



}
}

#endif
