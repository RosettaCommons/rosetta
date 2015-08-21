// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/BasePairStepLibrary.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/BasePairStepLibrary.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <string>

static thread_local basic::Tracer TR( "protocols.rna.BasePairStepLibrary" );

namespace protocols {
namespace farna {

//constructor
BasePairStepSequence::BasePairStepSequence( char const nt_i, char const nt_i_next, char const nt_j, char const nt_j_next ){
	base_pair_step_sequence_ = std::make_pair(
		std::make_pair( nt_i, nt_i_next ),
		std::make_pair( nt_j, nt_j_next ) );
}

//constructor
BasePairStepSequence::BasePairStepSequence( std::string const & sequence,
	Size const i, Size const i_next, Size const j, Size const j_next ){
	base_pair_step_sequence_ = std::make_pair(
		std::make_pair( sequence[ i-1 ], sequence[ i_next-1 ] ),
		std::make_pair( sequence[ j-1 ], sequence[ j_next-1 ] ) );
}

//constructor
BasePairStepSequence::BasePairStepSequence( std::string const & sequence, BasePairStep const & base_pair_step ){
	base_pair_step_sequence_ = std::make_pair(
		std::make_pair( sequence[ base_pair_step.i() - 1 ],
		sequence[ base_pair_step.i_next() - 1 ] ),
		std::make_pair( sequence[ base_pair_step.j() - 1 ],
		sequence[ base_pair_step.j_next() - 1 ] ) );
}

//Constructor
BasePairStepLibrary::BasePairStepLibrary():
	initialized_( false )
{}

//Destructor
BasePairStepLibrary::~BasePairStepLibrary(){}

//////////////////////////////////////////////////////////
void
BasePairStepLibrary::initialize(){
	using namespace core::pose;

	if ( initialized_ ) return;

	utility::vector1< std::string > rna_base_pairs = utility::tools::make_vector1( "au","ua","cg","gc","gu","ug" );
	for ( Size m = 1; m <= rna_base_pairs.size(); m++ ) {
		for ( Size n = 1; n <= rna_base_pairs.size(); n++ ) {

			std::string bps_seq;
			bps_seq += rna_base_pairs[m][0];
			bps_seq += rna_base_pairs[n][0];
			bps_seq += "_";
			bps_seq += rna_base_pairs[n][1];
			bps_seq += rna_base_pairs[m][1];

			std::string input_file( basic::database::full_name("sampling/rna/base_pair_steps/" + bps_seq + ".out" ) );

			if ( !utility::file::file_exists( input_file ) ) {
				TR << "WARNING! COULD NOT FIND:  "<< input_file << std::endl;
				continue;
			}

			utility::vector1< pose::PoseOP > pose_list;
			process_input_file( input_file, pose_list, false /*is_pdb*/ );

			PoseOP const & scratch_pose = pose_list[ 1 ];
			utility::vector1< pose::MiniPoseOP > mini_pose_list;
			for ( Size k = 1; k <= pose_list.size(); k++ ) {
				mini_pose_list.push_back( core::pose::MiniPoseOP( new core::pose::MiniPose( *(pose_list[k]) ) ) );
			}

			char nt_i      = rna_base_pairs[m][0];
			char nt_i_next = rna_base_pairs[n][0];
			char nt_j      = rna_base_pairs[n][1];
			char nt_j_next = rna_base_pairs[m][1];
			BasePairStepSequence base_pair_step_sequence( nt_i, nt_i_next, nt_j, nt_j_next );
			mini_pose_lists_[ base_pair_step_sequence ] = mini_pose_list;
			scratch_poses_[ base_pair_step_sequence ] = scratch_pose;

		}
	}

	initialized_ = true;
}

//////////////////////////////////////////////////////////
bool
BasePairStepLibrary::has_value( BasePairStepSequence const & base_pair_step_sequence ) const{
	return ( mini_pose_lists_.find( base_pair_step_sequence ) != mini_pose_lists_.end() );
}

//////////////////////////////////////////////////////////
// List of poses for each base pair step -- converted to
// 'mini-pose' format with just coordinates.
utility::vector1< core::pose::MiniPoseOP > const &
BasePairStepLibrary::mini_pose_list( BasePairStepSequence const & base_pair_step_sequence ){
	runtime_assert( initialized_ );
	runtime_assert( has_value( base_pair_step_sequence ) );
	return mini_pose_lists_[ base_pair_step_sequence ];
}

//////////////////////////////////////////////////////////
// example of a full pose for each base pair step
pose::PoseOP const &
BasePairStepLibrary::scratch_pose( BasePairStepSequence const & base_pair_step_sequence ){
	runtime_assert( initialized_ );
	runtime_assert( has_value( base_pair_step_sequence ) );
	return scratch_poses_[ base_pair_step_sequence ];
}

} //farna
} //protocols
