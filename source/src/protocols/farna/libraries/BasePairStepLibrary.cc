// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/libraries/BasePairStepLibrary.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/libraries/BasePairStepLibrary.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/MiniPose.hh>
#include <core/chemical/util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <string>

static THREAD_LOCAL basic::Tracer TR( "protocols.farna.libraries.BasePairStepLibrary" );

using namespace protocols::farna::base_pairs;
using namespace core;

namespace protocols {
namespace farna {
namespace libraries {

//constructor
BasePairStepSequence::BasePairStepSequence( char const nt_i,
	char const nt_i_next,
	char const nt_j,
	char const nt_j_next,
	Size num_bulge /* = 0 */){
	base_pair_step_sequence_ = std::make_pair(
		std::make_pair( nt_i, nt_i_next ),
		std::make_pair( nt_j, nt_j_next ) );
	num_bulge_ = num_bulge;
}

//constructor
BasePairStepSequence::BasePairStepSequence( std::string const & sequence,
	Size const i, Size const i_next, Size const j, Size const j_next ){
	base_pair_step_sequence_ = std::make_pair(
		std::make_pair( sequence[ i-1 ], sequence[ i_next-1 ] ),
		std::make_pair( sequence[ j-1 ], sequence[ j_next-1 ] ) );
	runtime_assert( j_next > j );
	num_bulge_ = j_next - j - 1;
}

//constructor
BasePairStepSequence::BasePairStepSequence( std::string const & sequence,
	BasePairStep const & base_pair_step ){
	base_pair_step_sequence_ = std::make_pair(
		std::make_pair( sequence[ base_pair_step.i() - 1 ],
		sequence[ base_pair_step.i_next() - 1 ] ),
		std::make_pair( sequence[ base_pair_step.j() - 1 ],
		sequence[ base_pair_step.j_next() - 1 ] ) );
	runtime_assert( base_pair_step.j_next() > base_pair_step.j() );
	num_bulge_ = base_pair_step.j_next()  - base_pair_step.j() - 1;
}

//////////////
std::string
BasePairStepSequence::subdir() const
{
	if ( num_bulge_ == 0 ) {
		return "standard/";
	} else if ( num_bulge_ >= 1 && num_bulge_ <= (MAX_BULGE_LENGTH) ) {
		return "bulge_"+ObjexxFCL::format::I( 1, num_bulge_ )+"nt/";
	} else {
		return "long_distance/";
	}
}


//Constructor
BasePairStepLibrary::BasePairStepLibrary( bool const canonical /* = true */ ):
	canonical_( canonical )
{
	initialize();
}

//Destructor
BasePairStepLibrary::~BasePairStepLibrary(){}

//////////////////////////////////////////////////////////
void
BasePairStepLibrary::initialize(){
	using namespace core::pose;
	using core::chemical::rna::rna_nts;
	utility::vector1< std::string > rna_base_pairs;
	if ( canonical_ ) {
		rna_base_pairs = utility::tools::make_vector1( "au","ua","cg","gc","gu","ug" );
	} else {
		// check all possible base pairs
		for ( Size m = 1; m <= 4; m++ ) {
			for ( Size n = 1; n <= 4; n++ ) {
				std::string bp;
				bp += rna_nts[m-1];
				bp += rna_nts[n-1];
				rna_base_pairs.push_back( bp );
			}
		}
	}

	Size const max_num_bulge = canonical_ ? 0 : MAX_BULGE_LENGTH;
	for ( Size q = 0; q <= max_num_bulge; q++ ) {
		for ( Size m = 1; m <= rna_base_pairs.size(); m++ ) {
			for ( Size n = 1; n <= rna_base_pairs.size(); n++ ) {
				char nt_i      = rna_base_pairs[m][0];
				char nt_i_next = rna_base_pairs[n][0];
				char nt_j      = rna_base_pairs[n][1];
				char nt_j_next = rna_base_pairs[m][1];
				BasePairStepSequence base_pair_step_sequence( nt_i, nt_i_next, nt_j, nt_j_next, q );
				initialize_data( base_pair_step_sequence, false /* lazy loading when actually needed*/ );
			}
		}
	}
}


//////////////////////////////////////////////////////////
// would be pretty easy to make this user-settable.
std::string
BasePairStepLibrary::database_dir() const {
	if ( canonical_ ) return basic::database::full_name( "sampling/rna/base_pair_steps/canonical/" );
	return  basic::database::full_name( "sampling/rna/base_pair_steps/general/" );
}

//////////////////////////////////////////////////////////
// although this is const, it can update the mutable list of mini_poses...
void
BasePairStepLibrary::initialize_data( BasePairStepSequence const & base_pair_step_sequence,
	bool const load_in_poses /* = true */ ) const
{
	// already initialized?
	if ( has_value( base_pair_step_sequence ) && mini_pose_lists_[ base_pair_step_sequence ].size() > 0 ) return; // already initialized.

	std::string input_file( database_dir() + base_pair_step_sequence.subdir() + base_pair_step_sequence.tag() + ".out" );
	if ( !utility::file::file_exists( input_file ) ) {
		input_file += ".gz";
		if ( !utility::file::file_exists( input_file ) ) {
			//   TR << "WARNING! COULD NOT FIND:  "<< input_file << std::endl;
			return;
		}
	}

	if ( load_in_poses ) {

		utility::vector1< pose::PoseOP > pose_list_raw, pose_list;
		process_input_file( input_file, pose_list_raw, false /*is_pdb*/ );

		// there appears to be an issue with some of the poses in the database. Oops!
		for ( core::pose::PoseOP const & pose : pose_list_raw ) {
			core::kinematics::FoldTree const & f( pose->fold_tree() );
			if ( ( f.is_cutpoint( 1 ) ) ||
					( !f.is_cutpoint( 2 ) ) ||
					( base_pair_step_sequence.num_bulge() == 0 && f.is_cutpoint( 3 ) ) ||
					( base_pair_step_sequence.num_bulge() >  0 && !f.is_cutpoint( 3 ) ) ) {
				TR.Warning << TR.Red << "screwed up fold tree " << pose->fold_tree() << " for base pair step sequence " << base_pair_step_sequence.tag() << " in pose " << tag_from_pose( *pose ) << " in input_file " << input_file << std::endl;
				continue;
			}
			pose_list.push_back( pose );
		}

		utility::vector1< pose::MiniPoseOP > mini_pose_list;
		for ( auto const & pose : pose_list ) {
			mini_pose_list.push_back( pose::MiniPoseOP( new core::pose::MiniPose( *pose ) ) );
		}

		// potential problem -- what if all poses are bad?
		runtime_assert( mini_pose_list.size() > 0 );
		mini_pose_lists_[ base_pair_step_sequence ] = mini_pose_list;
		if ( pose_list.size() > 0 ) scratch_poses_[   base_pair_step_sequence ] = pose_list[ 1 ];

	} else {
		mini_pose_lists_[ base_pair_step_sequence ] = utility::vector1< pose::MiniPoseOP >(); // empty
		scratch_poses_[   base_pair_step_sequence ] = 0;
	}
}

//////////////////////////////////////////////////////////
bool
BasePairStepLibrary::has_value( BasePairStepSequence const & base_pair_step_sequence ) const
{
	return ( mini_pose_lists_.find( base_pair_step_sequence ) != mini_pose_lists_.end() );
}

//////////////////////////////////////////////////////////
// List of poses for each base pair step -- converted to
// 'mini-pose' format with just coordinates.
utility::vector1< core::pose::MiniPoseOP > const &
BasePairStepLibrary::mini_pose_list( BasePairStepSequence const & base_pair_step_sequence ) const
{
	runtime_assert( has_value( base_pair_step_sequence ) );
	initialize_data( base_pair_step_sequence ); // makes sure initialized -- no op if already initialized
	// TR << base_pair_step_sequence << " has " <<  mini_pose_lists_[ base_pair_step_sequence ].size() << " members " << std::endl;
	return mini_pose_lists_[ base_pair_step_sequence ];
}

//////////////////////////////////////////////////////////
// example of a full pose for each base pair step
pose::PoseCOP const &
BasePairStepLibrary::scratch_pose( BasePairStepSequence const & base_pair_step_sequence ) const
{
	runtime_assert( has_value( base_pair_step_sequence ) );
	initialize_data( base_pair_step_sequence );
	return scratch_poses_[ base_pair_step_sequence ];
}

} //libraries
} //farna
} //protocols
