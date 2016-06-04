// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/fldsgn/CircularPermutation.cc
/// @brief  perform circularpermutation give a pose ( under development )
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// unit headers
#include <protocols/fldsgn/CircularPermutation.hh>
#include <protocols/fldsgn/CircularPermutationCreator.hh>

// projects headers
#include <core/conformation/Conformation.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>

#include <basic/Tracer.hh>

#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>


#include <protocols/forge/build/SegmentSwap.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/Interval.hh>

#include <protocols/forge/methods/fold_tree_functions.hh>

// Headers for Paser block
#include <utility/tag/Tag.hh>

#include <protocols/forge/build/SegmentRebuild.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
//#include <basic/datacache/DataMap.hh>


namespace protocols {
namespace fldsgn {

static THREAD_LOCAL basic::Tracer TR( "protocols.fldsgn.CircularPermutation" );

std::string
CircularPermutationCreator::keyname() const
{
	return CircularPermutationCreator::mover_name();
}

protocols::moves::MoverOP
CircularPermutationCreator::create_mover() const {
	return protocols::moves::MoverOP( new CircularPermutation );
}

std::string
CircularPermutationCreator::mover_name()
{
	return "CircularPermutation";
}


/// @brief default constructor
CircularPermutation::CircularPermutation() :
	Mover( "CircularPermutation" ),
	new_terminal_pos_( 0 ),
	ignore_chain_( false ),
	split_( 0 )
{}


/// @Brief copy constructor
CircularPermutation::CircularPermutation( CircularPermutation const & rval ) :
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	new_terminal_pos_( rval.new_terminal_pos_ ),
	split_( rval.split_ )
{}


/// @brief default destructor
CircularPermutation::~CircularPermutation() {}


/// @brief clone this object
CircularPermutation::MoverOP
CircularPermutation::clone() const
{
	return CircularPermutation::MoverOP( new CircularPermutation( *this ) );
}


/// @brief create this type of object
CircularPermutation::MoverOP
CircularPermutation::fresh_instance() const
{
	return CircularPermutation::MoverOP( new CircularPermutation() );
}


/// @brief new N- & C- terminal position
void
CircularPermutation::new_terminal_pos( Size const s )
{
	new_terminal_pos_ = s;
}


/// @brief total number of cycles
CircularPermutation::Size
CircularPermutation::new_terminal_pos() const
{
	return new_terminal_pos_;
}

/// @brief total number of cycles
CircularPermutation::Size
CircularPermutation::which_chain( Size const pos, Pose const & pose ) const
{
	for ( Size i=1; i<=pose.conformation().num_chains(); i++ ) {

		Size chain_begin( pose.conformation().chain_begin( i ) );
		Size chain_end( pose.conformation().chain_end( i ) );

		if ( chain_begin <= pos && chain_end >= pos ) {
			return i;
		}

	}
	return 0;
}

std::string
CircularPermutation::get_name() const {
	return CircularPermutationCreator::mover_name();
}


/// @brief
void
CircularPermutation::split_chains( Pose & pose, utility::vector1< Size > const & pos )
{

	using core::kinematics::FoldTree;
	using core::kinematics::MoveMap;
	using core::pose::Pose;
	using protocols::forge::methods::fold_tree_from_pose;

	// information of original fold tree completely overridded
	FoldTree ft = fold_tree_from_pose( pose, pose.fold_tree().root(), MoveMap() );

	// insert chain endings
	for ( Size i=1; i<=pos.size(); i++ ) {
		Size position = pos[ i ];
		runtime_assert( position >= 1 && position <= pose.total_residue() );
		pose.conformation().insert_chain_ending( position );
	}

	// add terminals
	for ( Size i=1; i<=pose.conformation().num_chains(); i++ ) {
		Size begin = pose.conformation().chain_begin( i );
		Size end = pose.conformation().chain_end( i );
		if ( pose.residue( begin ).is_protein() && !pose.residue( begin ).is_lower_terminus() ) {
			core::pose::add_lower_terminus_type_to_pose_residue( pose, begin );
		}
		if ( pose.residue( end ).is_protein() && !pose.residue( end ).is_upper_terminus() ) {
			core::pose::add_upper_terminus_type_to_pose_residue( pose, end );
		}
	}

	// add jumps
	for ( Size i=2; i<=pose.conformation().num_chains(); i++ ) {
		Size begin = pose.conformation().chain_begin( i );
		Size end = pose.conformation().chain_end( i-1 );
		ft.new_jump( end, begin, end );
	}

	pose.fold_tree( ft );

}


/// @brief main apply
void CircularPermutation::apply( Pose & pose )
{
	using protocols::forge::build::BuildManager;
	using protocols::forge::build::SegmentSwap;
	using protocols::forge::build::SegmentRebuild;
	using protocols::forge::build::Interval;
	//using core::chemical::ResidueTypeSet;
	using core::chemical::ResidueTypeSetCOP;
	using core::conformation::ResidueOP;
	using core::conformation::ResidueFactory;
	using core::kinematics::FoldTree;
	using core::pose::symmetry::is_symmetric;
	using core::pose::symmetry::make_asymmetric_pose;
	using core::pose::symmetry::get_asymmetric_pose_copy_from_symmetric_pose;

	// make sure position of new terminal within chain lengths
	if ( new_terminal_pos_ == 0 ) {
		TR << "No position to be permutated was specified. " << std::endl;
		return;
	}
	runtime_assert( new_terminal_pos_ > 1 && new_terminal_pos_ <= pose.total_residue() );

	// make pose asymmetric if pose is symmetric
	if ( is_symmetric( pose ) ) {
		Pose work_pose( pose );
		pose = get_asymmetric_pose_copy_from_symmetric_pose( work_pose );
	}

	// prepare pose for latter swapped in region
	Pose swap_in( pose );

	// determine begin and end positions of chain to be swapped
	Size chain_begin;
	Size chain( 0 );
	if ( ignore_chain_ ) {

		chain_begin = 1;
		// find final chains
		for ( Size i=1; i<=pose.total_residue(); i++ ) {
			if ( pose.residue( i ).is_protein() && static_cast< int >( chain )< pose.chain( i ) ) {
				chain = pose.chain( i );
			}
		}
		runtime_assert( chain != 0 );

	} else {
		chain = which_chain( new_terminal_pos_, pose );
		chain_begin = pose.conformation().chain_begin( chain );
	}
	Size chain_end = pose.conformation().chain_end( chain );

	// add 4 residues at the end of chain
	ResidueTypeSetCOP rsd_set( pose.residue(1).residue_type_set() );
	ResidueOP ala( ResidueFactory::create_residue( rsd_set->name_map( "ALA" ) ) );
	for ( Size i=1; i<=4; i++ ) {
		Size pos = pose.conformation().chain_end( chain );
		pose.conformation().safely_append_polymer_residue_after_seqpos( *ala, pos, true );
	}

	// Blow away the former parts of pose
	pose.conformation().delete_residue_range_slow( chain_begin, new_terminal_pos_ - 1 );

	// Blow away the latter parts of swap_in pose
	swap_in.conformation().delete_residue_range_slow( new_terminal_pos_, swap_in.total_residue() );
	if ( chain_begin > 1 ) {
		// Blow away the former parts of swap_in pose
		swap_in.conformation().delete_residue_range_slow( 1, chain_begin - 1 );
	}

	// swap_in swapped into the 3 residues added at the terminal of the chain
	MoveMap movemap;
	Size nt = chain_end - ( new_terminal_pos_ - chain_begin ) + 1;
	BuildManager manager;
	manager.add( forge::build::BuildInstructionOP( new SegmentSwap( Interval( nt, nt+2 ), movemap, swap_in ) ) );
	manager.modify( pose );

	// delete the final residues added at the terminal of the chain
	pose.conformation().delete_residue_range_slow( chain_end + 1,  chain_end + 1 );

	// re-add terminus
	if ( !pose.residue( chain_end - 1 ).is_upper_terminus() ) {
		core::pose::add_upper_terminus_type_to_pose_residue( pose, chain_end );
	}

	// re-add terminus
	if ( !pose.residue( chain_begin ).is_lower_terminus() ) {
		core::pose::add_lower_terminus_type_to_pose_residue( pose, chain_begin );
	}

	TR << "FoldTree after circular permutation: " << pose.fold_tree() << std::endl;

	if ( split_ != 0 ) {
		utility::vector1< Size > positions;
		positions.push_back( split_ );
		pose.conformation().reset_chain_endings();
		split_chains( pose, positions );
	}

	TR << pose.fold_tree() << std::endl;

}


/// @brief parse xml
void
CircularPermutation::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{

	// set positions of new N- and C- terminal
	new_terminal_pos_ = ( tag->getOption<Size>( "pos", 0 ) );

	// ignore chain
	ignore_chain_ = ( tag->getOption<bool>( "ignore_chain", 0 ) );

	// split chain
	split_ = ( tag->getOption<Size>( "split_chain", 0 ) );

}


} // namespace fldsgn
} // namespace protocols
