// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/movers/RNA_JumpMover.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/movers/RNA_JumpMover.hh>
#include <protocols/farna/libraries/RNA_JumpLibrary.hh>
#include <protocols/farna/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.farna.movers.RNA_JumpMover" );

using namespace core;
using namespace core::pose::rna;
using namespace ObjexxFCL::format;

namespace protocols {
namespace farna {
namespace movers {

//constructor
RNA_JumpMover::RNA_JumpMover( RNA_JumpLibraryCOP rna_jump_library,
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map ):
	rna_jump_library_( rna_jump_library ),
	atom_level_domain_map_( atom_level_domain_map )
{
}

//Destructor
RNA_JumpMover::~RNA_JumpMover()
{}

/////////////////////////////////////////////////////////
bool
RNA_JumpMover::random_jump_change( pose::Pose & pose ) const
{

	using namespace core::conformation;
	using namespace core::id;

	// Any jumps in here to tweak?
	Size const num_jump = pose.num_jump();
	if ( num_jump == 0 ) return false;

	Size ntries( 0 );
	Size const MAX_TRIES( 1000 );

	Size which_jump( 0 );

	while ( ntries++ < MAX_TRIES ) {
		// Randomly pick one.
		which_jump = static_cast<Size>( numeric::random::rg().uniform() * num_jump ) + 1 ;

		// Check that we can actually insert here. At least one of the jump partners
		// should allow moves. (I guess the other one can stay fixed).
		Size const jump_pos1( pose.fold_tree().upstream_jump_residue( which_jump ) );
		Size const jump_pos2( pose.fold_tree().downstream_jump_residue( which_jump ) ); // Unused variable causes warning.

		Residue const & rsd1 = pose.residue( jump_pos1 );
		AtomID jump_atom_id1( rsd1.atom_index( default_jump_atom( rsd1 ) ), jump_pos1 );

		Residue const & rsd2 = pose.residue( jump_pos2 ); // Unused variable causes warning.
		AtomID jump_atom_id2( rsd2.atom_index( default_jump_atom( rsd2 ) ), jump_pos2 ); // Unused variable causes warning.

		if ( moveable_jump( jump_atom_id1, jump_atom_id2, *atom_level_domain_map_ ) ) break;

	}

	if ( ntries >= MAX_TRIES ) return false;

	// std::cout << "will try to change jump " << which_jump << std::endl;

	// Tweak it -- for connections between chains for which the residue is unknown
	//  consider alternative connection -- if there aren't any jumps of these types,
	//  following does nothing.
	sample_alternative_chain_connection( pose, which_jump );

	// Shift jump.
	bool success( false );
	add_new_RNA_jump( pose, which_jump, success );

	return success;

}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_JumpMover::add_new_RNA_jump(
	pose::Pose & pose,
	Size const & which_jump,
	bool & success ) const
{

	kinematics::FoldTree fold_tree( pose.fold_tree() ); //Make a copy.

	Size const jump_pos1( fold_tree.upstream_jump_residue( which_jump ) );
	Size const jump_pos2( fold_tree.downstream_jump_residue( which_jump ) );

	// Later can be smart about virtual residue. (rigid body jumps)
	if ( !pose.residue( jump_pos1 ).is_RNA() || !pose.residue( jump_pos2 ).is_RNA() ) return;

	BaseEdge e1( ANY_BASE_EDGE ), e2( ANY_BASE_EDGE);
	BaseDoubletOrientation o( ANY_BASE_DOUBLET_ORIENTATION );
	bool found_pairing( false );
	for ( Size n = 1;  n <= rna_pairing_list_.size(); n++ ) {
		BasePair const & pairing = rna_pairing_list_[n];
		if ( pairing.res1() == jump_pos1 && pairing.res2() == jump_pos2 ) {
			e1 = pairing.edge1();
			e2 = pairing.edge2();
			o  = pairing.orientation();
			found_pairing = true;
			break;
		}
		if ( pairing.res1() == jump_pos2 && pairing.res2() == jump_pos1 ) {
			e1 = pairing.edge2();
			e2 = pairing.edge1();
			o  = pairing.orientation();
			found_pairing = true;
			break;
		}
	}

	//The jump may have come from a "chain_connection", which
	// doesn't know parallel vs. antiparallel..
	if ( !found_pairing ) {
		// To save time following could be as runtime_assert statement, only
		// active in debug builds.
		if ( check_in_chain_connections( jump_pos1, jump_pos2 ) > 0 ) found_pairing = true;
	}

	if ( !found_pairing ) {
		utility_exit_with_message(  "Trouble finding foldtree pairing in input pairing list? " + I(3,jump_pos1)+' '+I(3,jump_pos2) );
	}

	bool const forward1 = check_forward_backward( pose, jump_pos1 );
	bool const forward2 = check_forward_backward( pose, jump_pos2 );

	std::string atom_name1, atom_name2;
	runtime_assert( rna_jump_library_ != 0 );
	kinematics::Jump const new_jump = rna_jump_library_->get_random_base_pair_jump(
		pose.residue(jump_pos1).name1(),
		pose.residue(jump_pos2).name1(),
		e1, e2, o,
		atom_name1, atom_name2,
		success,
		forward1, forward2 );

	if ( !success ) return; //shh don't do anything.

	fold_tree.set_jump_atoms( which_jump, atom_name1, atom_name2 );

	pose.fold_tree( fold_tree );

	pose.set_jump( which_jump, new_jump );

}

/////////////////////////////////////////////////////////////////////////////////////////////////
Size
RNA_JumpMover::check_in_chain_connections( Size const & pos1, Size const & pos2 ) const
{
	Size const num_chain_connections( chain_connections_.size() );

	for ( Size n = 1; n <= num_chain_connections ; n++ ) {
		utility::vector1 < Size > const & res_list1( chain_connections_[n].first );
		utility::vector1 < Size > const & res_list2( chain_connections_[n].second);

		if ( res_list1.has_value( pos1 ) && res_list2.has_value( pos2 ) ) return n;
		if ( res_list1.has_value( pos2 ) && res_list2.has_value( pos1 ) ) return n;

	}

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_JumpMover::check_forward_backward(
	pose::Pose & pose,
	Size const jump_pos ) const
{
	using namespace core::kinematics::tree;

	if ( pose.residue( jump_pos ).is_coarse() || !pose.residue( jump_pos ).is_RNA() )  return true;

	Size atomno = pose.residue( jump_pos ).atom_index( " O5'" );
	AtomCOP current_atom ( pose.atom_tree().atom( id::AtomID( atomno,jump_pos) ).get_self_ptr() );
	id::AtomID const parent_id( current_atom->parent()->id() );
	std::string const & parent_name( pose.residue(parent_id.rsd()).atom_name(parent_id.atomno()) );
	if ( parent_name == " P  " ) {
		return true;
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_JumpMover::sample_alternative_chain_connection( pose::Pose & pose, Size const & which_jump ) const
{
	kinematics::FoldTree fold_tree( pose.fold_tree() ); //Make a copy.

	Size const jump_pos1( fold_tree.upstream_jump_residue( which_jump ) );
	Size const jump_pos2( fold_tree.downstream_jump_residue( which_jump ) );

	Size const chain_connection_number = check_in_chain_connections( jump_pos1, jump_pos2 );

	// Is this jump even part of a chain connection?
	if ( chain_connection_number < 1 ) return;

	utility::vector1 < Size > const & res_list1( chain_connections_[chain_connection_number].first );
	utility::vector1 < Size > const & res_list2( chain_connections_[chain_connection_number].second);

	//Now shift around the jump positions. Hmm, this is potentially dangerous.
	Size ntries( 0 );
	Size const MAX_TRIES( 10000 );
	bool success( false );

	// Get ready for a new fold-tree
	Size const num_jumps( fold_tree.num_jump() );
	ObjexxFCL::FArray2D <int> jump_points( 2, num_jumps );
	ObjexxFCL::FArray1D <int> cuts( num_jumps );
	for ( Size n = 1; n <= num_jumps; n++ ) {
		jump_points( 1, n ) = fold_tree.jump_point( 1, n );
		jump_points( 2, n ) = fold_tree.jump_point( 2, n );
		cuts( n ) = fold_tree.cutpoint( n );
	}

	// IS this necessary? Or does which_jump == which_jump_in_list?
	Size which_jump_in_list( 0 );
	for ( Size n = 1; n <= num_jumps; n++ ) {
		if ( Size( jump_points( 1, n ) ) == jump_pos1 &&
				Size( jump_points( 2, n ) ) == jump_pos2 ) {
			which_jump_in_list = n;
			break;
		}
		if ( Size( jump_points( 1, n ) ) == jump_pos2 &&
				Size( jump_points( 2, n ) ) == jump_pos1 ) {
			which_jump_in_list = n;
			break;
		}
	}
	if ( which_jump_in_list == 0 ) {
		utility_exit_with_message( "Problem with fold tree change --> Jump " + I(3, jump_pos1) + " " + I(3, jump_pos2) );
	}

	while ( !success && ntries++ < MAX_TRIES ) {
		Size jump_pos1 = numeric::random::rg().random_element( res_list1 );
		Size jump_pos2 = numeric::random::rg().random_element( res_list2 );
		jump_points( 1, which_jump_in_list ) = std::min( jump_pos1, jump_pos2 );
		jump_points( 2, which_jump_in_list ) = std::max( jump_pos1, jump_pos2 );
		success = fold_tree.tree_from_jumps_and_cuts( pose.total_residue(), num_jumps,
			jump_points, cuts, 1, false /*verbose*/ );
	}

	fill_in_default_jump_atoms( fold_tree, pose );

	// TR << "Changing fold_tree ==> " << fold_tree << std::endl;

	if ( success ) pose.fold_tree( fold_tree );


}

} //movers
} //farna
} //protocols
