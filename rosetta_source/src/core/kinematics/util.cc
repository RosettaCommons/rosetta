// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/util.cc
/// @brief  Kinematics utility functions
/// @author Phil Bradley

// Unit headers
#include <core/kinematics/util.hh>
//#include <core/kinematics/tree/IntraResidueStubJumpAtom.hh>
// Package headers
// AUTO-REMOVED #include <core/id/DOF_ID_Mask.hh>
#include <core/kinematics/constants.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/JumpAtom.hh>
#include <core/kinematics/tree/BondedAtom.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/types.hh>

// Project headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResConnID.hh>

#include <core/id/DOF_ID_Map.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <cassert>


namespace core {
namespace kinematics {

using namespace id;
static numeric::random::RandomGenerator RG(62457); // <- Magic number, do not change it!!!
static basic::Tracer TR( "core.kinematics.util");


///////////////////////////////////////////////////////////////////////////////
// wrapper for the add_atom recursive function which doesn't know anything
// about Residues
//
/// @details recursively called until all atoms in this reside are added.
/// @note The Atoms that are allocated within this function are owned by the atom_ptr
/// which is a container of AtomOP's (see AtomPointer.fwd.hh)
///
/// @note  Returns a raw pointer for internal tree links but memory management is handled by atom_ptr
///


tree::Atom*
add_atom(
	int const atomno,
	int const seqpos,
	Links const & links,
	AtomPointer1D & atom_ptr, // owns all the newly allocated atoms
	bool const add_jump_atom
)
{
	using namespace tree;

	// create new atom
	AtomOP const atom_p( add_jump_atom ? static_cast< Atom* >( new JumpAtom()) : static_cast< Atom* >(new BondedAtom()));

	// fill in the atom_ptr data
	assert( atom_ptr[ atomno ] == 0 );
	atom_ptr[ atomno ] = atom_p;

	// set the atom_id information
	atom_p->id( AtomID( atomno, seqpos ) );
	atom_p->parent( 0 );

	utility::vector1< Size > const & nbrs( links[ atomno ] );
	for ( Size i=1, i_end = nbrs.size(); i<= i_end; ++i ) {
		int const nbr( nbrs[i] );
		if ( atom_ptr[ nbr ] != 0 ) continue;
		atom_p->append_atom( add_atom( nbr, seqpos, links, atom_ptr, false ));
	}

	return atom_p();
}

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/// \brief pick a postion in n_res as the cutpoint
///
///  this is done based on probability info stored in cut_bias_sum.
/// This function is used during fold_tree construction.
int
pick_loopy_cutpoint(
	Size const n_res,
	ObjexxFCL::FArray1D_float const & cut_bias_sum
)
{

	float r = RG.uniform() * cut_bias_sum( n_res );

	int cutpoint( 0 );

	for ( Size i = 1; i <= n_res; ++i ) {
		if ( r > cut_bias_sum(i-1) && r <= cut_bias_sum(i) ) {
			cutpoint = i;
		}
	}

	if ( cutpoint == 0 ) {
		TR.Warning << "pick loopy cutpoint = 0! setting = 1" << std::endl;
		cutpoint = 1;
	}

	return cutpoint;
}



///////////////////////////////////////////////////////////////////////////////
/// helper function for setup_backrub_atom_tree
tree::Atom *
setup_cloned_atom(
									tree::Atom const *, // old_atom,
									utility::vector1< id::AtomID > const & // exclude
)
{
	utility_exit_with_message("needs to be refactored to meet new tree-building guidelines");
	tree::Atom* new_atom( new tree::BondedAtom() );
	/*
	new_atom->parent(0);
	new_atom->id ( old_atom->id () );
	new_atom->xyz( old_atom->xyz() );

	// add old_atom's (non-exclude-) children
	for ( Size i=0; i< Size(old_atom->n_children()); ++i ) { // 0-indexing AAARRGGHHH!!!
		AtomID const & child_id( old_atom->child(i)->id() );
		if ( std::find( exclude.begin(), exclude.end(), child_id ) == exclude.end() ) {
			new_atom->append_atom( old_atom->child(i)->clone(0) ); // recursively copies child's children
		}
	}
	*/
	return new_atom;
}


///////////////////////////////////////////////////////////////////////////////
/// in principal, this should be the inverse of the build_tree function at the beginning
///
// void
// fold_tree_from_atom_tree(
// 	conformation::ResidueCAPs const & residues
// 	AtomTree const & atom_tree,
// 	FoldTree & fold_tree
// )
// {

// 	/// get all inter-residue connections
// 	utility::vector1< Edge

tree::Atom*
setup_backrub_atom_tree(
												utility::vector1< AtomID >, // mainchain, // make our own local copy
												AtomID const &, // downstream_id, // mainchain child of last mainchain atom
												AtomPointer2D const &, // old_atom_pointer,
												utility::vector1< std::pair< Size, Size > > const &, // edges,
												Size const //first_new_pseudo_residue
)
{

	// In the interest of making the gameguys' deadline for refactoring atomtree I'm not thinking about this right now.
	utility_exit_with_message("This code needs to be refactored to use AtomOPs" );
	/*
	using utility::vector1;

	Vector const pseudo_offset( 0,0,0 );
// 	Vector const pseudo_offset( 0.25, 0.25, 0.25 );

	// constraints on the edge list (note that edges are ordered, ie i,j means edge *from* i *to* j (i-->j)
	//
	// -- 1,nbb is the first edge
	//
	// -- for each (i,j) after the first there exists an earlier edge (x,i) (ie, have to have seen first index already)
	//
	// -- no crossing: no pairs (i,j) and (k,l) with i<k<j<l
	//
	// -- for each edge (a,b), a must be downstream of all edges that contain (a,b)
	//    ie if (i,j) contains (a,b) (meaning that the interval [min(i,j),max(i,j)] contains the corresponding ab interval
	//    then there must be a path of edges leading from j to a (j,p1) (p1,p2) (p2,p3) .... (pN,a)
	//    so eg we cant have 1,10 10,7 7,4 1,4 (violated by (1,4) edge: no path from 10 --> 1 and [1,10] contains [1,4])
	//
	// -- no edges (i,j) with |i-j|<=1
	//

	Size const nbb( mainchain.size() );
	mainchain.push_back( downstream_id ); // so that we dont clone this guy
	assert( edges[1].first == 1 && edges[1].second == nbb );

	/// we're going to add pseudo residues, one for each edge in the edges vector
	Size n(0);

	// book-keeping

	// map from each mainchain_atom to the (possibly empty) list of pseudo-residues
	vector1< vector1< Size > > mainchain_to_pseudo( nbb );

	vector1< bool > seen( nbb, false );

	vector1< Atom* > mainchain_atom_pointer( nbb, 0 );
	vector1< Atom* > pseudo_atom_pointer( edges.size(), 0 );

	// create the root of the tree
	AtomOP root( setup_cloned_atom( old_atom_pointer[ mainchain[1] ], mainchain ) );
	mainchain_atom_pointer[ 1 ] = root;
	seen[ 1 ] = true;

	while ( n < edges.size() ) {
		++n;
		std::pair< Size, Size > const & edge( edges[n] );

		Size const a( edge.first );
		Size const b( edge.second );
		assert( seen[a] );

		// what should our anchor atom be?
		// look at all pseudo rsds associated to a
		Atom * anchor_atom;
		Size anchor(0);
		vector1< Size > const & pseudo_a( mainchain_to_pseudo[a] );
		for ( Size i=1; i<= pseudo_a.size(); ++i ) {
			Size const nn( pseudo_a[i] );
			assert( edges[nn].second == a );
			if ( ( edges[nn].first < b && b < a ) || ( edges[nn].first > b && b > a ) ) {
				anchor = nn; // dont break -- want the smallest segment containing a,b
			}
		}
		if ( anchor ) {
			anchor_atom = pseudo_atom_pointer[ anchor ];
		} else {
			// connect to the authentic a
			anchor_atom = mainchain_atom_pointer[ a ];
		}

		// has b been seen before?
		Atom const * old_b_atom( old_atom_pointer[ mainchain[b] ] );
		if ( !seen[ b ] ) {
			seen[b] = true;
			// add b and b's non-mainchain children
			Atom* b_atom( setup_cloned_atom( old_b_atom, mainchain ) );
			anchor_atom->insert_atom( b_atom ); // NOTE: insert atom
			mainchain_atom_pointer[ b ] = b_atom;
		}


		// add a new pseudo rsd at b's position
		mainchain_to_pseudo[ b ].push_back( n );
		Size const pseudo_b_seqpos( first_new_pseudo_residue + n - 1 );
		AtomID const pseudo_b_id( AtomID( 1, pseudo_b_seqpos ) );
		// delete the old one
		Atom* old_pseudo_b_atom( old_atom_pointer[ pseudo_b_id ] );
		old_pseudo_b_atom->parent()->delete_atom( old_pseudo_b_atom );
		delete old_pseudo_b_atom;
		// create a new one
		Atom* pseudo_b_atom( new BondedAtom() );
		pseudo_b_atom->id( pseudo_b_id );
		pseudo_b_atom->xyz( old_b_atom->xyz() + pseudo_offset );
		anchor_atom->append_atom( pseudo_b_atom );
		pseudo_atom_pointer[ n ] = pseudo_b_atom;
		std::cout << "new pseudo! " << n << ' ' << a << ' ' << b << ' ' << mainchain[a] << ' ' << mainchain[b] << ' ' <<
			pseudo_b_seqpos << ' ' << old_b_atom->id() << ' ' <<
			old_b_atom->xyz()[0] << ' ' <<
			old_b_atom->xyz()[1] << ' ' <<
			old_b_atom->xyz()[2] << std::endl;

		// check if this edge is terminal, if so add all intervening mainchain atoms and their children
		{
			bool terminal( true );
			for ( Size n2=n+1; n2<= edges.size(); ++n2 ) {
				if ( ( edges[n2].first == b ) &&
						 ( a < b && edges[n2].second < b || a > b && edges[n2].second > b ) ) {
					terminal = false;
				}
			}

			if ( terminal ) {
				int const dir( b<a ? 1 : -1 );
				Atom * parent_atom( pseudo_b_atom );
				for ( int c=b+dir; c != (int)a; c += dir ) {
					assert( !seen[c] );

					// add c and c's non-mainchain children
					Atom const * old_c_atom( old_atom_pointer[ mainchain[c] ] );
					Atom* c_atom( setup_cloned_atom( old_c_atom, mainchain ) );
					parent_atom->insert_atom( c_atom ); // at front of list since this is mainchain. may have already added kids
					mainchain_atom_pointer[ c ] = c_atom;
					seen[c] = true;

					parent_atom = c_atom;

				} // walk from b->a adding the mainchain atoms and their children to the tree

			} // if terminal

		} // scope


	} // loop over edges, add one pseudo rsd for each edge

	// confirm that all mainchain atoms have been seen
	for ( Size i=1; i<= nbb; ++i ) {
		assert( seen[i] );
	}
	return root;
	*/
	return 0;

}

/// @brief prints something like this ***1***C***1*********2***C********3****C****2********3*****
void
simple_visualize_fold_tree( FoldTree const & fold_tree, std::ostream& out ) {
	for ( Size pos = 1; pos <= fold_tree.nres(); pos++ ) {
		bool special( false );
		if ( fold_tree.is_jump_point( pos ) ) {
			for ( Size jnr=1; jnr<=fold_tree.num_jump(); jnr++ ) {
				if ( (Size) fold_tree.jump_edge( jnr ).start() == pos || (Size)fold_tree.jump_edge( jnr ).stop() == pos ) {
					if ( special ) out << "/";
					out << jnr;
					special = true;
				}
			}
		}
		if ( fold_tree.is_cutpoint( pos ) ) {
			if ( special ) out << "/";
			out << "C";
			special = true;
		}
		if (!special ) out << "*";
	}
	out << std::endl;
}


/// @brief prints something like this ***1***C***1*********2***C********3****C****2********3*****
///                                   **********xxxxxxxxxxxxx************************************
void
simple_visualize_fold_tree_and_movemap( FoldTree const & fold_tree, MoveMap const& mm, std::ostream& out ) {
	std::string move;
	out << "\n";
	for ( Size pos = 1; pos <= fold_tree.nres(); pos++ ) {
		bool special( false );
		if ( fold_tree.is_jump_point( pos ) ) {
			for ( Size jnr=1; jnr<=fold_tree.num_jump(); jnr++ ) {
				if ( (Size) fold_tree.jump_edge( jnr ).start() == pos || (Size)fold_tree.jump_edge( jnr ).stop() == pos ) {
					if ( special ) {
						out << "/";
						move.push_back( '.' );
					}
					out << jnr;
					move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
					special = true;
				}
			}
		}
		if ( fold_tree.is_cutpoint( pos ) ) {
			if ( special ) {
				out << "/";
				move.push_back( '.' );
			}
			out << "C";
			move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
			special = true;
		}
		if (!special ) {
			out << "*";
			move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
		}
	}
	out << "\n" << move;
	out << std::endl;
}

/// @brief prints something like this ***1***C***1*********2***C********3****C****2********3*****
///                                   **********xxxxxxxxxxxxx************************************
void
simple_visualize_fold_tree_and_movemap_bb_chi( FoldTree const & fold_tree, MoveMap const& mm, std::ostream& out ) {
	std::string move;
	std::string move_chi;
	out << "\n";
	for ( Size pos = 1; pos <= fold_tree.nres(); pos++ ) {
		bool special( false );
		if ( fold_tree.is_jump_point( pos ) ) {
			for ( Size jnr=1; jnr<=fold_tree.num_jump(); jnr++ ) {
				if ( (Size) fold_tree.jump_edge( jnr ).start() == pos || (Size)fold_tree.jump_edge( jnr ).stop() == pos ) {
					if ( special ) {
						out << "/";
						move.push_back( '.' );
					}
					out << jnr;
					move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
					move_chi.push_back( mm.get_chi( pos ) ? '*' : 'x' );
					special = true;
				}
			}
		}
		if ( fold_tree.is_cutpoint( pos ) ) {
			if ( special ) {
				out << "/";
				move.push_back( '.' );
				move.push_back( '.' );
			}
			out << "C";
			move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
			move_chi.push_back( mm.get_chi( pos ) ? '*' : 'x' );
			special = true;
		}
		if (!special ) {
			out << "*";
			move_chi.push_back( mm.get_chi( pos ) ? '*' : 'x' );
			move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
		}
	}
	out << "\n" << move << "\n" << move_chi;
	out << std::endl;
}


///@brief linearizes (or defoliates, if you prefer) a FoldTree.  "default" FoldTrees produced by the PDB reader have all chains (peptide edges) starting from jumps relative to residue 1.  This code modifies the tree to instead have all the jumps be relative to the preceding edge.  It is not tested with ligands and will not work with "functional" jumps.  From A to B:
///A:FOLD_TREE  EDGE 1 78 -1  EDGE 1 79 1   EDGE 79 454 -1  EDGE 1 455 2    EDGE 455 540 -1  EDGE 1 541 3    EDGE 541 697 -1
///B:FOLD_TREE  EDGE 1 78 -1  EDGE 78 79 1  EDGE 79 454 -1  EDGE 454 455 2  EDGE 455 540 -1  EDGE 540 541 3  EDGE 541 697 -1
core::kinematics::FoldTree
linearize_fold_tree( core::kinematics::FoldTree const & tree ) {
	core::kinematics::FoldTree newtree;
	for( core::kinematics::FoldTree::const_iterator it(tree.begin()), end(tree.end()); it != end; ++it){
		//if it is not a jump, we don't modify it
		if( !it->is_jump() ) newtree.add_edge(*it);
		//if it is a jump, we move start() to stop-1.  This is naive but works for the intended case.
		else newtree.add_edge(core::kinematics::Edge(it->stop()-1, it->stop(), it->label()));
	}
	return newtree;
}


} // namespace kinematics
} // namespace core
