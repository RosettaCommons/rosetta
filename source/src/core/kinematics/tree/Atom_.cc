// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/tree/Atom_.cc
/// @brief  Kinematics Atom abstract base class
/// @author Phil Bradley


// Unit headers
#include <core/kinematics/tree/Atom_.hh>

// Package headers
#include <core/kinematics/AtomWithDOFChange.hh>
#include <core/kinematics/ResidueCoordinateChangeList.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.hh>

// Utility headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <iostream>

#include <core/id/AtomID_Map.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/types.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.io.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector0.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace kinematics {
namespace tree {

static basic::Tracer TR( "core.kinematics.tree.Atom_" );

/////////////////////////////////////////////////////////////////////////////
/// @details get the input stub for building this atom first
void
Atom_::update_xyz_coords()
{
	Stub stub( get_input_stub() );
	update_xyz_coords( stub );
}


/////////////////////////////////////////////////////////////////////////////
/// @details get the input stub for building this atom first
void
Atom_::update_internal_coords(
	bool const recursive
)
{
	Stub stub( get_input_stub() );
	update_internal_coords( stub, recursive );
}


/////////////////////////////////////////////////////////////////////////////
/// @details first this atom, then its parent and then recursively all its children
void
Atom_::show() const
{
	std::cout << "ATOM: " << atom_id_ << std::endl;
	std::cout << "   PARENT: ";
	AtomCOP parent_op = parent();
	if ( parent_op ) {
		std::cout << parent_op->id() << std::endl;
		std::cout << "  GRAND PARENT: ";
		AtomCOP grand_parent_op = parent_op->parent();
		if ( grand_parent_op ) {
			std::cout << grand_parent_op->id() << std::endl;
		} else {
			std::cout << " NULL" << std::endl;
		}

	} else {
		std::cout << " NULL" << std::endl;
	}

	std::cout << "   CHILDREN: ";
	for ( auto const & atom : atoms_ ) {
		std::cout << atom->id() << ' ';
	}
	std::cout << std::endl;
	for ( auto const & atom : atoms_ ) {
		atom->show();
	}
}

/////////////////////////////////////////////////////////////////////////////
/// @details first this atom, then its parent and then recursively all its children up to n_level
void
Atom_::show(int const & n_level) const
{
	using namespace ObjexxFCL::format;
	TR << "ATOM: " << atom_id_ << std::endl;
	TR << "POSITION: " << F(8,3,x()) << F(8,3,y()) << F(8,3,z()) << std::endl;
	TR << "STUB: " << get_stub()  << std::endl;
	TR << "INPUT_STUB: " << get_input_stub()  << std::endl;
	TR << "   PARENT: ";
	AtomCOP parent_op = parent();
	if ( parent_op ) {
		TR << parent_op->id() << std::endl;
		TR << "  GRAND PARENT: ";
		AtomCOP grand_parent_op = parent_op->parent();
		if ( grand_parent_op ) {
			TR << grand_parent_op->id() << std::endl;
		} else {
			TR << " NULL" << std::endl;
		}

	} else {
		TR << " NULL" << std::endl;
	}

	TR << "   CHILDREN: ";
	for ( auto const & atom : atoms_ ) {
		TR << atom->id() << ' ';
	}
	TR << std::endl;

	//TR << "Yifan debug: " << n_level << std::endl;
	if ( n_level > 0 ) {
		int const next_level(n_level - 1);
		for ( auto const & atom : atoms_ ) {
			atom->show(next_level);
		}
	}
}


/////////////////////////////////////////////////////////////////////////////
//
// this should probably be debugged more thoroughly
//
/// @details update domain map for this atom and all its offspring
/// consider this like a graph coloring problem. we are recursively (depth-first)
/// assigning a color to each atom.
void
Atom_::update_domain_map(
	int & current_color,
	int & biggest_color,
	DomainMap & domain_map,
	AtomID_Mask const & dof_moved,
	AtomID_Mask const & atom_moved
) const
{
	int color( current_color );

	bool const my_dof_moved(  dof_moved[ atom_id_ ] );
	bool const my_xyz_moved( atom_moved[ atom_id_ ] );

	if ( my_dof_moved ) {
		++biggest_color;
		if ( !is_jump() ) {
			// propagates through atoms after me in the list of my parent's childrn
			// that's assuming its the phi-torsion that's changed
			// note that current_color is passed by reference
			current_color = biggest_color;
			// but if my D or THETA dof's are changing then myself
			// and my children will be on a different rigid-body from
			// my younger siblings
			++biggest_color;
		}
		color = biggest_color;
	}

	// update domain_map
	int const current_map( domain_map( atom_id_.rsd() ) );
	if ( !my_dof_moved && my_xyz_moved ) {
		// no propagating change, but my xyz coords have changed
		domain_map( atom_id_.rsd() ) = 0;
	} else if ( current_map < 0 ) {
		// unassigned
		domain_map( atom_id_.rsd() ) = color;
	} else if ( current_map == 0 || current_map == color ) {
		// leave the same
	} else {
		// color change within a residue
		domain_map( atom_id_.rsd() ) = 0;
	}


	for ( auto const & atom : atoms_ ) {
		atom->update_domain_map( color, biggest_color, domain_map, dof_moved, atom_moved );
	}
}


/////////////////////////////////////////////////////////////////////////////
/// @details if the child atom is a jump atom, put it before all the non-jump children
/// atoms but after all the children jump atoms. Otherwise, put it at the end
/// of the children atom list.
void
Atom_::append_atom(
	AtomOP const atom
)
{
	atom->parent( get_self_weak_ptr() ); // Note: You cannot give "this" to the child atom, or you will loose the reference-count information using boost::weak_ptr
	if ( atom->is_jump() ) {
		atoms_.insert( nonjump_atoms_begin(), atom );
	} else {
		atoms_.push_back(atom);
	}
}

/////////////////////////////////////////////////////////////////////////////
/// @details only unlink the child atom from this atom (parent). No recursive operation or
/// freeing memory( use Atom_::erase() to free up memory for an atom and all its
/// children
void
Atom_::delete_atom(
	AtomOP const child
)
{
	Atoms::iterator const iter( std::find( atoms_.begin(), atoms_.end(), child ) );
	if ( iter == atoms_.end() ) {
		std::cerr << "child not present in atoms list! " << atoms_.size() << std::endl;
		utility_exit();
	}
	atoms_.erase( iter );
}

/////////////////////////////////////////////////////////////////////////////
/// @details if the child atom is a jump atom, put it before all jump children atoms;
/// otherwise, put it before all non-jump children atoms; Different from
/// Atom_::append_atom in that it put the new child atom before instead after
/// the other existing child atoms.
void
Atom_::insert_atom(
	AtomOP const atom
)
{
	atom->parent( get_self_weak_ptr() ); // Note: you cannot give "this" to the child, or the reference-count information will be lost with boost::weak_ptr.
	if ( atom->is_jump() ) {
		atoms_.insert( atoms_.begin(), atom );
	} else {
		atoms_.insert( nonjump_atoms_begin(), atom );
	}
}

/////////////////////////////////////////////////////////////////////////////
/// @details insert the child atom in a specified position. If the specified postion is
/// out of range, put it either at the beginning or the end of the child atom list.
/// note that jump child atoms are always listed before non-jump atoms and the
/// index must be a non-negative integer.
void
Atom_::insert_atom(
	AtomOP const atom,
	int const index
)
{
	atom->parent( get_self_weak_ptr() ); // Note: you cannot give "this" to the child, or the reference-count information will be lost with boost::weak_ptr

	debug_assert( index >= 0 );
	auto pos( atoms_.begin() + index );
	if ( atom->is_jump() ) {
		if ( pos > nonjump_atoms_begin() ) pos = nonjump_atoms_begin();
	} else {
		if ( pos < nonjump_atoms_begin() ) {
			pos = nonjump_atoms_begin();
			//if ( keep_1st_child_pos() ) ++pos;
		}
		if ( pos > atoms_.end() ) pos = atoms_.end();
	}


	atoms_.insert( pos, atom );
}

/////////////////////////////////////////////////////////////////////////////
/// @details old atom and new atom need to belong to the same type ( either both jump
/// atoms or both non-jump atoms. New atom is inserted at the position of old
/// atom.
void
Atom_::replace_atom(
	AtomOP const old_atom,
	AtomOP const new_atom
)
{
	debug_assert( ( old_atom->is_jump() && new_atom->is_jump() ) ||
		( !old_atom->is_jump() && !new_atom->is_jump() ) );

	auto iter( std::find( atoms_.begin(), atoms_.end(), old_atom ) );
	if ( iter == atoms_.end() ) {
		TR.Fatal << "old_atom not present in atoms list! " << atoms_.size() << std::endl;
		utility_exit_with_message("'Replacing' an atom which doesn't currently exist.");
	}
	new_atom->parent( get_self_weak_ptr() ); // Note: you cannot give "this" to the child atom, or the reference-count data will be lost with boost::weak_ptr
	iter = atoms_.insert( iter, new_atom );
	//   std::cout << (*iter == new_atom) << ' ' <<
	//    (*iter == old_atom) << std::endl;
	++iter;
	debug_assert( *iter == old_atom );
	atoms_.erase( iter );

	new_atom->steal_inversion(old_atom);  //fpd steal the iversion state when replacing this atom
}


/////////////////////////////////////////////////////////////////////////////
/// @details returns 0 if atom doesnt exist
AtomCOP
Atom_::get_nonjump_atom(
	Size const index
) const
{
	auto iter( nonjump_atoms_begin() );
	iter += index;
	if ( iter >= atoms_.end() ) {
		return nullptr;
	} else {
		return *iter;
	}
}


/////////////////////////////////////////////////////////////////////////////
Size
Atom_::n_children() const
{
	return atoms_.size();
}

/////////////////////////////////////////////////////////////////////////////
AtomCOP
Atom_::child( Size const k ) const
{
	return atoms_[k];
}


/////////////////////////////////////////////////////////////////////////////
AtomOP
Atom_::child( Size const k )
{
	return atoms_[k];
}

/////////////////////////////////////////////////////////////////////////////
/// @details the atom-index of this child
Size
Atom_::child_index( AtomCOP child ) const
{
	debug_assert( child->parent().get() == this );
	for ( Size k=0; k< atoms_.size(); ++k ) {
		if ( atoms_[k] == child ) return k;
	}
	utility_exit_with_message( "problemo in Atom_'s atom list" );
	return Size(-1);
}

Size
Atom_::raw_child_index( Atom const * child ) const
{
	debug_assert( child->raw_parent() == this );
	for ( Size k=0; k < atoms_.size(); ++k ) {
		if ( atoms_[k].get() == child ) return k;
	}
	utility_exit_with_message( "problemo in Atom_'s atom list" );
	return Size(-1);
}

/////////////////////////////////////////////////////////////////////////////
/// @details the improper dihedral from child1 to child2 about my x-axis
/// (ie, axis defined by me and my parent for nonjump atoms). Since phi_ for
/// a non-first branched atom is defined as the improper angle offset with
/// respect to its previous sibling, we just need to add up all the offsets
/// between them.
Real
Atom_::dihedral_between_bonded_children(
	Atom const & child1,
	Atom const & child2
) const
{
	// debug args
	if ( child1.raw_parent() != this || child2.raw_parent() != this ) {
		utility_exit_with_message("Atom_::dihedral_between_bonded_children: atoms are not both my children!");
	}
	if ( child1.is_jump() || child2.is_jump() ) {
		utility_exit_with_message("Atom_::dihedral_between_bonded_children: one of the atoms is a JumpAtom!");
	}

	// keep track of which atoms we've seen and in what order
	Atom const * first_atom( nullptr );
	Atom const * second_atom( nullptr );

	Real phi_offset(0.0);

	for ( auto it=nonjump_atoms_begin(),
			it_end=atoms_end(); it != it_end; ++it ) {
		if ( first_atom ) phi_offset += (*it)->dof(PHI);

		if ( (*it).get() == &child1 || (*it).get() == &child2 ) {
			if ( first_atom ) {
				second_atom = (*it).get(); // seen both
				break;
			} else {
				first_atom = (*it).get();
			}
		}
	} // loop over nonjump atoms

	if ( !second_atom ) {
		utility_exit_with_message("Atom_::dihedral_between_bonded_children: atoms not found!");
	}
	if ( second_atom == &child1 ) phi_offset *= -1.0;
	return phi_offset;
}


/////////////////////////////////////////////////////////////////////////////
// I don't know where this routine is used...
//
// I don't think it even works in trunk.
//
bool
Atom_::downstream( AtomCOP atom1 ) const
{
	for ( int ii=0, ie= n_children(); ii < ie; ++ii ) {
		if ( child(ii) == atom1 ) {
			return true;
		} else {
			if ( child(ii)->downstream( atom1 ) ) return true;
		}
	}
	return false;
}


/////////////////////////////////////////////////////////////////////////////
/// @details my stub is center at myself. Normally for bonded atom, X direction is
/// from my parent to me; Z direction is perpendicular to the plane defined
/// by myself, my parent and my parent's parent
Stub
Atom_::get_stub() const
{
	try {
		return Stub(
			position(),
			stub_atom1()->position(),
			stub_atom2()->position(),
			stub_atom3()->position()
		);
	} catch ( utility::excn::Exception const & ) {
		TR.Error << "Issue getting stub for atom " << atom_id() << " -- possibly due to degenerate/colinear atoms:" << std::endl;
		TR.Error << "\t " << input_stub_atom1_id() << " -- " << input_stub_atom1()->position() << std::endl;
		TR.Error << "\t " << input_stub_atom2_id() << " -- " << input_stub_atom2()->position() << std::endl;
		TR.Error << "\t " << input_stub_atom3_id() << " -- " << input_stub_atom3()->position() << std::endl;
		throw; // Make sure to re-throw error. This is not recoverable.
	}
}


/////////////////////////////////////////////////////////////////////////////
/// @details the stub that is passed to me during folding, which is normally the stub
/// centered at the parent atom.
Stub
Atom_::get_input_stub() const
{
	if ( !parent_.expired() ) {
		//   std::cout << "Get input stub: ";
		//   std::cout << "(0: " << input_stub_atom0()->atom_id().rsd() << ", "<< input_stub_atom0()->atom_id().atomno() << ") ";
		//   std::cout << "(1: " << input_stub_atom1()->atom_id().rsd() << ", "<< input_stub_atom1()->atom_id().atomno() << ") ";
		//   std::cout << "(2: " << input_stub_atom2()->atom_id().rsd() << ", "<< input_stub_atom2()->atom_id().atomno() << ") ";
		//   std::cout << "(3: " << input_stub_atom3()->atom_id().rsd() << ", "<< input_stub_atom3()->atom_id().atomno() << ") " << std::endl;
		try {
			return Stub(
				input_stub_atom0()->position(),
				input_stub_atom1()->position(),
				input_stub_atom2()->position(),
				input_stub_atom3()->position()
			);
		} catch ( utility::excn::Exception const & ) {
			TR.Error << "Issue getting stub for atom " << atom_id() << " -- possibly due to degenerate/colinear atoms:" << std::endl;
			TR.Error << "\t " << input_stub_atom1_id() << " -- " << input_stub_atom1()->position() << std::endl;
			TR.Error << "\t " << input_stub_atom2_id() << " -- " << input_stub_atom2()->position() << std::endl;
			TR.Error << "\t " << input_stub_atom3_id() << " -- " << input_stub_atom3()->position() << std::endl;
			throw; // Make sure to re-throw error. This is not recoverable.
		}
	} else {
		return default_stub;
	}
}


/////////////////////////////////////////////////////////////////////////////
/// @details call parent's previous_child method to get its previous sibling; return 0
/// if no parent is present.
AtomCOP
Atom_::previous_sibling() const
{
	AtomCOP parent_op = parent();
	if ( parent_op ) {
		return parent_op->previous_child( get_self_ptr() );
	} else {
		return nullptr;
	}
}

/////////////////////////////////////////////////////////////////////////////
/// @details return 0 if the input child is the first child in the list
AtomCOP
Atom_::previous_child(
	AtomCOP child
) const
{
	//std::cout << "atoms_.size() = " << atoms_.size() << std::endl;
	auto iter( std::find( atoms_.begin(), atoms_.end(), child ) );
	if ( iter == atoms_.end() ) {
		std::cerr << "child not present in atoms list! " << atoms_.size() << std::endl;
		utility_exit();
	}
	if ( iter == atoms_.begin() ) {
		return nullptr;
	} else {
		--iter;
		return *iter;
	}
}

/////////////////////////////////////////////////////////////////////////////
/// @details return 0 if the input child is the last child in the list
AtomOP
Atom_::next_child(
	AtomCOP child
)
{
	//std::cout << "atoms_.size() = " << atoms_.size() << std::endl;
	Atoms::const_iterator iter( std::find( atoms_.begin(), atoms_.end(), child ) );
	if ( iter == atoms_.end() ) {
		std::cerr << "child not present in atoms list! " << atoms_.size() << std::endl;
		utility_exit();
	}
	++iter;
	if ( iter == atoms_.end() ) {
		return nullptr;
	} else {
		return *iter;
	}
}


/////////////////////////////////////////////////////////////////////////////
/// @details
/// - bonded-atom always has its stub defined (from its parent)
/// - jump-atom needs to have at least three atoms (including itself) at its
///   end to define a stub, such as three linearly bonded atoms or two children
///   atoms bonded to itself.
bool
Atom_::stub_defined() const
{
	// have to handle a couple of cases here:

	// note -- in counting dependent atoms, exclude JumpAtom's


	// 1. no dependent atoms --> no way to define new coord sys
	//    on this end. ergo take parent's M and my xyz
	//
	// 2. one dependent atom --> no way to define unique coord
	//    on this end, still take parent's M and my xyz
	//
	// 3. two or more dependent atoms
	//    a) if my first atom has a dependent atom, use
	//       myself, my first atom, and his first atom
	//
	//    b) otherwise, use
	//       myself, my first atom, my second atom


	if ( is_jump() ) {
		AtomCOP first = get_nonjump_atom(0);
		if ( first != nullptr &&
				( first->get_nonjump_atom(0) != nullptr || get_nonjump_atom(1) != nullptr ) ) {
			return true;
		} else {
			return false;
		}
	} else {
		return true;
	}
}


/////////////////////////////////////////////////////////////////////////////
/// @details the coordinates of the atom "*child" -- one of my children -- have
/// changed. This routine updates the torsions to reflect this. Useful
/// if we have just repacked or rotamer-trialed, ie sidechain atoms_
/// have moved but the backbone is still the same, more efficient than
/// calling update_internal_coords on the entire tree...
///
/// @note update_internal_coords is called recursively on all of *child's
/// children.
void
Atom_::update_child_torsions(
	AtomOP const child
)
{
	Stub my_stub( this->get_stub() );
	bool found( false );
	for ( auto it = atoms_.begin(), it_end = atoms_.end(); it != it_end; ++it ) {
		if ( *it == child ) {
			(*it)->update_internal_coords( my_stub );
			// phi torsion for the atom after child may have changed
			++it;
			if ( it != it_end ) {
				(*it)->update_internal_coords( my_stub, false /* not recursive */ );
			}
			found = true;
			break;
		} else {
			// just advances the stub as if we had called update_internal_coords
			(*it)->update_stub( my_stub );
		}
	}
	if ( !found ) {
		std::cerr << "update_child_torsions:: child not in atoms list" <<
			atom_id_ << ' ' << child->id() << std::endl;
		utility_exit();
	}
}


/////////////////////////////////////////////////////////////////////////////
/// @details Keep track of any residue that moves so that the Conformation object
/// may be correctly updated.
void
Atom_::transform_Ax_plus_b_recursive(
	Matrix const & A,
	Vector const & b,
	ResidueCoordinateChangeList & res_change_list
)
{
	position_ = A * position_ + b;
	debug_assert( position_.is_finite() );

	for ( auto & atom : atoms_ ) {
		atom->transform_Ax_plus_b_recursive( A, b, res_change_list );
	}
	res_change_list.mark_residue_moved( atom_id_ );
}


/////////////////////////////////////////////////////////////////////////////
void
Atom_::get_path_from_root( utility::vector1< AtomCAP > & path ) const
{
	AtomCOP parent = parent_.lock();
	if ( parent ) {
		parent->get_path_from_root( path );
	}
	path.push_back( get_self_weak_ptr() );
}


/////////////////////////////////////////////////////////////////////////////
bool
Atom_::atom_is_on_path_from_root( AtomCOP atm ) const
{
	if ( atm.get() == this ) {
		return true;
	}
	AtomCOP parent = parent_.lock();
	return ( parent && parent->atom_is_on_path_from_root( atm ) );
}


/////////////////////////////////////////////////////////////////////////////
Atom_::Atoms_ConstIterator
Atom_::nonjump_atoms_begin() const
{
	auto iter( atoms_.begin() );
	while ( iter != atoms_.end() && (*iter)->is_jump() ) ++iter;
	return iter;
}


/////////////////////////////////////////////////////////////////////////////
Atom_::Atoms_Iterator
Atom_::nonjump_atoms_begin()
{
	auto iter( atoms_.begin() );
	while ( iter != atoms_.end() && (*iter)->is_jump() ) ++iter;
	return iter;
}

/// @details the "start atom index" is used to determine which subtrees have already
/// been examined.  The atom tree guarantees that the dfs's occur in increasing order
/// by dof refold indices.  If this atom has a dof refold index less than the start
/// atom index, then a dfs has previously been launched from this node and the recursion
/// must stop: subtrees must be visited at most once, or running time grows quadratically
/// in the number of atoms in the tree.  The atom is responsible for handing the
/// start atom index down to its children in the recursion (if it continues recursing).
void
Atom_::dfs(
	AtomDOFChangeSet & changeset,
	ResidueCoordinateChangeList & res_change_list,
	Size const start_atom_index
) const
{
	res_change_list.mark_residue_moved( atom_id_ );
	if ( start_atom_index == dof_refold_index_ ) {
		/// This is the root atom for a new subtree dfs.
		for ( auto const & atom : atoms_ ) {
			atom->dfs( changeset, res_change_list, start_atom_index );
		}

	} else if ( dof_refold_index_ != 0 ) {
		changeset[ dof_refold_index_ ].reached_ = true;
		if ( dof_refold_index_ > start_atom_index ) {
			//recurse -- the subtree rooted from here won't be reached otherwise.
			for ( auto const & atom : atoms_ ) {
				atom->dfs( changeset, res_change_list, start_atom_index );
			}
		}
	} else { // dof_refold_index() == 0 -- recurse.
		for ( auto iter = atoms_begin(), iter_e = atoms_end(); iter != iter_e; ++iter ) {
			(*iter)->dfs( changeset, res_change_list, start_atom_index );
		}
	}
}

/// @details Records this atom as an atom with one-or-more changed dofs in the input
/// AtomDOFChangeSet if this is first DOF on this atom to change; if it is not the
/// first DOF on this atom that has changed, then this atom is already in set and
/// this atom is recorded at position dof_refold_index_.
void
Atom_::note_dof_change(
	AtomDOFChangeSet & changset
)
{
	if ( dof_refold_index_ == 0 ) {
		/// This is the first dof on this atom to change.  Add it to the list of atoms that
		/// have had DOFs change, and record the position of this atom in that list.

		changset.push_back( AtomWithDOFChange( atom_id_ ) );
		dof_refold_index_ = changset.size();
	} else {
		debug_assert( changset[ dof_refold_index_].atomid_ == atom_id_ );
	}
}


void
Atom_::abort_bad_call() const
{
	std::cerr << "kinematics::Atom bad method call in Atom hierarchy!" << std::endl;
	utility_exit();
}

Atom const *
Atom_::raw_parent() const
{
	return raw_parent_;
}

Atom const *
Atom_::raw_previous_sibling() const
{
	if ( raw_parent_ ) {
		return raw_parent_->raw_previous_child( this );
	} else {
		return nullptr;
	}
}

Atom const *
Atom_::raw_previous_child(
	Atom const * child
) const
{
	for ( Size ii = 0; ii < atoms_.size(); ++ii ) {
		if ( atoms_[ ii ].get() == child ) {
			if ( ii == 0 ) {
				return nullptr;
			} else {
				return atoms_[ ii-1 ].get();
			}
		}
	}
	std::cerr << "child not present in atoms list! " << atoms_.size() << std::endl;
	utility_exit();
	return nullptr;
}


Atom const *
Atom_::raw_input_stub_atom0() const {
	return raw_parent();
}


Atom const *
Atom_::raw_input_stub_atom1() const
{
	return raw_parent()->raw_stub_atom1();
}


Atom const *
Atom_::raw_input_stub_atom2() const {
	return raw_parent()->raw_stub_atom2();
}


Atom const *
Atom_::raw_input_stub_atom3() const
{
	Atom const * parent_ptr = raw_parent();
	Atom const * sibling_ptr = raw_previous_sibling();
	if ( is_jump() || ! sibling_ptr || sibling_ptr->is_jump() ||
			is_collinear( *(raw_parent()->raw_stub_atom1()), *(raw_parent()->raw_stub_atom2()), *sibling_ptr) ||
			( parent_ptr->is_jump() && sibling_ptr->id() == parent_ptr->stub_atom2_id() ) ) {
		return parent_ptr->raw_stub_atom3();
	} else {
		return sibling_ptr;
	}
}

Atom const *
Atom_::raw_get_nonjump_atom(
	Size const i
) const
{
	auto iter( nonjump_atoms_begin() );
	iter += i;
	if ( iter >= atoms_.end() ) {
		return nullptr;
	} else {
		return iter->get();
	}
}


}
} // namespace kinematics
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::kinematics::tree::Atom_::save( Archive & arc ) const {
	arc( CEREAL_NVP( atom_id_ ) ); // AtomID
	arc( CEREAL_NVP( parent_ ) ); // AtomAP
	// don't serialize the raw pointer to the parent
	// EXEMPT raw_parent_
	arc( CEREAL_NVP( position_ ) ); // PointPosition
	arc( CEREAL_NVP( atoms_ ) ); // Atoms
	arc( CEREAL_NVP( dof_refold_index_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::kinematics::tree::Atom_::load( Archive & arc ) {

	arc( atom_id_ ); // AtomID
	arc( parent_ ); // AtomAP
	AtomOP parent( parent_.lock() );
	if ( parent ) {
		raw_parent_ = parent.get();
	}
	arc( position_ ); // PointPosition
	arc( atoms_ ); // Atoms
	arc( dof_refold_index_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::kinematics::tree::Atom_ );
CEREAL_REGISTER_TYPE( core::kinematics::tree::Atom_ )

CEREAL_REGISTER_DYNAMIC_INIT( core_kinematics_tree_Atom_ )
#endif // SERIALIZATION
