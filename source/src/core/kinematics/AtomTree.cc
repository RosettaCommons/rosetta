// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/AtomTree.cc
/// @brief  Atom tree class
/// @author Phil Bradley


// Unit headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/ResidueCoordinateChangeList.hh>

// Package headers
#include <core/kinematics/AtomWithDOFChange.hh>
#include <core/kinematics/types.hh>
#include <core/kinematics/tree/Atom.hh>

// Basic headers
#include <basic/basic.hh> // periodic_range
#include <basic/prof.hh> // profiling
#include <basic/Tracer.hh> // profiling

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1D.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/assert.hh>
#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

#include <core/id/AtomID_Map.srlz.hh>

// Cereal headers
#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace kinematics {

static THREAD_LOCAL basic::Tracer TR( "core.kinematics.AtomTree" );

/////////////////////////////////////////////////////////////////////////////
/// @details this will claim the tree as our own. new_root has information about its children,
/// and they have information about their children. From those atom positions, internal
/// coordinates can be updated and atom pointers will be added into the AtomTree map.
/// @note that we steal the atoms, ie it's incorporated into the AtomTree (by recording
/// their pointers),not cloned
AtomTree::AtomTree(
	AtomPointer2D const & new_atom_pointer,
	bool const from_xyz // = true
):
	this_weak_ptr_(/* 0 */),
	root_( /* 0 */ ),
	atom_pointer_(), // default_setting_ = null pointer
	internal_coords_need_updating_( false ),
	xyz_coords_need_updating_( false ),
	topological_match_to_( /* 0 */ ),
	external_coordinate_residues_changed_( ResidueCoordinateChangeListOP( new ResidueCoordinateChangeList ) )
{
	replace_tree( new_atom_pointer, from_xyz );
	external_coordinate_residues_changed_->total_residue( new_atom_pointer.size() );
}

AtomTree::AtomTree():
	this_weak_ptr_(/* 0 */),
	root_(/* 0 */),
	atom_pointer_(),
	internal_coords_need_updating_( false ),
	xyz_coords_need_updating_( false ),
	topological_match_to_( /* 0 */ ),
	external_coordinate_residues_changed_( ResidueCoordinateChangeListOP( new ResidueCoordinateChangeList ) )
{}

/// @brief Destructor
AtomTree::~AtomTree()
{
	clear();
}


/////////////////////////////////////////////////////////////////////////////
/// @details copy ctor, uses operator=
AtomTree::AtomTree( AtomTree const & src ) :
	utility::pointer::ReferenceCount(),
	this_weak_ptr_( /* 0 */ ),
	root_( /* 0 */ ), /// without this initialization, the destruction of this
	/// uninitialized pointer might have disasterous consequences
	atom_pointer_(), // default_setting_ = null pointer
	internal_coords_need_updating_( false ),
	xyz_coords_need_updating_( false ),
	topological_match_to_( /* 0 */ ),
	external_coordinate_residues_changed_( ResidueCoordinateChangeListOP( new ResidueCoordinateChangeList ) )
{
	*this = src;
}


void AtomTree::set_weak_pointer_to_self( AtomTreeCAP self_pointer )
{
	debug_assert( utility::pointer::equal(self_pointer, this) );
	this_weak_ptr_ = self_pointer;
}


/////////////////////////////////////////////////////////////////////////////
void
AtomTree::find_root_from_atom_pointer()
{
	root_.reset();
	for ( Size i=1; i<= atom_pointer_.size(); ++i ) {
		for ( Size j=1; j<= atom_pointer_[i].size(); ++j ) {
			debug_assert( atom_pointer_[i][j] && atom_pointer_[i][j]->id() == AtomID( j,i ) );
			if ( atom_pointer_[i][j]->parent() == nullptr ) {
				debug_assert( !root_ );
				root_ = atom_pointer_[i][j];
			}
		}
	}
}


/////////////////////////////////////////////////////////////////////////////
///
/// @details fill the AtomTree with a new tree of atoms by recording their pointers in
/// the map. Sync internal and xyz coords.
void
AtomTree::replace_tree(
	AtomPointer2D const & new_atom_pointer,
	bool const from_xyz // = true
)
{
	clear();

	atom_pointer_ = new_atom_pointer;
	find_root_from_atom_pointer();

	external_coordinate_residues_changed_->total_residue( atom_pointer_.size() );

	if ( from_xyz ) {
		internal_coords_need_updating_ = true;
		xyz_coords_need_updating_ = false;
		update_internal_coords();
	} else {
		xyz_coords_need_updating_ = true;
		internal_coords_need_updating_ = false;

		update_xyz_coords();
	}


	// anything that depends on the tree topology needs to be updated
	set_new_topology();
}


/// @details  This is a helper function to find a linear transform that when applied to the downstream stub
/// has the effect that RT( instub, transformed-downstream-stub) == target_rt
void
find_stub_transform(
	Stub const & stub1, // upstream stub
	Stub const & stub2, // downstream stub
	RT const & rt, // the target RT
	Stub::Matrix & A,
	Vector & b
)
{
	Stub::Matrix const & M1( stub1.M ), M2( stub2.M ), R( rt.get_rotation() );
	Vector const & v1( stub1.v ), v2( stub2.v ), t( rt.get_translation() );

	// look for a transformation of the form x |----> A*x + b
	//
	// this will change stub2 to stub2' with M2' = A * M2, v2' = A*v2 + b
	//
	// if we let (R,t) be the target RT, then we want
	//
	//  R = M1^T * M2' = M1^T * A * M2  ==> A = M1 * R * M2^T
	//
	//  t = M1^T * ( v2' - v1 ) ==> v2' = M1 * t + v1, which with b = v2' - A*v2 gives b = M1 * t + v1 - A * v2


	A = M1 * R * M2.transposed();
	b = M1 * t + v1 - A * v2;
}


/////////////////////////////////////////////////////////////////////////////
/// assumes one incoming and at most one outgoing
/// Need fancier fxn for general case
///
void
AtomTree::delete_seqpos( Size const seqpos )
{
	Size const old_size( size() ), new_size( old_size - 1 );

	// find the anchor, root, and perhaps child atoms
	Size const natoms( atom_pointer_[seqpos].size() );
	AtomOP anchor(nullptr), root(nullptr), child(nullptr);
	for ( Size i=1; i<= natoms; ++i ) {
		AtomOP atom( atom_pointer_[seqpos][i] );
		if ( !atom ) continue;
		if ( Size(atom->parent()->id().rsd()) != seqpos ) {
			debug_assert( !anchor );
			root = atom;
			anchor = atom->parent();
			// could break but debug 1 incoming connxn by continuing
		}
		for ( Atom::Atoms_ConstIterator iter=atom->atoms_begin(), iter_end = atom->atoms_end(); iter!= iter_end; ++iter ) {
			if ( Size((*iter)->id().rsd()) != seqpos ) {
				debug_assert( !child );
				child = *iter;
				// could break but debug at most 1 outgoing connxn by continuing
			}
		}
	}

	if ( !anchor ) {
		utility_exit_with_message( "AtomTree::delete_seqpos can't handle deleting the root residue");
	}

	// rewire connxns
	if ( child ) {
		anchor->replace_atom( root, child );
	} else {
		anchor->delete_atom( root );
	}

	//  // now delete atoms
	//  for ( Size i=1; i<= natoms; ++i ) {
	//   delete atom_pointer_[seqpos][i];
	//  }

	atom_pointer_[ seqpos ].resize(0);

	// now renumber
	utility::vector1< int > old2new( old_size );
	for ( Size i=1; i<= old_size; ++i ) {
		if      ( i <  seqpos ) old2new[i] = i;
		else if ( i == seqpos ) old2new[i] = 0;
		else                    old2new[i] = i-1;
	}

	update_sequence_numbering( new_size, old2new );

	// anything that depends on the tree topology needs to be updated
	set_new_topology();
}


/// @note  This can also handle the case of inserting a residue, proved that we've already renumbered atomtree leaving an empty slot at seqpos.
/// @note  If the new residue contains the root of the tree, incoming.atom1 should be BOGUS_ATOM_ID
/// @note  Note that for all BondID's atom1 should be the parent and atom2 should be the child

/// Couple of special cases:
/// -- appending a residue
/// -- inserting a residue
/// -- replacing the root residue
/// -- inserting a new root residue

void
AtomTree::replace_residue_subtree(
	id::BondID const & incoming,
	utility::vector1< id::BondID > const & outgoing,
	AtomPointer1D const & new_atoms // this will be the new slice of our atom_pointer
)
{
	// general rule:
	//
	// before we set or get a type of coordinate, have to
	// do an update
	//
	// in this case, we are "setting" the xyz's for the new subtree
	// this will put the internal coords out of date, so we want to
	// ensure that the xyz's are up to date at the start, otherwise
	// we can end up in the disastrous situation of having both
	// internal and xyz coords out of data -- and thus no "reference"
	// state from which to update!
	//
	update_xyz_coords();

	Size const seqpos( incoming.atom2.rsd() ); // should be Size but AtomID::rsd() returns int
	AtomPointer1D const & old_atoms( atom_pointer_[ seqpos ] );

	// confirm that bond id's go from parent to child, atom1 to atom2
	for ( Size i=1; i<= outgoing.size(); ++i ) debug_assert( outgoing[i].atom1.rsd() == seqpos );
	for ( Size i=1; i<= new_atoms.size(); ++i ) debug_assert( new_atoms[i]->id() == AtomID( i, seqpos ) );


	AtomOP anchor_atom(nullptr);
	AtomOP old_root_atom(nullptr);
	AtomOP new_root_atom( new_atoms[ incoming.atom2.atomno() ] );

	if ( incoming.atom1.valid() ) anchor_atom = atom_pointer( incoming.atom1 );

	// rewire the outgoing connections -- these are added to new atoms in the order that they come in the outgoing list
	// the children in these connections will all have parents in seqpos unless we are inserting into an empty slot,
	// ie unless old_atoms.empty()
	for ( Size i=1; i<= outgoing.size(); ++i ) {
		AtomOP child( atom_pointer( outgoing[i].atom2 ) );
		AtomOP old_parent( child->parent() );
		debug_assert( child->id().rsd() != seqpos );
		if ( !old_parent ) {
			// we're becoming the new root residue
			debug_assert( old_atoms.empty() && !incoming.atom1.valid() ); // implies anchor_atom == 0
			debug_assert( !old_root_atom );
			old_root_atom = child;
		} else if ( old_parent->id().rsd() != seqpos ) {
			// we're inserting into a bond
			debug_assert( old_atoms.empty() && old_parent->id() == incoming.atom1 );
			debug_assert( !old_root_atom );
			old_root_atom = child;
		} else {
			// only necessary for debugging purposes
			old_parent->delete_atom( child );
		}
		AtomOP new_parent( new_atoms[ outgoing[i].atom1.atomno() ] );
		new_parent->insert_atom( child );
	}

	// potentially have to look for old_root_atom
	if ( old_root_atom ) {
		debug_assert( old_atoms.empty() );
	} else {
		for ( Size i=1; i<= old_atoms.size(); ++i ) {
			AtomOP old_atom( old_atoms[i] );
			debug_assert( old_atom ); // atom_pointer_ is ragged, always keep dimension equal to actual number of atoms
			if ( ! old_atom->parent() ) {
				// this was the root of the atomtree
				debug_assert( !incoming.atom1.valid() );
				debug_assert( !old_root_atom );
				old_root_atom = old_atom;
			} else if ( old_atom->parent()->id().rsd() != seqpos ) {
				// this is the root of the old tree
				debug_assert( incoming.atom1 == old_atom->parent()->id() );
				debug_assert( !old_root_atom );
				old_root_atom = old_atom;
			}
			// this is just debugging to confirm that outgoing vector actually contains all the outgoing connections
			for ( Size i=0; i< old_atom->n_children(); ++i ) debug_assert( old_atom->child(i)->id().rsd() == seqpos );
		}
	}

	// rewire the incoming connection
	if ( anchor_atom ) {
		if ( old_root_atom ) {
			atom_pointer( incoming.atom1 )->replace_atom( old_root_atom, new_root_atom );
		} else {
			debug_assert( outgoing.empty() && old_atoms.empty() );
			atom_pointer( incoming.atom1 )->insert_atom( new_root_atom );
		}

	} else {
		debug_assert( root_ == old_root_atom );
		root_ = new_root_atom;
		new_root_atom->parent(tree::AtomAP());
	}


	// now nuke the old atoms
	atom_pointer_[ seqpos ].clear();
	atom_pointer_[ seqpos ] = new_atoms;


	// we've added the new atoms assuming that their xyz coords are
	// valid, but their internal coords (especially at the junctions
	// of the old and new) are likely to be messed up. This is also
	// true, eg, for any younger siblings of subtree_root
	//
	internal_coords_need_updating_ = true;

	// anything that depends on the tree topology needs to be updated
	set_new_topology();

}


/////////////////////////////////////////////////////////////////////////////
/// @brief retrieve a specific DOF by its ID.
Real
AtomTree::dof( DOF_ID const & id ) const
{
	update_internal_coords();
	return atom_pointer_[ id.atom_id() ]->dof( id.type() );
}

/////////////////////////////////////////////////////////////////////////////
/// @brief retrieve the xyz position of an atom in the tree by its AtomID
PointPosition const &
AtomTree::xyz( AtomID const & id ) const
{
	update_xyz_coords();
	// if ( !has(id) ) {
	//   std::cerr << "AtomTree::atom_pointer_ has not the atom " << id << std::endl;
	//  }
	//  runtime_assert( has( id ) );
	return atom_pointer_[ id ]->position();
}

/////////////////////////////////////////////////////////////////////////////
/// @brief retreive a kinematic Atom in the tree by its AtomID
tree::Atom const &
AtomTree::atom( AtomID const & id ) const
{
	// we don't know what kind of access the caller may perform:
	update_internal_coords();
	update_xyz_coords();

	return *(atom_pointer_[ id ]);
}

/////////////////////////////////////////////////////////////////////////////
/// @brief retreive a kinematic Atom in the tree by its AtomID -- no update!
tree::Atom const &
AtomTree::atom_dont_do_update( AtomID const & id ) const
{
	return *(atom_pointer_[ id ]);
}

/////////////////////////////////////////////////////////////////////////////
/// @brief retrieve a Jump in the tree by that JumpAtom's AtomID.
/// @note will abort if a BondedAtom's AtomID is passed in.
Jump const &
AtomTree::jump( AtomID const & id ) const
{
	update_internal_coords();
	return atom_pointer_[ id ]->jump();
}


/////////////////////////////////////////////////////////////////////////////
//get the DOF_ID of a torsion angle given those four atoms which define this torsion
/// @details an "offset" value is also calculated that torsion(id1,id2,id3,id4) = dof( dof_id ) + offset.
/// A BOGUS_DOF_ID will be returned if no proper DOF can be found for these four atoms.
/// offset is mainly for an atom with a previous sibling as the torsion(PHI) attached
/// to it is calculated as improper angle with respect to its sibling.
id::DOF_ID
AtomTree::torsion_angle_dof_id(
	AtomID const & atom1_in_id,
	AtomID const & atom2_in_id,
	AtomID const & atom3_in_id,
	AtomID const & atom4_in_id,
	Real & offset,
	bool const quiet
) const
{
	using numeric::conversions::degrees;
	using numeric::constants::d::pi_2;
	using numeric::constants::d::pi;
	using numeric::dihedral;
	using numeric::dihedral_radians;

	//TR << "Getting a torsion_angle_dof_id for " << atom1_in_id << "-" << atom2_in_id << "-" << atom3_in_id << "-" << atom4_in_id << std::endl;
	//bool const debug( false );

	// We use the internal DoFs if necessary to calculate the offset.
	update_internal_coords();

	//if ( debug ) update_xyz_coords();

	debug_assert( atom_pointer( atom1_in_id ) && atom_pointer( atom2_in_id ) &&
		atom_pointer( atom3_in_id ) && atom_pointer( atom4_in_id ) );

	// TODO: STUART -- (low priority) I'd like to be able to cache the results of this calculation
	// to allow faster access.

	Atom const * atom1_in( atom_pointer_( atom1_in_id ).get() ); // not atom_pointer_( ).get()
	Atom const * atom2_in( atom_pointer_( atom2_in_id ).get() );
	Atom const * atom3_in( atom_pointer_( atom3_in_id ).get() );
	Atom const * atom4_in( atom_pointer_( atom4_in_id ).get() );

	Atom const * atom1( atom1_in );
	Atom const * atom2( atom2_in );
	Atom const * atom3( atom3_in );
	Atom const * atom4( atom4_in );

	// Reorder the atoms if necessary.
	// We want it to be the case that atom4 has input_stub_atom1 == atom3 and input_stub_atom2 == atom2.
	// AMW TODO: keep_dof_fixed is super slow?
	if ( !atom4->is_jump() &&
			!atom4->keep_dof_fixed( id::PHI ) && // not stub atom3 of a jump
			atom4->raw_input_stub_atom1() == atom3 &&
			atom4->raw_input_stub_atom2() == atom2 ) {
		// pass -- this is what we want, not quite as perfect as 1st case though
	} else if ( !atom1->is_jump() &&
			!atom1->keep_dof_fixed( id::PHI ) && // not stub atom3 of a jump
			atom1->raw_input_stub_atom1() == atom2 &&
			atom1->raw_input_stub_atom2() == atom3 ) {
		// reverse the order of the atoms
		atom1 = atom4_in;
		atom2 = atom3_in;
		atom3 = atom2_in;
		atom4 = atom1_in;
	} else { /* else if (atom1_in_id.rsd() == atom2_in_id.rsd() &&
		atom1_in_id.rsd() != atom3_in_id.rsd() &&
		atom3_in_id.rsd() == atom4_in_id.rsd() ) {
		// BRANCH case: torsion includes 2 atoms from each of 2 residues
		TR << "Noncanonical connection torsion does not have a corresponding DOF_ID" << std::endl;
		return id::BOGUS_DOF_ID;
		}*/

		// AMW: We should note that this is FINE. This doesn't mean that we cannot have a fine
		// time with the data we need. If we return id:BOGUS_DOF_ID, then conformation will
		// catch that and just calculate the torsion manually (in atom_tree_torsion) or return
		// said BOGUS_DOF_ID (in dof_id_from_torsion_id)
		// and that, on this level, is what's desired.
		if ( !quiet ) {
			TR.Error << "No proper DoF can be found for these four atoms: ";
			TR.Error << atom1_in_id.rsd() << "-" << atom1_in_id.atomno() << ", ";
			TR.Error << atom2_in_id.rsd() << "-" << atom2_in_id.atomno() << ", ";
			TR.Error << atom3_in_id.rsd() << "-" << atom3_in_id.atomno() << ", ";
			TR.Error << atom4_in_id.rsd() << "-" << atom4_in_id.atomno() << "!" << std::endl;

			/*
			TR.Error << "What condition failed?" << std::endl;
			TR.Error << "!atom4->is_jump() = " << !atom4->is_jump() << std::endl;
			TR.Error << "!atom4->keep_dof_fixed( id::PHI ) = " << !atom4->keep_dof_fixed( id::PHI ) << std::endl;
			TR.Error << "atom4->raw_input_stub_atom1() == ( atom3 = " << atom3_in_id.rsd() << "-" << atom3_in_id.atomno() << " ) = " << ( atom4->raw_input_stub_atom1() == atom3 ) << std::endl;
			TR.Error << "atom4->raw_input_stub_atom2() == ( atom2 = " << atom2_in_id.rsd() << "-" << atom2_in_id.atomno() << " ) = " << ( atom4->raw_input_stub_atom2() == atom2 ) << std::endl;
			TR.Error << " or... " << std::endl;

			TR.Error << "!atom1->is_jump() = " << !atom1->is_jump() << std::endl;
			TR.Error << "!atom1->keep_dof_fixed( id::PHI ) = " << !atom1->keep_dof_fixed( id::PHI ) << std::endl;
			TR.Error << "atom1->raw_input_stub_atom1() == ( atom2 = " << atom2_in_id.rsd() << "-" << atom2_in_id.atomno() << " ) = " << ( atom1->raw_input_stub_atom1() == atom2 ) << std::endl;
			TR.Error << "atom1->raw_input_stub_atom2() == ( atom3 = " << atom3_in_id.rsd() << "-" << atom3_in_id.atomno() << " ) = " << ( atom1->raw_input_stub_atom2() == atom3 ) << std::endl;
			*/
		}
		return id::BOGUS_DOF_ID;
	}

	offset = 0.0; // initialize

	if ( atom4->raw_input_stub_atom3() == atom1 ) {
		// perfect match
		return DOF_ID( atom4->id(), id::PHI );
	}

	debug_assert( !atom4->is_jump() &&
		atom4->raw_input_stub_atom0() == atom3 &&
		atom4->raw_input_stub_atom1() == atom3 &&
		atom4->raw_input_stub_atom2() == atom2 &&
		atom4->raw_parent() == atom3 );

	// special case if atoms 1 and 4 are siblings -- not really well defined!
	if ( atom1->raw_parent() == atom3 ) {
		Size const atom1_index( atom3->raw_child_index( atom1 ) ),
			atom4_index( atom3->raw_child_index( atom4 ) );
		if ( atom1_index < atom4_index ) {
			Real const current_value( atom3->dihedral_between_bonded_children( *atom1, *atom4 ) );
			offset = current_value - atom4->dof(id::PHI); // since torsion(id1...id4) = dof + offset;

			ASSERT_ONLY( Real const actual_current_value
				( dihedral_radians( atom4->xyz(), atom3->xyz(), atom2->xyz(), atom1->xyz() ) ); );
			debug_assert( std::abs( basic::subtract_radian_angles( actual_current_value, current_value ) ) < 1e-3 );
			debug_assert( std::abs( basic::subtract_radian_angles( current_value, atom4->dof( id::PHI ) + offset ) ) < 1e-3 );

			return DOF_ID( atom4->id(), id::PHI );
		} else {
			Real const current_value( atom3->dihedral_between_bonded_children( *atom1, *atom4 ) );
			offset = current_value - atom1->dof( id::PHI ); // since torsion(id1...id4) = dof + offset;
			ASSERT_ONLY( Real const actual_current_value
				( dihedral_radians( atom4->xyz(), atom3->xyz(), atom2->xyz(), atom1->xyz() ) ););
			debug_assert( std::abs( basic::subtract_radian_angles( actual_current_value, current_value ) ) < 1e-3 );
			debug_assert( std::abs( basic::subtract_radian_angles( current_value, atom1->dof( id::PHI ) + offset ) ) < 1e-3 );
			return DOF_ID( atom1->id(), id::PHI );
		}
	}
	// atom4 is not the first sibling of atom3, get offset for that.
	if ( atom4 != atom3->raw_get_nonjump_atom(0) ) {
		Atom const * new_atom4( atom3->raw_get_nonjump_atom(0) );
		offset += atom3->dihedral_between_bonded_children( *new_atom4, *atom4 );

		/* if ( debug ) { // debugging
		ASSERT_ONLY( Real const actual_dihedral
		( dihedral_radians( atom4->xyz(), atom3->xyz(), atom2->xyz(),
		new_atom4->xyz() ) );)
		debug_assert( std::abs( basic::subtract_radian_angles
		( actual_dihedral, offset ) ) < 1e-3 );
		} */

		atom4 = new_atom4;
	}

	DOF_ID dof_id( atom4->id(), id::PHI );

	Atom const * dof_atom1( atom4->raw_input_stub_atom3() );

	if ( dof_atom1 == atom1 ) {
		return dof_id;
	} else if ( atom1->raw_parent() == atom2 && !atom1->is_jump() &&
			atom3->raw_parent() == atom2 && !atom3->is_jump() ) {

		// handle offset between atom1 and dof_atom1
		// the only case we can do is if atom1->parent() == atom2
		//
		// this will happen, eg for chi1 if we are folding c2n:
		//
		// atom1 =  N    while    dof_atom1 = parent(atom2) =  C
		// atom2 = CA
		// atom3 = CB
		// atom4 = CG


		// a little tricky: we don't want to access the positions,
		// since we can't be sure that they are up to date in a routine
		// that would be called inside dof setting and getting routines
		//
		//
		// need to use some spherical geometry
		//
		// we've got three unit vectors:
		//
		// u1 -- from atom2 to atom1
		// u2 -- from atom2 to dof_atom1
		// u3 -- from atom2 to atom3
		//
		// the angle between u2 and u1 = theta3 = pi - atom1(THETA)
		// the angle between u2 and u3 = theta1 = pi - atom3(THETA)
		// the torsion between u1 and u3 along u2 = phi2 = atom2->bonded_dih(1,3)
		Real
			theta1( pi - atom3->dof( id::THETA ) ),
			theta3( pi - atom1->dof( id::THETA ) ),
			phi2( atom2->dihedral_between_bonded_children( *atom1, *atom3 ) ),
			sign_factor( 1.0 );
		phi2 = basic::periodic_range( phi2, pi_2 );
		if ( phi2 < 0 ) {
			sign_factor = -1;
			phi2 *= -1;
		}
		//If bond angle is varied, these angles can sometimes go out of range,
		// during exploration by minimizer:
		theta3 = basic::periodic_range( theta3, pi_2 );
		if ( theta3 < 0 ) {
			sign_factor *= -1;
			theta3 *= -1;
		}
		theta1 = basic::periodic_range( theta1, pi_2 );
		if ( theta1 < 0 ) {
			sign_factor *= -1;
			theta1 *= -1;
		}
		if ( !( 0 <= theta1 && theta1 <= pi &&
				0 <= theta3 && theta3 <= pi &&
				0 <=   phi2 &&   phi2 <= pi ) ) {
			TR.Error << "dof_atom1 " << dof_atom1->id() << std::endl;
			TR.Error << "atom1 " << atom1->id() << std::endl;
			TR.Error << "atom2 " << atom2->id() << std::endl;
			TR.Error << "atom3 " << atom3->id() << std::endl;
			TR.Error << "atom4 " << atom4->id() << std::endl;
			TR.Error << "THETA1 " << theta1 << std::endl;
			TR.Error << "THETA3 " << theta3 << std::endl;
			TR.Error << "PHI2 " << phi2 << std::endl;
			TR.Error << "Atom1 position " << ObjexxFCL::format::F(8,3,atom1->xyz().x()) << ObjexxFCL::format::F(8,3,atom1->xyz().y()) << ObjexxFCL::format::F(8,3,atom1->xyz().z()) << std::endl;
			TR.Error << "Atom2 position " << ObjexxFCL::format::F(8,3,atom2->xyz().x()) << ObjexxFCL::format::F(8,3,atom2->xyz().y()) << ObjexxFCL::format::F(8,3,atom2->xyz().z()) << std::endl;
			TR.Error << "Atom3 position " << ObjexxFCL::format::F(8,3,atom3->xyz().x()) << ObjexxFCL::format::F(8,3,atom3->xyz().y()) << ObjexxFCL::format::F(8,3,atom3->xyz().z()) << std::endl;
			TR.Error << "Atom4 position " << ObjexxFCL::format::F(8,3,atom4->xyz().x()) << ObjexxFCL::format::F(8,3,atom4->xyz().y()) << ObjexxFCL::format::F(8,3,atom4->xyz().z()) << std::endl;

			utility_exit_with_message( "AtomTree::torsion_angle_dof_id: angle range error" );

		}
		// want to solve for phi3 -- dihedral between u2 and u1 along axis defined by u3
		//
		// spherical law of cosines says:
		//
		// cos(theta2) = cos(theta1) * cos(theta3) +
		//               sin(theta1) * sin(theta3) * cos(phi2)
		//
		// so we can solve for cos(theta2); then use formula again to solve for
		// cos(phi3).
		Real const cos_theta2 = std::cos( theta1 ) * std::cos( theta3 ) +
			std::sin( theta1 ) * std::sin( theta3 ) * std::cos( phi2 );
		Real const theta2 = std::acos( cos_theta2 );
		debug_assert( 0 <= theta2 && theta2 <= pi );

		// use formula again interchanging 2 and 3
		Real cos_phi3 = ( cos( theta3 ) - cos( theta1 ) * cos_theta2 ) /
			( sin( theta1 ) * sin( theta2 ) );

		// PHIL! really should handle degenerate cases more carefully
		// PHIL! also could reuse some of the trig calculations

		// dof-torsion: from dof_atom1 to atom4  ~  atom4 - dof_atom1
		// desired:     from atom1 to atom4  ~  atom4 - atom1
		//
		// atom4 - atom1 = atom4 - dof_atom1 + ( dof_atom1 - atom1 )
		//
		// so offset == dof_atom1 - atom1 ~ dihedral from atom1 to dof_atom1

		Real const dihedral_from_atom1_to_dof_atom1( sign_factor * std::acos( cos_phi3 ) );

		offset += dihedral_from_atom1_to_dof_atom1;

		/*if ( debug ) { // debugging
		using basic::periodic_range;
		using basic::subtract_radian_angles;
		Real const actual_dihedral_from_atom1_to_dof_atom1
		( dihedral_radians( atom1->xyz(), atom2->xyz(), atom3->xyz(),
		dof_atom1->xyz()) );

		Real const actual_dihedral_of_interest
		( dihedral_radians( atom1_in->xyz(), atom2_in->xyz(),
		atom3_in->xyz(), atom4_in->xyz()) );

		ASSERT_ONLY(Real const dev1
		( std::abs( subtract_radian_angles
		( actual_dihedral_from_atom1_to_dof_atom1,
		dihedral_from_atom1_to_dof_atom1 ) ) );)

		ASSERT_ONLY(Real const dev2
		( std::abs( subtract_radian_angles( actual_dihedral_of_interest,
		dof( dof_id ) + offset ) ) );)

		ASSERT_ONLY(Real const dev3
		( std::abs( subtract_radian_angles
		( dof( dof_id ),
		dihedral_radians( atom4->xyz(),
		atom4->input_stub_atom1()->xyz(),
		atom4->input_stub_atom2()->xyz(),
		atom4->input_stub_atom3()->xyz()))));)

		debug_assert( dev1 < 1e-3 && dev2 < 1e-3 && dev3 < 1e-3 );

		if ( false ) {
		TR.Trace << "offset: " << offset << " sgn_fac: " << sign_factor <<
		' ' << dihedral_from_atom1_to_dof_atom1 << " =?= " <<
		actual_dihedral_from_atom1_to_dof_atom1 <<
		' ' << periodic_range( dof( dof_id ) + offset, pi_2 ) <<
		" =?= " << actual_dihedral_of_interest << std::endl;
		}
		}*/

		return dof_id;
	}

	// failure
	return id::BOGUS_DOF_ID;
}


/////////////////////////////////////////////////////////////////////////////
// set a specific DOF in the tree
void
AtomTree::set_dof(
	DOF_ID const & id,
	Real const setting
)
{
	update_internal_coords();
	atom_pointer_[ id.atom_id() ]->set_dof( id.type(), setting, dof_changeset_ );
	xyz_coords_need_updating_ = true;
}

/////////////////////////////////////////////////////////////////////////////

void
AtomTree::set_xyz(
	AtomID const & id,
	PointPosition const & xyz
)
{
	update_xyz_coords();
	atom_pointer_[ id ]->position( xyz );
	internal_coords_need_updating_ = true;

}

/////////////////////////////////////////////////////////////////////////////

void
AtomTree::batch_set_xyz(
	utility::vector1<AtomID> const & ids,
	utility::vector1<PointPosition> const & xyzs
)
{
	runtime_assert( ids.size() == xyzs.size() );
	update_xyz_coords();
	for ( core::Size i=1; i<=ids.size(); ++i ) {
		atom_pointer_[ ids[i] ]->position( xyzs[i] );
	}
	internal_coords_need_updating_ = true;
}


/////////////////////////////////////////////////////////////////////////////
/// @note AtomID must point to a JumpAtom, otherwise it will die.
void
AtomTree::set_jump(
	AtomID const & id,
	Jump const & jump
)
{
	update_internal_coords();
	atom_pointer_[ id ]->jump( jump, dof_changeset_ );
	xyz_coords_need_updating_ = true;
}


/////////////////////////////////////////////////////////////////////////////
/// @note DEPRECATED.  AtomID must point to a JumpAtom, otherwise it will die.
/// Hopefully, the new output-sensitive refold will make this method unnecessary
void
AtomTree::set_jump_now(
	AtomID const & id,
	Jump const & jump
)
{
	// set_jump_now no longer necessary.  Commented-out implementation below is buggy
	// because the Conformation does not get its list of changed residues updated; so while
	// the coordinates in the atom tree are updated, the coordinates in the Conformation are not.
	// Instead of fixing this bug to handle this special case "set now", rely on the functional
	// and more efficient output-sensitive "set_jump"
	set_jump( id, jump );

	// We use our parent atoms' positions to find our stub,
	// which we can't do if their positions aren't correct.
	//if( xyz_coords_need_updating_ ) {
	// set_jump(id, jump);
	// return;
	//}
	//update_internal_coords();
	//atom_pointer_[ id ]->jump( jump );
	//Stub stub( atom_pointer_[ id ]->get_input_stub() );
	//atom_pointer_[ id ]->update_xyz_coords( stub );
	// No change to xyz_coords_need_updating_ -- we've done our part,
	// but previous changes may have invalidated other coords.
	//xyz_coords_need_updating_ = true;
}


/////////////////////////////////////////////////////////////////////////////
/// @details It is possible that no DOF can be obtained given these four atoms or the torsion
/// does not match exactly the DOF(PHI) angle. In the former case, a BOGUS_DOF_ID is
/// returned as indicator and in the latter case, an offset value is deducted from
/// input setting to set the real DOF value properly.
id::DOF_ID
AtomTree::set_torsion_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3,
	AtomID const & atom4,
	Real const setting,
	bool const quiet
)
{
	//TR << "amw can we get a dof for " << atom1.atomno() << "-" << atom2.atomno() << "-" << atom3.atomno() << "-" << atom4.atomno() << "?" << std::endl;
	Real offset;
	DOF_ID const & id( torsion_angle_dof_id( atom1, atom2, atom3, atom4, offset, quiet ) );

	if ( !id.valid() ) {
		if ( !quiet ) TR.Warning << "DOF for this torsion angle could not be found in AtomTree." << std::endl;
		return id;
	}

	// note: offset is defined (see torsion_angle_dof_id) so that
	//
	// torsion(atom1,atom2,atom3,atom4) = dof(id) + offset

	set_dof( id, setting - offset );

	return id;
}

/////////////////////////////////////////////////////////////////////////////

id::DOF_ID
AtomTree::set_bond_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3,
	Real const setting
)
{

	Real offset(0.0);
	DOF_ID const dof_id( bond_angle_dof_id( atom1, atom2, atom3, offset ) );

	if ( dof_id.valid() ) {
		debug_assert( offset == 0.0 ); // not handling this case right now

		set_dof( dof_id, numeric::constants::d::pi - setting );
	}

	return dof_id;
}

/////////////////////////////////////////////////////////////////////////////

Real
AtomTree::bond_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3
) const
{

	Real offset(0.0);
	DOF_ID const dof_id( bond_angle_dof_id( atom1, atom2, atom3, offset ) );

	if ( dof_id.valid() ) {
		debug_assert( offset == 0.0 ); // not handling this case right now
		return numeric::constants::d::pi - dof( dof_id );
	} else {
		TR << "unable to find DOF_ID for bond_angle: " << atom1 << ' ' << atom2 << ' ' << atom3 << std::endl;
	}

	return 0.0;
}

/////////////////////////////////////////////////////////////////////////////

id::DOF_ID
AtomTree::bond_angle_dof_id(
	AtomID const & atom1_in_id,
	AtomID const & atom2_in_id,
	AtomID const & atom3_in_id,
	Real & offset
) const
{
	offset = 0.0;

	// CASES
	// I.   the bond angle is the THETA dof of atom3 (atom3's parent is atom2 and atom2's parent is atom1)
	// II.  the bond angle is the THETA dof of atom1 (atom1's parent is atom2 and atom2's parent is atom3)
	// III. the bond angle is set by a torsion offset (atom1's parent is atom2 and atom3's parent is atom2)

	debug_assert( atom_pointer( atom1_in_id ) && atom_pointer( atom2_in_id ) && atom_pointer( atom3_in_id ) );

	AtomCOP
		atom1_in( atom_pointer( atom1_in_id ) ),
		atom2_in( atom_pointer( atom2_in_id ) ),
		atom3_in( atom_pointer( atom3_in_id ) );

	AtomCOP atom1( atom1_in ), atom2( atom2_in ), atom3( atom3_in );

	// reorder the atoms if necessary
	// not necessary at all in the current logic
	DOF_ID dof_id( id::BOGUS_DOF_ID );

	if ( !atom3->is_jump() &&
			!atom3->keep_dof_fixed( THETA ) && // not stub atom2 of a jump
			atom3->input_stub_atom1() == atom2 &&
			atom3->input_stub_atom2() == atom1 ) {
		// case I
		dof_id = DOF_ID( atom3->id(), THETA );

	} else if ( !atom1->is_jump() &&
			!atom1->keep_dof_fixed( THETA ) && // not stub atom2 of a jump
			atom1->input_stub_atom1() == atom2 &&
			atom1->input_stub_atom2() == atom3 ) {
		// case II, perfect case in reverse, want to reverse the order of the atoms
		atom1 = atom3_in;
		atom2 = atom2_in;
		atom3 = atom1_in;
		dof_id = DOF_ID( atom3->id(), THETA );

	} else {
		// not handling other cases right now
	}

	return dof_id;
}

/////////////////////////////////////////////////////////////////////////////

id::DOF_ID
AtomTree::set_bond_length(
	AtomID const & atom1,
	AtomID const & atom2,
	Real const setting
)
{

	DOF_ID const dof_id( bond_length_dof_id( atom1, atom2 ) );

	if ( dof_id.valid() ) {
		set_dof( dof_id, setting );
	}

	return dof_id;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Real
AtomTree::bond_length(
	AtomID const & atom1,
	AtomID const & atom2
) const
{

	DOF_ID const dof_id( bond_length_dof_id( atom1, atom2 ) );

	if ( dof_id.valid() ) {
		return dof( dof_id );
	} else {
		TR << "unable to find DOF_ID for bond_length: " << atom1 << ' ' << atom2 << std::endl;
	}

	return 0.0;
}

/////////////////////////////////////////////////////////////////////////////

id::DOF_ID
AtomTree::bond_length_dof_id(
	AtomID const & atom1_id,
	AtomID const & atom2_id
) const
{

	debug_assert( atom_pointer( atom1_id ) && atom_pointer( atom2_id ) );

	AtomCOP atom1( atom_pointer( atom1_id ) ), atom2( atom_pointer( atom2_id ) );

	if ( !atom2->is_jump() &&
			atom2->input_stub_atom1() == atom1 ) {
		// case I
		return DOF_ID( atom2->id(), D );

	} else if ( !atom1->is_jump() &&
			atom1->input_stub_atom1() == atom2 ) {
		// case II
		return DOF_ID( atom1->id(), D );

	} else {
		// not handling other cases right now
	}

	return id::BOGUS_DOF_ID;
}

/////////////////////////////////////////////////////////////////////////////
///
/// @note this is NOT calculated straight from atom positions.Instead, it is
/// calculated from internal DOFs. If no DOF can be found from these four atoms,
/// 0.0 will be returned and a warning message is printed.
Real
AtomTree::torsion_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3,
	AtomID const & atom4
) const
{
	update_internal_coords();

	// find the atomtree degree of freedom that corresponds to this torsion
	// angle
	Real offset;
	DOF_ID const & id
		( torsion_angle_dof_id( atom1, atom2, atom3, atom4, offset ) );

	if ( !id.valid() ) {
		// couldn't find this angle
		TR << "AtomTree::torsion_angle() cant find dof! " <<
			atom1 << ' ' << atom2 << ' ' << atom3 << ' ' << atom4 << std::endl;
		return 0.0;
	}

	// note: offset is defined (see torsion_angle_dof_id) so that
	//
	// torsion(atom1,atom2,atom3,atom4) = dof(id) + offset
	//
	//TR << "amw in torsion_angle" << std::endl;

	return atom_pointer_[ id.atom_id() ]->dof( id.type() ) + offset;
}

/////////////////////////////////////////////////////////////////////////////


/// @details This is done by releasing memories of all atoms including the root atom
///, and clearing atom_pointer map
void
AtomTree::clear()
{
	atom_pointer_.clear();
	root_.reset();
	external_coordinate_residues_changed_->total_residue(0);

	// anything that depends on the tree topology needs to be updated
	set_new_topology();

}

/////////////////////////////////////////////////////////////////////////////
/// @note this may trigger a coordinate update of src (access to root atom)


void
AtomTree::copy_coords(
	AtomTree const & src
)
{
	if ( !root_ ) {
		if ( src.root() ) utility_exit_with_message("AtomTree::copy_coords: I'm empty but src is not!");
		return;
	}
	internal_coords_need_updating_ = src.internal_coords_need_updating_;
	xyz_coords_need_updating_ = src.xyz_coords_need_updating_;

	root_->copy_coords( *(src.root()) );
	dof_changeset_ = src.dof_changeset_;
	(*external_coordinate_residues_changed_) = (*src.external_coordinate_residues_changed_);

}

////////////////////////////////////////////////////////////////////////////////////
/// @details domain map is residue-based and indicate which residues stay relatively rigid
/// to each other in one domain since last move and therefore their interaction
/// energy does not have to be recalculated. This function calls atom->update_domain_map
/// from the root atom of the tree.
void
AtomTree::update_domain_map(
	DomainMap & domain_map,
	AtomID_Mask const & dof_moved,
	AtomID_Mask const & xyz_moved
) const
{
	domain_map = -1;

	int current_color(1), biggest_color(1);
	root_->update_domain_map( current_color, biggest_color, domain_map,
		dof_moved, xyz_moved );
}


/////////////////////////////////////////////////////////////////////////////
// private
void
AtomTree::update_atom_ids_from_atom_pointer()
{
	for ( Size i=1, i_end = atom_pointer_.size(); i<= i_end; ++i ) {
		for ( Size j=1, j_end = atom_pointer_[i].size(); j<= j_end; ++j ) {
			debug_assert( atom_pointer_[i][j] );
			if ( atom_pointer_[i][j] ) atom_pointer_[i][j]->id( AtomID(j,i) );
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
void
AtomTree::update_sequence_numbering(
	Size const new_size,
	utility::vector1< int > const & old2new
)
{
	/// ResidueCoordinateChangeList is not setup to handle a remapping.  It must
	/// be empty, which means the Conformation its tracking moved data for must have
	/// already retrieved its moved data.
	debug_assert( external_coordinate_residues_changed_->empty() );
	external_coordinate_residues_changed_->total_residue( new_size );

	atom_pointer_.update_sequence_numbering( new_size, old2new );

	update_atom_ids_from_atom_pointer();
}


/////////////////////////////////////////////////////////////////////////////
/// @details makes a complete copy of another AtomTree src by cloning the root_atom and
/// all its offspring atoms. Also update atom_pointer map.
/// @note clear the content of tree first to release memory before doing the assignment.
AtomTree &
AtomTree::operator=( AtomTree const & src )
{
	if ( this == &src ) {
		return *this;
	}

	if ( utility::pointer::equal(topological_match_to_, &src) || utility::pointer::equal(src.topological_match_to_, this) ) {
		/// Why include the second condition: src.topological_match_to_ == this ?
		/// As an optimization for the common case when atom trees are being copied
		/// back and forth into each other as happens in repeated calls to MC::boltzman.
		/// Moreover, "topological match to" is a commutative property.
		copy_coords( src );
	} else {
		// no memory leaks -- notify observers that the topology is changing.
		clear();

		// just to get the dimensions right:
		atom_pointer_ = src.atom_pointer_;

		// copy tree (recursive)
		if ( src.root_ ) {
			root_ = src.root_->clone( tree::AtomAP() /*parent=0*/, atom_pointer_ );
		}

		internal_coords_need_updating_ = src.internal_coords_need_updating_;
		xyz_coords_need_updating_ = src.xyz_coords_need_updating_;
		dof_changeset_ = src.dof_changeset_;

		(*external_coordinate_residues_changed_) = (*src.external_coordinate_residues_changed_);

	}

	if ( !utility::pointer::equal(topological_match_to_, (&src)) ) {
		AtomTreeCOP topological_match_to( topological_match_to_.lock() );
		if ( topological_match_to ) {
			debug_assert( !this_weak_ptr_.expired() );
			topological_match_to->detach_topological_observer( this_weak_ptr_ );
		}
		/// topological observation only allowed if both this, and src hold weak pointers to themselves.
		/// If either AtomTree had been declared on the stack, or if the code that instantiated either one
		/// never gave them their weak-pointer-to-selves, then the observer system is bypassed.
		if ( !this_weak_ptr_.expired() && !src.this_weak_ptr_.expired() ) {
			src.attach_topological_observer( this_weak_ptr_ );
		}
	}
	return *this;

}

/// @details This works by invoking this AtomTree's opeator = and then resetting
/// the topological_match_to_ pointer by calling set_new_topology();
void
AtomTree::detached_copy( AtomTree const & src ) {
	*this = src;
	set_new_topology();
}

/////////////////////////////////////////////////////////////////////////////
/// @details update_internal_coords would be called in situations where we want
/// to ensure that the internal degrees are valid, e.g., when we are about
/// to set a new internal DOF.
/// @note private method, const to allow lazy updating of
/// @note see usage in AtomTree.cc
/// @note these guys are const because of the lazy updating scheme we are using:
/// they need to be called from within const accessing functions.
void
AtomTree::update_internal_coords() const
{
	// this would be bad:

	debug_assert( ! ( xyz_coords_need_updating_ && internal_coords_need_updating_ ) );

	if ( internal_coords_need_updating_ ) {
		if ( !root_ ) utility_exit_with_message("Attempting to update an AtomTree with no root!");

		PROF_START( basic::ATOM_TREE_UPDATE_INTERNAL_COORDS );  // profiling
		root_->update_internal_coords( default_stub );
		PROF_STOP ( basic::ATOM_TREE_UPDATE_INTERNAL_COORDS );

		internal_coords_need_updating_ = false;
	}
}

id::AtomID
AtomTree::get_jump_atom_id( StubID const& stub_id1,
	StubID const& stub_id2,
	int& direction ) const {
	AtomCOP jump_atom(nullptr);
	for ( Size i=1; i<= 3; ++i ) {
		AtomCOP atom1( atom_pointer( stub_id1.atom( i ) ) );
		for ( Size j=1; j<= 3; ++j ) {
			AtomCOP atom2( atom_pointer( stub_id2.atom( j ) ) );
			if ( atom1->is_jump() && atom1->parent() == atom2 ) {
				debug_assert( !jump_atom );
				jump_atom = atom1;
				direction = -1;
			} else if ( atom2->is_jump() && atom2->parent() == atom1 ) {
				debug_assert( !jump_atom );
				jump_atom = atom2;
				direction = 1;
			}
		}
	}

	return jump_atom->id();
}

/// @details  Set the transform between two stubs
/// Returns the atomid of the jump atom which moved.
/// @note  Requires that there be a jump in the atomtree between a pair of atoms in the stubs
id::AtomID
AtomTree::set_stub_transform(
	StubID const & stub_id1,
	StubID const & stub_id2,
	RT const & target_rt
)
{
	int dir = 0;
	AtomOP jump_atom = atom_pointer( get_jump_atom_id( stub_id1, stub_id2, dir ) );
	if ( !jump_atom ) {
		utility_exit_with_message( "AtomTree::set_stub_transform: No jump between these atoms!" );
	}

	Stub const  instub( dir == 1 ? stub_from_id( stub_id1 ) : stub_from_id( stub_id2 ) );
	Stub const outstub( dir == 1 ? stub_from_id( stub_id2 ) : stub_from_id( stub_id1 ) );
	RT rt( target_rt );
	if ( dir == -1 ) rt.reverse();
	// now solve for a linear transform that will move outstub so that RT( instub, outstub ) == rt
	Stub::Matrix A;
	Vector b;
	find_stub_transform( instub, outstub, rt, A, b );
	jump_atom->transform_Ax_plus_b_recursive( A, b, *external_coordinate_residues_changed_ );
	internal_coords_need_updating_ = true; // this could be more efficient!

	// confirm that things worked
	debug_assert( RT( stub_from_id( stub_id1 ), stub_from_id( stub_id2 ) ).distance_squared( target_rt ) < 1e-3 );
	return jump_atom->id();
}

/// @brief  get the transform between two stubs
RT
AtomTree::get_stub_transform(
	StubID const & stub_id1,
	StubID const & stub_id2
) const
{
	return RT( stub_from_id( stub_id1 ), stub_from_id( stub_id2 ) );
}
/////////////////////////////////////////////////////////////////////////////
// private
///\brief update xyz coordinates from internal cooridnates
void
AtomTree::update_xyz_coords() const
{
	// this would be bad:
	debug_assert( ! ( xyz_coords_need_updating_ && internal_coords_need_updating_ ) );


	if ( xyz_coords_need_updating_ ) {
		if ( !root_ ) utility_exit_with_message("phil how did we get here-2?");

		PROF_START( basic::ATOM_TREE_UPDATE_XYZ_COORDS ); // profiling
		for ( Size ii = 1; ii <= dof_changeset_.size(); ++ii ) {
			if ( dof_changeset_[ ii ].reached_ ) continue;
			atom_pointer_[ dof_changeset_[ ii ].atomid_ ]->dfs( dof_changeset_, *external_coordinate_residues_changed_, ii );
		}

		for ( Size ii = 1; ii <= dof_changeset_.size(); ++ii ) {
			if ( dof_changeset_[ ii ].reached_ ) continue;
			//std::cout << "Refold from " << dof_changeset_[ ii ].atomid_.rsd() << std::endl; // << " " << dof_changeset_[ ii ].atomid_.atomno() << " " << dof_changeset_[ ii ].reached_ << std::endl;
			atom_pointer_[ dof_changeset_[ ii ].atomid_ ]->update_xyz_coords(); // it must find its own stub.
		}
		dof_changeset_.clear();
		PROF_STOP ( basic::ATOM_TREE_UPDATE_XYZ_COORDS );

		xyz_coords_need_updating_ = false;

		//std::cout << "REFOLD END" << std::endl;

	}
}

/// @brief The AtomTree provides to the Conformation object a list of residues
/// whose xyz coordinates have changed.  When the Conformation has finished reading off
/// residues that have changed from the AtomTree, and has copied the coordinates of
/// those residues into its conformation::Residue objects, it informs the AtomTree
/// to reset this list by a call to mark_changed_residues_registered
///
/// @details The list of which residues have had coordinate changes is unknown until
/// the DFS has completed.  The DFS must be triggered before iterators are given to the
/// Conformation object.
ResidueListIterator
AtomTree::residue_xyz_change_list_begin() const
{
	if ( xyz_coords_need_updating_ ) update_xyz_coords();
	return external_coordinate_residues_changed_->residues_moved_begin();
}

/// @details The list of which residues have had coordinate changes is unknown until
/// the DFS has completed.  The DFS must be triggered before iterators are given to the
/// Conformation object.
ResidueListIterator
AtomTree::residue_xyz_change_list_end() const
{
	if ( xyz_coords_need_updating_ ) update_xyz_coords();
	return external_coordinate_residues_changed_->residues_moved_end();
}


/// @brief The AtomTree provides a list of residues who's xyz coordinates have changed
/// to the Conformation object.  When the Conformation has finished reading off residues
/// that have changed from the AtomTree, and has copied the coordinates of those residues
/// into its conformation::Residue objects, it informs the AtomTree to reset this list
/// by a call to mark_changed_residues_registered
void
AtomTree::note_coordinate_change_registered() const
{
	external_coordinate_residues_changed_->clear();
}


/// @brief When an atom tree copies the topology of another atom tree, it must
/// register itself as a topological observer of that other tree.  When the other
/// tree changes its topology, then the tree being observed must notify its
/// observers that they are no longer a topological copy of this tree.  An atom
/// tree may only be the topological copy of a single other atom tree, though several
/// atom trees may be copies of a single atom tree.
void
AtomTree::attach_topological_observer( AtomTreeCAP observer_ap ) const
{
	AtomTreeCOP observer( observer_ap );
	debug_assert( !observer->this_weak_ptr_.expired() && !this_weak_ptr_.expired() );
	debug_assert( observer->topological_match_to_.expired() );
	bool resize_topo_observers_array_( false );
	for ( Size ii = 1; ii <= topological_observers_.size(); ++ii ) {
		if ( topological_observers_[ ii ].expired() ) {
			resize_topo_observers_array_ = true;
			/// Make sure the observer is not already observing me
			debug_assert( !utility::pointer::equal(topological_observers_[ ii ], observer) );
		}
	}

	if ( resize_topo_observers_array_ ) {
		Size n_valid( 0 );
		for ( Size ii = 1; ii <= topological_observers_.size(); ++ii ) {
			if ( !topological_observers_[ ii ].expired() ) ++n_valid;
		}
		utility::vector1< AtomTreeCAP > valid_observers;
		valid_observers.reserve( n_valid + 1 );
		for ( Size ii = 1; ii <= topological_observers_.size(); ++ii ) {
			if ( !topological_observers_[ ii ].expired() ) valid_observers.push_back( topological_observers_[ ii ] );
		}
		topological_observers_.swap( valid_observers );
	}

	topological_observers_.push_back( observer );
	observer->topological_match_to_ = this_weak_ptr_;
}

/// @details The AtomTree being observed calls this notify_topological_change on all
/// of the AtomeTrees that are observing it.  The this object "detaches" itself from
/// the observee in this function.  In debug mode, this function ensures that this
/// object was actually observing the observee -- if it wasn't, then an internal error
/// has occurred.  The observee and the observer have gotten out of sync.
void
AtomTree::notify_topological_change( AtomTreeCAP ASSERT_ONLY( observee ) ) const
{
	debug_assert( utility::pointer::equal(observee, topological_match_to_) );
	topological_match_to_.reset();
}

/// @details When an AtomTree which was observing this AtomTree changes its
/// topology or is deleted, it must invoke this method on this AtomTree.
/// When this happens, this AtomTree marks the observer's position in
/// its list of observers as null -- and it also resets the
/// "topological_match_to_" pointer in the oberserver_ap AtomTree.
void
AtomTree::detach_topological_observer( AtomTreeCAP observer_ap ) const
{
	AtomTreeCOP observer( observer_ap );
	debug_assert( utility::pointer::equal(observer->topological_match_to_, this) );
#ifndef NDEBUG
	bool found( false );
#endif
	for ( Size ii = 1; ii <= topological_observers_.size(); ++ii ) {
		if ( utility::pointer::equal(topological_observers_[ ii ], observer) ) {
#ifndef NDEBUG
			found = true;
#endif
			topological_observers_[ ii ].reset();
			observer->topological_match_to_.reset();
			break;
		}
	}
	debug_assert( found );
}

/// @details If this tree changes its topology, then its no longer a match in topology
/// to the tree it was copied from, nor are the trees that are a copy of it.  Detach
/// this as an observer of the AtomTree it is observing (if any) and inform all the
/// trees observing this tree that the topology has changed before clearing the
/// topological_observers_ array.
void
AtomTree::set_new_topology()
{
	AtomTreeCOP topological_match_to = topological_match_to_.lock();
	if ( topological_match_to ) {
		topological_match_to->detach_topological_observer( this_weak_ptr_ );
	}
	for ( Size ii = 1; ii <= topological_observers_.size(); ++ii ) {
		AtomTreeCOP topological_observers_ii = topological_observers_[ ii ].lock();
		if ( topological_observers_ii ) {
			topological_observers_ii->notify_topological_change( this_weak_ptr_ );
		}
	}
	topological_observers_.clear();
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// CARTESIAN-COORDINATE FRAGMENT INSERTION ROUTINES ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void
AtomTree::get_frag_atoms(
	StubID const & id,
	FragXYZ const & frag_xyz,
	AtomCOP & frag_atom,
	AtomCOP & nonfrag_atom // could be zero if atom1 is root of tree
) const
{
	bool const atom1_in_frag( frag_xyz.count( id.atom1 ) );
	bool const atom2_in_frag( frag_xyz.count( id.atom2 ) );

	AtomCOP atom1( atom_pointer( id.atom1 ) );

	if ( atom1->is_jump() ) {
		if ( !atom1_in_frag && atom1->parent() && frag_xyz.count( atom1->parent()->atom_id() ) ) {
			frag_atom    = atom1->parent();
			nonfrag_atom = atom1;
			return;
		} else if ( atom1_in_frag && ( !atom1->parent() || !frag_xyz.count( atom1->parent()->atom_id() ) ) ) {
			frag_atom    = atom1;
			nonfrag_atom = atom1->parent();
			return;
		}
	}

	if ( atom1_in_frag && !atom2_in_frag ) {
		frag_atom    = atom_pointer( id.atom1 );
		nonfrag_atom = atom_pointer( id.atom2 );
		return;

	} else if ( atom2_in_frag && !atom1_in_frag ) {
		frag_atom    = atom_pointer( id.atom2 );
		nonfrag_atom = atom_pointer( id.atom1 );
		return;

	}
	utility_exit_with_message( "AtomTree::get_frag_atoms failed" );

}


/////////////////////////////////////////////////////////////////////////////
/// @details id is a frag atom
///
/// look for two more nearby frag atoms to define a pseudo-stub for getting coords for parents or children of this atom
id::StubID
AtomTree::get_frag_pseudo_stub_id(
	AtomID const & id,
	FragXYZ const & frag_xyz,
	bool & fail
) const
{
	debug_assert( frag_xyz.count( id ) );

	utility::vector1< AtomID > ids;

	AtomCOP atom( atom_pointer( id ) );
	AtomCOP parent( atom->parent() ), child1( atom->get_nonjump_atom(0) ), child2( atom->get_nonjump_atom(1) );

	if ( !atom->is_jump() && parent && frag_xyz.count( parent->atom_id() ) ) {
		ids.push_back( parent->atom_id() );
	}
	if ( child1 && frag_xyz.count( child1->atom_id() ) ) {
		ids.push_back( child1->atom_id() );
	}
	if ( child2 && frag_xyz.count( child2->atom_id() ) ) {
		ids.push_back( child2->atom_id() );
	}
	if ( child1 && frag_xyz.count( child1->atom_id() ) ) {
		AtomCOP gchild1( child1->get_nonjump_atom(0) );
		if ( gchild1 && frag_xyz.count( gchild1->atom_id() ) ) {
			ids.push_back( gchild1->atom_id() );
		}
	}
	if ( child2 && frag_xyz.count( child2->atom_id() ) ) {
		AtomCOP gchild2( child2->get_nonjump_atom(0) );
		if ( gchild2 && frag_xyz.count( gchild2->atom_id() ) ) {
			ids.push_back( gchild2->atom_id() );
		}
	}
	if ( !atom->is_jump() && parent && frag_xyz.count( parent->atom_id() ) ) {
		AtomCOP gparent( parent->parent() );
		if ( !parent->is_jump() && gparent && frag_xyz.count( gparent->atom_id() ) ) {
			ids.push_back( gparent->atom_id() );
		}
		for ( Size i=0; i<parent->n_children(); ++i ) {
			AtomCOP sibling( parent->child(i) );
			if ( sibling != atom && !sibling->is_jump() && frag_xyz.count( sibling->atom_id() ) ) {
				ids.push_back( sibling->atom_id() );
			}
		}
	}
	if ( ids.size() < 2 ) {
		fail = true;
		return StubID( id, id, id );
	}
	return StubID( id, ids[1], ids[2] );
}


/////////////////////////////////////////////////////////////////////////////
/// @details private helper for fragment insertion routines
Stub
AtomTree::get_frag_local_stub(
	StubID const & stubid,
	FragXYZ const & frag_xyz,
	bool & fail
) const
{
	// first look for special case of stub of frag jump-child:
	AtomID const stub_atom1_id( stubid.atom1 ), stub_atom2_id( stubid.atom2 ), stub_atom3_id( stubid.atom3 );
	AtomCOP stub_atom1( atom_pointer( stub_atom1_id ) );

	if ( !frag_xyz.count( stub_atom1_id ) && stub_atom1->is_jump() &&
			stub_atom1->parent() && frag_xyz.count( stub_atom1->parent()->atom_id()) &&
			stub_atom1->stub_atom1_id() == stub_atom1_id &&
			stub_atom1->stub_atom2_id() == stub_atom2_id &&
			stub_atom1->stub_atom3_id() == stub_atom3_id ) {

		/// special case: handled more easily this way
		StubID const instub_id( stub_atom1->input_stub_atom1_id(),
			stub_atom1->input_stub_atom2_id(),
			stub_atom1->input_stub_atom3_id() );
		debug_assert( stub_atom1->input_stub_atom0_id() == stub_atom1->input_stub_atom1_id() );

		// current xyz transform:
		RT const current_rt( stub_from_id( instub_id ), stub_from_id( stubid ) );
		Stub const instub( get_frag_local_stub( instub_id, frag_xyz, fail ) ); // recursive call
		Stub outstub;
		TR.Trace << "get_frag_local_stub:: making jump: " <<
			stubid.atom1 << ' ' << stubid.atom2 << ' ' << stubid.atom3 << ' ' <<
			instub_id.atom1 << ' ' << instub_id.atom2 << ' ' << instub_id.atom3 << std::endl;

		current_rt.make_jump( instub, outstub );
		return outstub;
	}


	return Stub( get_frag_local_xyz( stub_atom1_id, frag_xyz, fail ),
		get_frag_local_xyz( stub_atom2_id, frag_xyz, fail ),
		get_frag_local_xyz( stub_atom3_id, frag_xyz, fail ) );
}

/////////////////////////////////////////////////////////////////////////////
/// @details private helper for fragment insertion routines
///
/// id is either in frag or a child or a gchild or a parent
Vector
AtomTree::get_frag_local_xyz(
	AtomID const & id,
	FragXYZ const & frag_xyz,
	bool & fail
) const
{
	// easiest case:
	if ( frag_xyz.count( id ) ) return frag_xyz.find( id )->second;

	AtomCOP atom( atom_pointer(id) );

	if ( ( atom->parent() && frag_xyz.count( atom->parent()->atom_id() ) ) ||
			( atom->parent() && atom->parent()->parent() && frag_xyz.count( atom->parent()->parent()->atom_id() ) ) ) {
		// child or grand child
		return get_frag_descendant_local_xyz( atom, frag_xyz, fail );
	}

	// now should be parent of frag
	AtomCOP child( nullptr );
	for ( Size i=0; i< atom->n_children(); ++i ) {
		if ( frag_xyz.count( atom->child(i)->atom_id() ) ) {
			debug_assert( child == nullptr );
			child = atom->child(i);
		}
	}
	debug_assert( child ); // this is a known hole: could be asked for a grandparent of a frag, not just a parent. fix this phil
	return get_frag_parent_local_xyz( child, frag_xyz, fail );


}

/////////////////////////////////////////////////////////////////////////////
/// @details private helper for fragment insertion routines
Vector
AtomTree::get_frag_descendant_local_xyz(
	AtomCOP atom,
	FragXYZ const & frag_xyz,
	bool & fail
) const
{
	AtomID const id( atom->atom_id() );
	debug_assert( !frag_xyz.count( id ) );
	bool const frag_child( atom->parent() && frag_xyz.count( atom->parent()->atom_id() ) );
	ASSERT_ONLY(bool const frag_gchild
		( atom->parent() && atom->parent()->parent() && frag_xyz.count( atom->parent()->parent()->atom_id() ) );)
		debug_assert( frag_child || frag_gchild );

	AtomID id1( atom->input_stub_atom1_id() );
	AtomID id2( atom->input_stub_atom2_id() );
	AtomID id3( atom->input_stub_atom3_id() );
	debug_assert( atom->input_stub_atom0_id() == id1 ); // cant handle this case yet

	if ( id == id1 || id == id2 || id == id3 ) {
		/// circular!!! potential for infinite loop!
		if ( frag_child ) {
			StubID tmp( get_frag_pseudo_stub_id( atom->parent()->atom_id(), frag_xyz, fail ) );
			id1 = tmp.atom1;
			id2 = tmp.atom2;
			id3 = tmp.atom3;
		} else {
			StubID tmp( get_frag_pseudo_stub_id( atom->parent()->parent()->atom_id(), frag_xyz, fail ) );
			id1 = tmp.atom1;
			id2 = tmp.atom2;
			id3 = tmp.atom3;
		}
		if ( fail ) return Vector(0.0);
	}

	Stub const current_stub( stub_from_id( StubID( id1, id2, id3) ) );
	Stub const   local_stub( get_frag_local_stub( StubID( id1, id2, id3 ), frag_xyz, fail ) ); // potentially recursive
	return local_stub.local2global( current_stub.global2local( xyz( id ) ) );
}

/////////////////////////////////////////////////////////////////////////////
/// @details private helper for fragment insertion routines
Vector
AtomTree::get_frag_parent_local_xyz(
	AtomCOP child,
	FragXYZ const & frag_xyz,
	bool & fail
) const
{
	AtomCOP parent( child->parent() );
	//debug_assert( !parent->is_jump() ); // dont think we have to handle this case...
	debug_assert( frag_xyz.count( child->atom_id() ) && ! frag_xyz.count( parent->atom_id() ) );

	// build a pseudo-stub for the parent
	StubID const pseudo_stubid( get_frag_pseudo_stub_id( child->atom_id(), frag_xyz, fail ) );
	if ( fail ) return Vector(0.0);
	Stub const current_stub( stub_from_id( pseudo_stubid ) );
	Stub const   local_stub( get_frag_local_stub( pseudo_stubid, frag_xyz, fail ) );
	debug_assert( !fail ); // since pseudo stub atoms are all in the fragment
	return local_stub.local2global( current_stub.global2local( xyz( parent->atom_id() ) ) );
}


/////////////////////////////////////////////////////////////////////////////
// a private function
//
// the incoming stub should in fact be the true incoming connection
// and all the outgoing connections should be accounted for.
//
// This is all checked in assert statements.
void
AtomTree::insert_single_fragment(
	StubID const & instub_id,
	FragRT const & outstub_transforms,
	FragXYZ const & frag_xyz,
	utility::vector1< AtomID > & moving_atoms
)
{
	// get the incoming stub:
	Stub const instub( stub_from_id( instub_id ) );

	AtomCOP instub_frag_atom(nullptr), instub_nonfrag_atom(nullptr);
	get_frag_atoms( instub_id, frag_xyz, instub_frag_atom, instub_nonfrag_atom );

	debug_assert( instub_frag_atom->parent() == instub_nonfrag_atom ); // sanity check

	utility::vector1< AtomCOP > outstub_nonfrag_atoms; // just for debugging

	 for ( auto const & outstub_transform : outstub_transforms ) {
		StubID const & outstub_id( outstub_transform.first );
		AtomCOP outstub_frag_atom(nullptr), outstub_nonfrag_atom(nullptr);
		get_frag_atoms( outstub_id, frag_xyz, outstub_frag_atom, outstub_nonfrag_atom );
		outstub_nonfrag_atoms.push_back( outstub_nonfrag_atom ); // for debugging
		debug_assert( outstub_nonfrag_atom->parent() == outstub_frag_atom );

		// now transform outstub_nonfrag_atom so that current transform becomes desired transform
		//
		Stub const & stub1( instub ); // to match names below
		Stub const   stub2( stub_from_id( outstub_id ) );
		RT const & rt( outstub_transform.second );

		// now things are arranged so that when we transform outstub_nonfrag_atom and children, stub1 stays fixed and
		// stub2 moves, and we want RT( stub1, new_stub2 ) == target_rt

		Stub::Matrix const & M1( stub1.M ), M2( stub2.M ), R( rt.get_rotation() );
		Vector const & v1( stub1.v ), v2( stub2.v ), t( rt.get_translation() );

		// look for a transformation of the form x |----> A*x + b
		//
		// this will change stub2 to stub2' with M2' = A * M2, v2' = A*v2 + b
		//
		// if we let (R,t) be the target RT, then we want
		//
		//  R = M1^T * M2' = M1^T * A * M2  ==> A = M1 * R * M2^T
		//
		//  t = M1^T * ( v2' - v1 ) ==> v2' = M1 * t + v1, which with b = v2' - A*v2 gives b = M1 * t + v1 - A * v2


		Stub::Matrix const A( M1 * R * M2.transposed() );
		Vector const b( M1 * t + v1 - A * v2 );

		atom_pointer_[ outstub_nonfrag_atom->atom_id() ]->
			transform_Ax_plus_b_recursive( A, b, *external_coordinate_residues_changed_ ); // get nonconst version
		internal_coords_need_updating_ = true;

		moving_atoms.push_back( outstub_nonfrag_atom->atom_id() );
	}


	for ( auto it=frag_xyz.begin(), ite= frag_xyz.end(); it != ite; ++it ) {
		AtomID const & id( it->first );
		AtomOP atom( atom_pointer( id ) );

		// update xyz using instub
		atom->xyz( instub.local2global( it->second ) );

		{ // sanity check/debug
			// if parent not in frag, assert this is instub_frag_atom
			debug_assert( ( atom->parent() && frag_xyz.count( atom->parent()->atom_id() ) ) || atom == instub_frag_atom );

			// if has a child not in frag, assert child is in outstub_nonfrag_atoms
			for ( Size i=0; i< atom->n_children(); ++i ) {
				debug_assert( frag_xyz.count( atom->child(i)->atom_id() ) ||
					( std::find( outstub_nonfrag_atoms.begin(), outstub_nonfrag_atoms.end(), atom->child(i) ) !=
					outstub_nonfrag_atoms.end()));
			}
		}
		moving_atoms.push_back( id );

	}
	internal_coords_need_updating_ = true;


	// now debug transforms
	for ( auto it= outstub_transforms.begin(), ite= outstub_transforms.end(); it != ite; ++it ) {
		debug_assert( RT( instub, stub_from_id( it->first )).distance_squared( it->second ) < 1e-3 );
	}


}


/////////////////////////////////////////////////////////////////////////////
void
AtomTree::set_jump_atom_stub_id(
	StubID const & id
)
{
	update_xyz_coords(); // since we will need to recalculate all the internal dofs when we're done

	AtomOP atom1( atom_pointer( id.atom1 ) );
	AtomOP atom2( atom_pointer( id.atom2 ) );
	AtomOP atom3( atom_pointer( id.atom3 ) );
	if ( !atom1->is_jump() || atom2->is_jump() || atom3->is_jump() ||
			atom2->parent() != atom1 || atom3->parent() != atom2 ) {
		utility_exit_with_message( "set_jump_atom_stub_id failed!" );
	}

	atom1->delete_atom( atom2 );
	atom1->insert_atom( atom2 ); // goes to head of line
	atom2->delete_atom( atom3 );
	atom2->insert_atom( atom3 ); // goes to head of line

	internal_coords_need_updating_ = true;

	debug_assert( atom1->stub_atom1() == atom1 && atom1->stub_atom2() == atom2 && atom1->stub_atom3() == atom3 );

}


/////////////////////////////////////////////////////////////////////////////
// completely new plan
void
AtomTree::insert_fragment(
	StubID const & instub_id,
	FragRT const & outstub_transforms,
	FragXYZ const & frag_xyz,
	utility::vector1< AtomID > & moving_atoms
)
{
	// first compile a list of incoming/outgoing connections

	// look for more incoming and/or outgoing connections:

	utility::vector1< AtomCOP > incoming_stub_frag_atoms, outgoing_stub_nonfrag_atoms;
	utility::vector1< Stub > incoming_local_stubs, outgoing_local_stubs;
	utility::vector1< StubID > incoming_stub_ids, outgoing_stub_ids;

	{
		AtomCOP frag_atom, nonfrag_atom;
		get_frag_atoms( instub_id, frag_xyz, frag_atom, nonfrag_atom );
		if ( frag_atom->parent() == nonfrag_atom ) {
			TR.Trace << "AtomTree::insert_fragment: instub is incoming" << std::endl;
			incoming_stub_frag_atoms.push_back( frag_atom );
			incoming_local_stubs.push_back( Stub() );
			incoming_stub_ids.push_back( instub_id );
		} else if ( nonfrag_atom && nonfrag_atom->parent() == frag_atom ) {
			TR.Trace << "AtomTree::insert_fragment: instub is outgoing" << std::endl;
			outgoing_stub_nonfrag_atoms.push_back( nonfrag_atom );
			outgoing_local_stubs.push_back( Stub() );
			outgoing_stub_ids.push_back( instub_id );
		}
	} // scope

	 for ( auto const & outstub_transform : outstub_transforms ) {
		StubID const & outstub_id( outstub_transform.first );
		AtomCOP frag_atom, nonfrag_atom;
		get_frag_atoms( outstub_id, frag_xyz, frag_atom, nonfrag_atom );
		if ( nonfrag_atom && nonfrag_atom->parent() == frag_atom ) {
			TR.Trace << "AtomTree::insert_fragment: outstub is outgoing" << std::endl;
			outgoing_stub_nonfrag_atoms.push_back( nonfrag_atom );
			outgoing_local_stubs.push_back( Stub( outstub_transform.second ) );
			outgoing_stub_ids.push_back( outstub_id );
		} else if ( frag_atom->parent() == nonfrag_atom ) {
			TR.Trace << "AtomTree::insert_fragment: outstub is incoming" << std::endl;
			incoming_stub_frag_atoms.push_back( frag_atom );
			incoming_local_stubs.push_back( Stub( outstub_transform.second ) );
			incoming_stub_ids.push_back( outstub_id );
		}
	}

	for ( auto it=frag_xyz.begin(), ite= frag_xyz.end(); it != ite; ++it ) {
		//AtomID const & id( it->first );
		AtomOP atom( atom_pointer( it->first ) );
		if ( ( atom->parent() == nullptr || !frag_xyz.count( atom->parent()->atom_id() ) ) &&
				( std::find( incoming_stub_frag_atoms.begin(), incoming_stub_frag_atoms.end(), atom ) ==
				incoming_stub_frag_atoms.end() ) ) {
			// found a new incoming connection!
			TR.Trace << "AtomTree::insert_fragment: found new incoming connection: " << atom->atom_id() << std::endl;
			// now want to generate a StubID for this guy as well as a local stub
			if ( atom->is_jump() ) {
				// incoming jump-atom connection -- or root of tree??
				// get local coords for all the stub atoms
				bool fail( false );
				StubID const stubid( atom->stub_atom1_id(), atom->stub_atom2_id(), atom->stub_atom3_id() );
				Stub const stub( get_frag_local_stub( stubid, frag_xyz, fail ) );
				if ( !fail ) {
					incoming_stub_frag_atoms.push_back( atom );
					incoming_local_stubs.push_back( stub );
					incoming_stub_ids.push_back( stubid );
				}
				debug_assert( !fail ); // could fail
			} else {
				// incoming bonded-atom connection
				bool fail( false );
				StubID const stubid( get_frag_pseudo_stub_id( atom->atom_id(), frag_xyz, fail ) );
				if ( !fail ) {
					incoming_stub_frag_atoms.push_back( atom );
					incoming_stub_ids.push_back( stubid );
					incoming_local_stubs.push_back( get_frag_local_stub( stubid, frag_xyz, fail ) );
					debug_assert( !fail ); // shouldnt fail for a frag_pseudo_stub -- all atoms are already in frag
				}
				debug_assert( !fail ); // could fail
			}
		}

		// look for outgoing connections:
		for ( Size i=0; i< atom->n_children(); ++i ) {
			AtomOP child( atom->child(i) );
			if ( !frag_xyz.count( child->atom_id() ) &&
					( std::find( outgoing_stub_nonfrag_atoms.begin(), outgoing_stub_nonfrag_atoms.end(), child ) ==
					outgoing_stub_nonfrag_atoms.end() ) ) {
				// new outgoing connection!
				TR.Trace << "AtomTree::insert_fragment: found new outgoing connection: " << atom->atom_id() << std::endl;
				bool fail( false );
				StubID const stubid( child->stub_atom1_id(), child->stub_atom2_id(), child->stub_atom3_id() );
				Stub const outstub( get_frag_local_stub( stubid, frag_xyz, fail ) );
				if ( !fail ) {
					outgoing_stub_nonfrag_atoms.push_back( child );
					outgoing_stub_ids.push_back( stubid );
					outgoing_local_stubs.push_back( outstub );
				}
				debug_assert( !fail );
			}
		}
	}


	/// Now associate outgoing stubs and frag atoms with the incoming stubs, call insert_single_fragment once for each
	/// incoming stub.
	///
	for ( Size i=1; i<= incoming_stub_frag_atoms.size(); ++i ) {
		AtomCOP instub_atom( incoming_stub_frag_atoms[i] );
		Stub const & local_instub( incoming_local_stubs[ i ] );
		FragXYZ new_frag_xyz;
		FragRT new_outstub_transforms;
		// get frag atoms that depend on this guy:
		 for ( auto const & it : frag_xyz ) {
			AtomID const & id( it.first );
			AtomCOP frag_atom( atom_pointer( id ) );
			if ( frag_atom->atom_is_on_path_from_root( instub_atom ) ) {
				new_frag_xyz[ id ] = local_instub.global2local( it.second );
			}
		}
		for ( Size j=1; j<= outgoing_stub_nonfrag_atoms.size(); ++j ) {
			AtomCOP outstub_atom( outgoing_stub_nonfrag_atoms[j] );
			if ( outstub_atom->atom_is_on_path_from_root( instub_atom ) ) {
				new_outstub_transforms[ outgoing_stub_ids[j] ] = RT( local_instub, outgoing_local_stubs[j] );
			}
		}

		TR.Trace << "AtomTree::insert_fragment: inserting single fragment nout= " << new_outstub_transforms.size() <<
			" natoms= " << new_frag_xyz.size() << std::endl;

		insert_single_fragment( incoming_stub_ids[i], new_outstub_transforms, new_frag_xyz, moving_atoms );
	}


}


/// @details  Useful for guaranteeing that a stub remains within a single residue
void
AtomTree::promote_sameresidue_nonjump_child( AtomID const & parent_atom_id )
{
	update_xyz_coords(); // since we are about to invalidate the internal coords

	AtomOP parent( atom_pointer( parent_atom_id ) );
	tree::AtomOP sameresidue_child( nullptr );
	for ( Size i=0; i< parent->n_nonjump_children(); ++i ) {
		tree::AtomOP child( atom_pointer( parent->get_nonjump_atom( i )->id() ) ); // want nonconst, use atom_pointer
		if ( child->id().rsd() == parent->id().rsd() ) {
			sameresidue_child = child;
			break;
		}
	}
	if ( sameresidue_child ) {
		debug_assert( !sameresidue_child->is_jump() );
		parent->delete_atom( sameresidue_child );
		parent->insert_atom( sameresidue_child );
	} else {
		TR.Warning << "promote_sameresidue_nonjump_child failed, parent has no non-jump, same-residue children!" <<
			std::endl;
	}

	internal_coords_need_updating_ = true;
	set_new_topology();
}

} // namespace kinematics
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::kinematics::AtomTree::save( Archive & arc ) const {
	//arc( CEREAL_NVP( this_weak_ptr_ ) ); // AtomTreeCAP
	// EXEMPT this_weak_ptr_ topological_match_to_ topological_observers_
	arc( CEREAL_NVP( root_ ) ); // AtomOP
	arc( CEREAL_NVP( atom_pointer_ ) ); // AtomPointer2D
	arc( CEREAL_NVP( internal_coords_need_updating_ ) ); // _Bool
	arc( CEREAL_NVP( xyz_coords_need_updating_ ) ); // _Bool
	//arc( CEREAL_NVP( topological_match_to_ ) ); // AtomTreeCAP
	//arc( CEREAL_NVP( topological_observers_ ) ); // utility::vector1<AtomTreeCAP>
	arc( CEREAL_NVP( dof_changeset_ ) ); // AtomDOFChangeSet
	arc( CEREAL_NVP( external_coordinate_residues_changed_ ) ); // ResidueCoordinateChangeListOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::kinematics::AtomTree::load( Archive & arc ) {
	//arc( this_weak_ptr_ ); // AtomTreeCAP
	// EXEMPT this_weak_ptr_ topological_match_to_ topological_observers_
	arc( root_ ); // AtomOP
	arc( atom_pointer_ ); // AtomPointer2D
	arc( internal_coords_need_updating_ ); // _Bool
	arc( xyz_coords_need_updating_ ); // _Bool
	//arc( topological_match_to_ ); // AtomTreeCAP
	//arc( topological_observers_ ); // utility::vector1<AtomTreeCAP>
	arc( dof_changeset_ ); // AtomDOFChangeSet
	arc( external_coordinate_residues_changed_ ); // ResidueCoordinateChangeListOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::kinematics::AtomTree );
CEREAL_REGISTER_TYPE( core::kinematics::AtomTree )

CEREAL_REGISTER_DYNAMIC_INIT( core_kinematics_AtomTree )
#endif // SERIALIZATION
