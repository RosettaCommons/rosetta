// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/tree/BondedAtom.cc
/// @brief  Kinematics
/// @author Phil Bradley


// Unit headers
#include <core/kinematics/tree/BondedAtom.hh>

// Package headers
#include <core/id/DOF_ID.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.hh>

// Project headers
#include <core/id/AtomID_Map.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/kinematics/types.hh>
#include <core/types.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// C++ headers
#include <cassert>
#include <iostream>


namespace core {
namespace kinematics {
namespace tree {

/// @details Invokes Atom_ dfs function before, optionally recursing to younger siblings
/// should those younger siblings be effected by a DOF change on this node (e.g. phi_ change).
void
BondedAtom::dfs(
	AtomDOFChangeSet & changeset,
	ResidueCoordinateChangeList & res_change_list,
	Size const start_atom_index
) const
{

	Atom_::dfs( changeset, res_change_list, start_atom_index );
	if ( start_atom_index == dof_refold_index() ) {
		if ( dof_change_propagates_to_younger_siblings_ ) {
			AtomOP parent( parent_ ); // must have parent
			Atoms_Iterator iter = parent->atoms_begin();
			/// you had better find yourself in your parent's atom list.
			while ( (*iter).get() != this ) { ++iter; assert( iter != parent->atoms_end() );}
			++iter; // point to your next-youngest sibling.
			while ( iter != parent->atoms_end() ) {
				(*iter)->dfs( changeset, res_change_list, start_atom_index );
				++iter;
			}
		}
	}
}

/// @details Relies on get_input_stub, which will use the coordinates of this atom's ancestors to
/// build a stub.  If this function has been invoked by AtomTree update_xyz_coords(), then these
/// coordinates are guaranteed correct, since this atom is the root of a tree which needs to be
/// refolded. Ergo, nothing in the tree above this atom needs to be refolded.
void
BondedAtom::update_xyz_coords()
{


	/// dof_change_propagates_to_younger_siblings_ will be set to false inside
	/// update_xyz_coords -- keep a local copy.
	bool local_dof_change_propagates_to_younger_siblings( dof_change_propagates_to_younger_siblings_ );

	/// Ancestral coordinates are up-to-date since this node
	/// is the root of a subtree that needs refolding.
	/// The stub is passed to update_xyz_coords, and this atom will modify it;
	/// after which the stub is ready to be passed to the younger siblings.
	Stub stub( get_input_stub() );

	BondedAtom::update_xyz_coords( stub );

	if ( local_dof_change_propagates_to_younger_siblings ) {
		AtomOP parent( parent_ ); // must have parent
		Atoms_Iterator iter = parent->atoms_begin();
		/// you had better find yourself in your parent's atom list.
		while ( (*iter).get() != this ) { ++iter; assert( iter != parent->atoms_end() );}
		++iter; // point to your next-youngest sibling.
		while ( iter != parent->atoms_end() ) {
			(*iter)->update_xyz_coords( stub );
			++iter;
		}
	}
}


/////////////////////////////////////////////////////////////////////////////
/// @details starting from the input stub, calculate xyz position of this atom from
/// its internal coordinates d_, theta_ and phi_. If recusrvie is true,
/// obtain the new stub centered at this atom and pass the new stub to
/// all its children atoms to update their xyz positions recursively.
/// @note stub passed in is modified by rotating phi_ around x in the stub frame
void
BondedAtom::update_xyz_coords(
	Stub & stub
)
{
	using numeric::x_rotation_matrix_radians;
	using numeric::z_rotation_matrix_radians;
	using numeric::constants::d::pi;

	assert( stub.is_orthogonal( 1e-3 ) );

	stub.M *= x_rotation_matrix_radians( phi_ ); // this gets passed out

	Stub new_stub( stub.M * z_rotation_matrix_radians( theta_ ), stub.v );

	if ( std::abs( theta_ - pi ) < 1e-6 ) {
		// very special case
		if ( keep_dof_fixed( THETA ) ) {
			new_stub.M *= x_rotation_matrix_radians( pi );
		}
	}

	new_stub.v += d_ * new_stub.M.col_x();

	position( new_stub.v );

	for ( Atoms_Iterator it=atoms_begin(), it_end = atoms_end();
				it != it_end; ++it ) {
		(*it)->update_xyz_coords( new_stub );
	}

	/// Reset the output-sensitive refold information
	dof_change_propagates_to_younger_siblings_ = false;
	note_xyz_uptodate();
}


/////////////////////////////////////////////////////////////////////////////
/// @details starting from the input stub, calculate the internal coordinates d_, theta_ and
/// phi_ for this atom. If recursive is true, obtain the new stub centered at
/// this atom and pass the new stub to all its children atoms to update their
/// internal coordinates recursively.
/// @note stub passed in is modified by rotating phi_ around x in the stub frame
void
BondedAtom::update_internal_coords(
	Stub & stub,
	bool const recursive // = true
)
{
	using numeric::x_rotation_matrix_radians;
	using numeric::z_rotation_matrix_radians;
	using numeric::constants::d::pi;

	assert( stub.is_orthogonal( 1e-3 ) );

	numeric::xyzVector< core::Real > w( position() - stub.v );

	d_ = w.length();

	bool flip_stub( false );
	if ( d_ < 1e-2 ) {
		// phi, theta don't make much sense
		//	std::cerr << "WARNING:: very small d= " << d_ << ' ' << id() << std::endl;
		phi_ = 0.0;
		theta_ = 0.0;
	} else {
		//if ( d_ < 1e-1 ) {
			//std::cerr << "WARNING:: small d but we are calculating phi,theta: " << d_ << std::endl;
		//}
		w.normalize();
		Real const x( dot( w, stub.M.col_x() ) );
		Real const y( dot( w, stub.M.col_y() ) );
		Real const z( dot( w, stub.M.col_z() ) );

		Real const tol( 1e-6 );
		if ( x < -1.0 + tol ) {
			// very special case:
			// confirm that we are the stub_atom2 of a jump:
			if ( keep_dof_fixed( THETA ) ) {
				theta_ = pi;
				phi_ = 0.0;
				flip_stub = true; // very special case
			} else {
				theta_ = pi;
				phi_ = 0.0;
			}
		} else if ( x > 1.0 - tol ) {
			//std::cout << "WARNING:: update_internal_coords: exactly parallel? " << id() << std::endl;
			theta_ = 0.0;
			phi_ = 0.0;
		} else {
			theta_ = std::acos( x ); // DANGER
			//if ( theta_ < 1e-2 || pi - theta_ < 1e-2 ) {
				// less than 0.57 degrees
				//std::cout << "WARNING:: small theta but we are calculating phi: " <<
				//	theta_ << std::endl;
			//}
			phi_  = std::atan2( z, y );
		} // small theta
	} // small d

	stub.M *= x_rotation_matrix_radians( phi_ );


	if ( recursive ) {
		Stub new_stub( stub.M * z_rotation_matrix_radians( theta_ ), position() );

		if ( flip_stub ) {
			// special case if I'm stub_atom2 of my parent (who is a jump)
			new_stub.M *= x_rotation_matrix_radians( pi );
		}

		for ( Atoms_Iterator it=atoms_begin(), it_end = atoms_end();
					it != it_end; ++it ) {
			(*it)->update_internal_coords( new_stub );
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
/// @details pass out the stub we would pass out if we were actually updating
/// coords or DOFs, which is rotate the stub around its x axis by phi_.
void
BondedAtom::update_stub(
	Stub & stub
) const
{
	stub.M *= numeric::x_rotation_matrix_radians( phi_ ); // this gets passed out
}


/////////////////////////////////////////////////////////////////////////////
void
BondedAtom::set_dof(
	DOF_Type const type,
	Real const value
)
{
	if ( type == id::PHI ) {
		phi_ = value;
	} else if ( type == id::THETA ) {
		theta_ = value;
	} else if ( type == id::D ) {
		d_ = value;
	} else {
		std::cout << "bad torsion type for Atom: " << type << std::endl;
	}
}

/// @details calls set_dof non-polymorphically: assumption is that BondedAtom is not subclassed,
/// or, that if it is, that the derived class implements this overloaded set_dof function.
void
BondedAtom::set_dof(
	DOF_Type const type,
	Real const value,
	AtomDOFChangeSet & set
)
{
	BondedAtom::set_dof( type, value );
	Atom_::note_dof_change( set );
	dof_change_propagates_to_younger_siblings_ |= ( type == id::PHI );
}


/////////////////////////////////////////////////////////////////////////////
Real
BondedAtom::dof(
	DOF_Type const type
) const
{
	if ( type == PHI ) {
		return phi_;
	} else if ( type == THETA ) {
		return theta_;
	} else if ( type == D ) {
		return d_;
	} else {
		std::cout << "bad torsion type for Atom: " << type << std::endl;
		utility_exit();
	}
	return 0.0;
}

/////////////////////////////////////////////////////////////////////////////
/// @note  This will recursively clone all this atom's offspring atoms
/// @note  Assumes that atom_pointer has already been properly dimensioned

AtomOP
BondedAtom::clone( AtomAP parent_in, AtomPointer2D & atom_pointer ) const
{

	BondedAtomOP new_me( new BondedAtom( *this ) );

	atom_pointer[ id() ] = new_me; // handles memory management

	new_me->id( id() );
	new_me->parent( parent_in );

	// copy DOFs
	new_me->set_dof(   PHI, dof(   PHI ) );
	new_me->set_dof( THETA, dof( THETA ) );
	new_me->set_dof(     D, dof(     D ) );

	// copy coords
	new_me->position( position() );

	new_me->dof_change_propagates_to_younger_siblings_ = dof_change_propagates_to_younger_siblings_;

	// copy atoms
	for ( Atoms_ConstIterator a=atoms_begin(), a_end = atoms_end();
				a != a_end; ++a ) {
		new_me->append_atom( (*a)->clone( AtomAP(new_me) /*the new parent*/, atom_pointer ) );
	}

	return new_me;
}


/////////////////////////////////////////////////////////////////////////////
/// @details last torsion is the torsion( Phi for BondedAtom and RB for JumpAtom) of the
/// parent atom or previous bonded sibling. Since BondedAtom's PHI is dependent
/// on its previous sibling BondedAtom,this may modify last_torsion, if our bond
/// torsion angle is changing and this atom has other sibling atoms after itself.
/// recursively done all its offspring
void
BondedAtom::setup_min_map(
	DOF_ID & last_torsion,
	DOF_ID_Mask const & allow_move,
	MinimizerMapBase & min_map
) const
{

	DOF_ID phi_torsion  ( id(), PHI   );
	DOF_ID theta_torsion( id(), THETA );
	DOF_ID d_torsion    ( id(), D     );

	if ( allow_move[ phi_torsion ] && !keep_dof_fixed( PHI ) ) {
		min_map.add_torsion( phi_torsion, last_torsion );
		last_torsion = phi_torsion;
	}

	// no more changes to last_torsion from now on //

	DOF_ID last_torsion_local( last_torsion );

	if ( allow_move[ theta_torsion ] && !keep_dof_fixed( THETA ) ) {
		min_map.add_torsion( theta_torsion, last_torsion_local );
		last_torsion_local = theta_torsion;
	}

	if ( allow_move[ d_torsion ] && !keep_dof_fixed( D ) ) {
		min_map.add_torsion( d_torsion, last_torsion_local );
		last_torsion_local = d_torsion;
	}

	// add me to the min_map
	min_map.add_atom( id(), last_torsion_local );

	for ( Atoms_ConstIterator it=atoms_begin(), it_end = atoms_end();
				it != it_end; ++it ) {
		(*it)->setup_min_map( last_torsion_local, allow_move, min_map );
	}
}

/////////////////////////////////////////////////////////////////////////////
///@li axis is the unit vector along the rotation axis for that DOF(Eab).
///@li end_pos is the ending point of this unit vector(Vb).
///
/// @details consider simple case like, A->B->C->D and we want to know the change of D position
/// with respect to the rotation along B->C bond( dr/dphi). In this case, Eab is
/// the unit vector along B->C, end_positon is C, and dr/dphi is given as Eab x (D-C).\n
/// For B->C bond rotation is the DOF(PHI) of Atom D, so end_pos is input_stub.v which is
/// the position of Atom C and axis is input_stub.M.col(1) which is the unit vector pointing
/// from B to C.\n
/// For DOF(THETA) of Atom D, it is a rotation around a unit vector (ending at C)
/// which is perpendicular to B->C->D plane. So the end_pos is the position of Atom C and the
/// axix is the my_stub.M.col(3).
/// For DOF(D) of Atom D, it is a translation along C->D axis and that is my_stub.M.col(1)

void
BondedAtom::get_dof_axis_and_end_pos(
	Vector & axis,
	Position & end_pos,
	DOF_Type const type
) const
{
	Stub const my_stub( get_stub() );
	Stub const input_stub( get_input_stub() );

	if ( type == PHI ) {
		end_pos = input_stub.v;
		axis = input_stub.M.col(1);
	} else if ( type == THETA ) {
		end_pos = input_stub.v;
		axis = my_stub.M.col(3);
	} else if ( type == D ) {
		axis = my_stub.M.col(1);
	} else {
		std::cout << "Bad torsion type for Atom" << type << std::endl;
		utility_exit();
	}
}

/////////////////////////////////////////////////////////////////////////////
/// @details
/// - bond distance D can always be flexible
/// - bond angle THETA is not meaningful if it is defined by A-B-A in some
///   special cases, for example N-CA-C in which CA is a jump atom, then THETA
///   of the N atom is defined by N-CA-N.
/// - torsion angle PHI is not meaningful if it is defined by A-B-A-D or A-B-D-A
///   or A-D-B-A in some special cases.
bool
BondedAtom::keep_dof_fixed(
	DOF_Type const type
) const
{
	if ( type == D ) {
		return false;
	}

	if ( type == THETA ) {
		AtomCOP parent( parent_ ); // must have parent
		return ( parent->is_jump() && id() == parent->stub_atom2_id() );
	} else if ( type == PHI ) {
		AtomCOP parent( parent_ ); // must have parent
		if( ( parent->is_jump() &&
				( id() == parent->stub_atom2_id() ||
					id() == parent->stub_atom3_id() ) ) ) {
				return true;
		}

		AtomCOP parent_parent = parent->parent();
		if( parent_parent && parent_parent->is_jump() && id() == parent_parent->stub_atom3_id() ) {
			return true;
		}
		return false;
		
	} else {
		std::cout << "BondedAtom::keep_dof_fixed: BAD_TYPE: " <<type <<
			std::endl;
		assert( false );
		utility_exit();
	}
	return false;
}


///////////////////////////////////////////////////////////////////////////
/// @details copy DOFs, xyz's.
/// this asserts equal topology.
/// do recursively to copy for all its children
void
BondedAtom::copy_coords( Atom const & src )
{
	phi_   = src.dof(   PHI );
	theta_ = src.dof( THETA );
	d_     = src.dof(     D );

	// copy xyz as well as the dof_change_index_
	Super::operator= ( static_cast< Atom_ const & > ( src ));

	// check for topology mismatch
	assert( atom_id() == src.atom_id() && n_children() == src.n_children() );

	// call recursively for my children
	int i(0);
	for ( Atoms_Iterator it=atoms_begin(), it_end = atoms_end();
				it != it_end; ++it, ++i ) {
		(*it)->copy_coords( *(src.child(i)) );
	}
}

///
Jump BOGUS_JUMP;

} // tree
} // namespace kinematics
} // namespace core
