// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/tree/JumpAtom.cc
/// @brief  Jump atom
/// @author Phil Bradley


// Unit headers
#include <core/kinematics/tree/JumpAtom.hh>

// Package headers
#include <core/kinematics/MinimizerMapBase.hh>

// Project headers
#include <core/id/AtomID_Map.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace kinematics {
namespace tree {

/////////////////////////////////////////////////////////////////////////////
/// @details the root jump is not flexible (currently)
bool
JumpAtom::keep_dof_fixed(
	DOF_Type const //type
) const
{
	// the root jump is not flexible (currently)
	return ! parent();
}

/////////////////////////////////////////////////////////////////////////////
/// @note along n2c direction
Real
JumpAtom::dof(
	DOF_Type const type
) const
{
	int const n2c(1);
	int const rb_no( get_rb_number( type ) );
	if ( rb_no == 0 ) {
		std::cerr << "bad torsion type for JumpAtom: " << type << std::endl;
		utility_exit();
	}
	return jump_.get_rb_delta( rb_no, n2c );
}

/// @details Invokes Atom_ dfs function (which visits the subtree routed at this node).
/// No younger siblings are effected by a DOF change at a JumpAtom.
void
JumpAtom::dfs(
	AtomDOFChangeSet & changeset,
	ResidueCoordinateChangeList & res_change_list,
	Size const start_atom_index
) const
{
	Atom_::dfs( changeset, res_change_list, start_atom_index );
}

/////////////////////////////////////////////////////////////////////////////
/// @note along n2c direction
void
JumpAtom::set_dof(
	DOF_Type const type,
	Real const value
)
{
	debug_assert( parent().get() != nullptr );
	int const n2c(1);
	int const rb_no( get_rb_number( type ) );
	if ( rb_no == 0 ) {
		std::cerr << "bad torsion type for JumpAtom: " << type << std::endl;
		utility_exit();
	}
	jump_.set_rb_delta( rb_no, n2c, value );
}

/// @details calls set_dof non-polymorphically: assumption is that JumpAtom is not subclassed,
/// or, that if it is, that the derived class implements this overloaded set_dof function.
void
JumpAtom::set_dof(
	DOF_Type const type,
	Real const value,
	AtomDOFChangeSet & set
)
{
	JumpAtom::set_dof( type, value );
	Atom_::note_dof_change( set );
}

/////////////////////////////////////////////////////////////////////////////
Jump const &
JumpAtom::jump() const
{
	return jump_;
}

/////////////////////////////////////////////////////////////////////////////
void
JumpAtom::jump( Jump const & jump_in )
{
	jump_ = jump_in;
}

/////////////////////////////////////////////////////////////////////////////
void
JumpAtom::jump(
	Jump const & jump_in,
	AtomDOFChangeSet & set
)
{
	jump_ = jump_in; // could invoke 1-argument jump() method to avoid duplication...
	Atom_::note_dof_change( set );
}


/////////////////////////////////////////////////////////////////////////////
/// @details copy DOFs, xyz's.
/// this asserts equal topology.
/// do recursively to copy for all its children
void
JumpAtom::copy_coords(
	Atom const & src
)
{
	jump_ = src.jump();

	// copy xyz as well as the dof_change_index_
	Super::operator= ( static_cast< Atom_ const & > ( src ));

	if ( id() != src.id() || n_children() != src.n_children() ) {
		std::cerr << "JumpAtom:: copy_coords: topology_mismatch!" <<
			id()  << ' ' << src.id() << ' ' << n_children() << ' ' <<
			src.n_children() << std::endl;
		utility_exit();
	}

	int i(0);
	for ( auto it= atoms_begin(), it_end= atoms_end();
			it != it_end; ++it, ++i ) {
		(*it)->copy_coords( *(src.child(i) ) );
	}
}

/// @details Relies on get_input_stub, which will use the coordinates of this atom's ancestors to
/// build a stub.  If this function has been invoked by AtomTree update_xyz_coords(), then these
/// coordinates are guaranteed correct, since this atom is the root of a tree which needs to be
/// refolded. Ergo, nothing in the tree above this atom needs to be refolded.
void
JumpAtom::update_xyz_coords()
{
	Stub stub( get_input_stub() );

	// this code initializes refolding from an inverted jump
	if ( jump_.get_invert_upstream() ) {
		stub.M = stub.M * Jump::mirror_z_transform;
		//stub.v = stub.v; //?
	}

	JumpAtom::update_xyz_coords( stub );
}

/////////////////////////////////////////////////////////////////////////////
/// @details call make_jump to "jump" from parent to this atom.
/// Will recursively update xyz positions for all its offspring atoms.
/// @note the input stub is not changed
void
JumpAtom::update_xyz_coords(
	Stub & stub // in fact is const
)
{
	debug_assert( stub.is_orthogonal( 1e-3 ) );

	Stub new_stub;
	jump_.make_jump( stub, new_stub );
	position( new_stub.v );

	for ( auto it= atoms_begin(), it_end= atoms_end();
			it != it_end; ++it ) {
		(*it)->update_xyz_coords( new_stub );
	}
	note_xyz_uptodate();
}

/////////////////////////////////////////////////////////////////////////////
/// @details update the jump from the input stub and this atom's own stub. If defined,
/// will recursively update internal coords for all its offspring atoms.
/// @note The input stub is not changed
void
JumpAtom::update_internal_coords(
	Stub & stub,
	bool const recursive // = true
)
{
	debug_assert( stub.is_orthogonal( 1e-3 ) );

	Stub new_stub( get_stub() );


	// 8 cases dependending on: 1) stub handedness, 2) inv_up, and 3) inv_down
	// new_stub is always right-handed
	// jump (internally) is always right-handed
	// rules
	//  inv_stub inv_up inv_down   stub  new_stub  new_down
	//  ---------------------------------------------------
	//  0        0      0          no    no        no
	//  0        0      1          no    no        yes
	//  0        1      0          yes   yes       no
	//  0        1      1          yes   yes       yes
	//  1        0      0          yes   no        no
	//  1        0      1          yes   no        yes
	//  1        1      0          no    yes       no
	//  1        1      1          no    yes       yes
	// perhaps this logic should live in Jump?

	bool inv_stub = (stub.M.det()<0);
	bool inv_up=jump_.get_invert_upstream(), inv_down=jump_.get_invert_downstream();

	// it is safe to edit this in place since upstream is inverted for all
	// child atoms (since an atom is either left- or right- handed
	if ( inv_stub != inv_up ) {
		stub.M = stub.M * Jump::mirror_z_transform;
	}
	if ( inv_up ) {
		new_stub.M =  new_stub.M * Jump::mirror_z_transform;
	}
	jump_.from_stubs( stub, new_stub );

	if ( inv_up != inv_down ) {
		new_stub.M =  new_stub.M * Jump::mirror_z_transform;
	}

	if ( recursive ) {
		for ( auto it= atoms_begin(), it_end= atoms_end();
				it != it_end; ++it ) {
			(*it)->update_internal_coords( new_stub, recursive );
		}
	}
}

void
JumpAtom::steal_inversion(AtomOP steal_from) {
	JumpAtomOP steal_fromJ = utility::pointer::dynamic_pointer_cast< JumpAtom >(steal_from);
	std::pair<bool,bool> inv = steal_fromJ->get_inversion();
	jump_.set_invert(inv.first, inv.second);
}


/////////////////////////////////////////////////////////////////////////////
/// @note this will recursively clone all this atom's offspring atoms
AtomOP
JumpAtom::clone( AtomAP parent_in, AtomPointer2D & atom_pointer ) const
{

	JumpAtomOP new_me( new JumpAtom(*this) );

	atom_pointer[ id() ] = new_me;

	new_me->id( id() );
	new_me->parent( parent_in );

	// copy DOFs
	new_me->jump_ = jump_;

	// copy coords
	new_me->position( position() );

	// copy atoms
	for ( auto it= atoms_begin(), it_end= atoms_end();
			it != it_end; ++it ) {
		new_me->append_atom( (*it)->clone( AtomAP(new_me), atom_pointer ) );
	}

	return new_me;
}


/////////////////////////////////////////////////////////////////////////////
/// @details last torsion is the torsion( Phi for BondedAtom and RB for JumpAtom) of the
/// parent atom or the previous bonded sibling. Since unlike BondedAtom, JumpAtom's RB is independent from
/// other sibling atoms, this will not modify last_torsion, unlike Atom::setup_min_map.
/// recursively done all its offspring
void
JumpAtom::setup_min_map(
	DOF_ID & last_torsion, // const!!!
	DOF_ID_Mask const & allow_move,
	MinimizerMapBase & min_map
) const
{
	DOF_ID last_torsion_local( last_torsion );

	for ( int k=1; k<=6; ++k ) {
		DOF_Type const & type( id::get_rb_type(k) );
		DOF_ID rb_torsion( id(), type );
		if ( allow_move[ rb_torsion ] && !keep_dof_fixed( type ) ) {
			debug_assert( parent().get() != nullptr ); // root DOFs don't move
			min_map.add_torsion( rb_torsion, last_torsion_local );
			last_torsion_local = rb_torsion;
		}
	}

	// add me to the min_map
	min_map.add_atom( id(), last_torsion_local );

	for ( auto it= atoms_begin(), it_end= atoms_end();
			it != it_end; ++it ) {
		(*it)->setup_min_map( last_torsion_local, allow_move, min_map );
	}
}


/////////////////////////////////////////////////////////////////////////////
/// @li axis is the unit vector along the rotation axis for that DOF(Eab).
///2li end_pos is the ending point of this unit vector(Vb).
///
/// @details RB1, RB2 and RB3 are translation along x, y and z axis in the jump_start
///(input_stub) frame, which are input_stub.col(1), input_stub.col(2) and
///input_stub.col(3).\n
///RB4, RB5 and RB6 are rotations along x, y and z axis in the input_stub frame, and
///the rotation is applied around the point of "rb_center" written in the jump_end
///(my_stub) frame. So "end_pos"for these 3 DOFs is the rb_center rewritten in xyz
///frame, which is my_stub.V+my_stub.M*rb_center. The axis is not simply x, y and z
///in the input_stub because these 3 DOFs are not independently applied. Here are how
///they are derived (by me and there might be other smarter way to think of it):
///
/// @li the jump rotation matrix is R = Rz(RB6) * Ry(RB5) * Rx(RB4) * rt.rotation.
///
/// @li if RB6 is perturbed by "d", the new rotation matrix is R' = Rz(d+RB6) * Ry(RB5) * Rx(RB4) * rt.rotation
/// = Rz(d) * Rz(RB6) * Ry(RB5) * Rx(RB4) * rt.rotation = Rz(d) * R. This is to say perturbing RB6 by d is equivalent
/// to making an extra rotation around Z axis in the input_stub frame and therefore the axis is just input_stub.col(3).
///
/// @li if RB5 is perturbed by "d", the new rotation matrix is R' = Rz(RB6) * Ry(d+RB5) * Rx(RB4) * rt.rotation
/// = Rz(RB6) * Ry(d) * Ry(RB5) * Rx(RB4) * rt.rotation = Rz(RB6) * Ry(d) * Rz(RB6)^ * R. This is to say perturbing
/// RB5 by d is equivalent to making an extra rotation Rz(RB6) * Ry(d) * Rz(RB6)^ in the input_stub frame and
/// the axis of rotation is the one we try to get. what is it? Remember this rotation matrix is in the context of
/// the input_stub frame M, so rewritting it in the xyz lab frame results in M * Rz(RB6) * Ry(d) * Rz(RB6)^ * M^. Since
/// for matrices (A*B)^ = B^ * A^, we get (M *Rz(RB6)) * Ry(d) * (M * Rz(RB6))^ and this is eqivanelent to making
/// a rotation around y axis in the frame of (M *Rz(RB6)) and therefore the axis is (input_stub.M * Z_rotation ).col(2)
///
/// @li similarly, one can the axis for a perturbation to RB4 DOF.
void
JumpAtom::get_dof_axis_and_end_pos(
	Vector & axis,
	Position & end_pos,
	DOF_Type const type
) const
{
	using numeric::y_rotation_matrix_degrees;
	using numeric::z_rotation_matrix_degrees;

	Stub input_stub( get_input_stub() );

	// special case #1: inverted upstream jump
	if ( jump_.get_invert_upstream() ) {
		assert (input_stub.M.det()>0);
		input_stub.M = input_stub.M * Jump::mirror_z_transform;
	}

	int const rb_no( get_rb_number( type ) );
	if ( rb_no <= 3 ) {
		axis = input_stub.M.col( rb_no );
	} else {
		Stub my_stub( get_stub() );

		int const n2c( 1 );
		utility::vector1< Real > const & rb_delta( jump_.get_rb_delta( n2c ) );
		numeric::xyzVector< Real > const & rb_center( jump_.get_rb_center( n2c ));

		// special case #2: inverted downstream jump
		//   recall the rotation parameters describe the uninverted jump
		//   so we work with an uninverted input_stub
		if ( jump_.get_invert_downstream() ) {
			assert (my_stub.M.det()>0);

			my_stub.M = my_stub.M * Jump::mirror_z_transform;
			end_pos = my_stub.v + my_stub.M * rb_center;

			if ( type == id::RB6 ) {
				axis = -1.0 * (input_stub.M).col(3);
			} else if ( type == id::RB5 ) {
				Real const theta_z = rb_delta[6];
				axis = -1.0 * ( input_stub.M * z_rotation_matrix_degrees( theta_z ) ).col(2);
			} else if ( type == id::RB4 ) {
				Real const theta_z = rb_delta[6], theta_y = rb_delta[5];
				axis = -1.0 * ( input_stub.M * z_rotation_matrix_degrees( theta_z ) *
					y_rotation_matrix_degrees( theta_y ) ).col(1);
			} else {
				std::cerr << "Bad torsion type for Atom" << type << std::endl;
				utility_exit();
			}
		} else {
			end_pos = my_stub.v + my_stub.M * rb_center;

			if ( type == id::RB6 ) {
				axis = input_stub.M.col(3);
			} else if ( type == id::RB5 ) {
				Real const theta_z = rb_delta[6];
				axis = ( input_stub.M * z_rotation_matrix_degrees( theta_z ) ).col(2);
			} else if ( type == id::RB4 ) {
				Real const theta_z = rb_delta[6], theta_y = rb_delta[5];
				axis = ( input_stub.M * z_rotation_matrix_degrees( theta_z ) *
					y_rotation_matrix_degrees( theta_y ) ).col(1);
			} else {
				std::cerr << "Bad torsion type for Atom" << type << std::endl;
				utility_exit();
			}
		}
	}
}

Atom const *
JumpAtom::raw_stub_atom1() const
{
	return ( stub_defined() ? this : raw_parent() );
}

Atom const *
JumpAtom::raw_stub_atom2() const
{
	if ( stub_defined() ) {
		return raw_get_nonjump_atom(0);
	}
	Atom const * parent_ptr = raw_parent();
	if ( parent_ptr ) {
		return parent_ptr->raw_stub_atom2();
	}
	return nullptr;
}

Atom const *
JumpAtom::raw_stub_atom3() const
{
	//std::cout << "stub_atom3: " << this << ' ' << parent_ << std::endl();
	if ( stub_defined() ) {
		Atom const * first( raw_get_nonjump_atom(0) );
		Atom const * second( first->raw_get_nonjump_atom(0) );
		if ( second ) {
			return second;
		} else {
			return raw_get_nonjump_atom(1);
		}
	}
	Atom const * parent_ptr = raw_parent();
	if ( parent_ptr ) {
		return parent_ptr->raw_stub_atom3();
	}
	return nullptr;
}


} // namespace tree
} // namespace kinematics
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::kinematics::tree::JumpAtom::save( Archive & arc ) const {
	arc( cereal::base_class< Atom_ >( this ) );
	arc( CEREAL_NVP( jump_ ) ); // class core::kinematics::Jump
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::kinematics::tree::JumpAtom::load( Archive & arc ) {
	arc( cereal::base_class< Atom_ >( this ) );
	arc( jump_ ); // class core::kinematics::Jump
}

SAVE_AND_LOAD_SERIALIZABLE( core::kinematics::tree::JumpAtom );
CEREAL_REGISTER_TYPE( core::kinematics::tree::JumpAtom )

CEREAL_REGISTER_DYNAMIC_INIT( core_kinematics_tree_JumpAtom )
#endif // SERIALIZATION
