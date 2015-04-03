// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

#ifndef INCLUDED_devel_cartesian_frags_CartesianFragment_hh
#define INCLUDED_devel_cartesian_frags_CartesianFragment_hh


// libRosetta headers

#include <devel/cartesian_frags/CartesianFragment.fwd.hh>
#include <devel/cartesian_frags/SafeID.hh>

#include <core/types.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh> // atomindices

#include <core/kinematics/RT.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/AtomTree.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <utility/vector1.hh>


namespace devel {
namespace cartesian_frags {


class CartesianFragment {

public:
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector Vector;
	typedef core::id::AtomID AtomID;
	typedef core::conformation::Conformation Conformation;
	typedef core::kinematics::RT RT;
	typedef core::kinematics::Stub Stub;
	typedef core::chemical::AtomIndices AtomIndices;
	typedef utility::vector1< Size > ComponentOffsets;
	typedef core::id::AtomID_Map< int > SeqposOffsetMap;

public:
	///////////////
	// construction
	///////////////

	CartesianFragment()
	{}


	/// c-tor -- only 1 outgoing, no new bonds
	CartesianFragment(
		TorsionStubID const & incoming,
		SafeAtomID const & root,
		TorsionStubID const & outgoing,
		Conformation const & conf
	)
	{
		utility::vector1< TorsionStubID > v;
		v.push_back( outgoing );
		initialize( incoming, root, v, utility::vector1< SafeBondID >(), conf );
	}


	/// c-tor -- >=1 outgoing, no new bonds
	CartesianFragment(
		TorsionStubID const & incoming,
		SafeAtomID const & root,
		utility::vector1< TorsionStubID > const & outgoing,
		Conformation const & conf
	)
	{
		initialize( incoming, root, outgoing, utility::vector1< SafeBondID >(), conf );
	}


	/// c-tor -- >= 1 outgoing, 1 new bond
 	CartesianFragment(
 		TorsionStubID const & incoming,
 		SafeAtomID const & root,
 		TorsionStubID const & outgoing,
 		SafeBondID const & new_bond,
 		Conformation const & conf
 	)
	{
		initialize( incoming, root, utility::vector1< TorsionStubID >( 1, outgoing), utility::vector1< SafeBondID >( 1, new_bond ), conf );
	}

 	CartesianFragment(
 		TorsionStubID const & incoming,
 		SafeAtomID const & root,
 		utility::vector1< TorsionStubID > const & outgoing,
 		utility::vector1< SafeBondID > const & new_bonds,
 		Conformation const & conf
 	)
 	{
 		initialize( incoming, root, outgoing, new_bonds, conf );
 	}


	/////////
	// access
	/////////


	RT const &
	stub_transform( int const i ) const
	{
		return stub_transforms_[i];
	}


	RT const &
	rt( int const i ) const
	{
		return stub_transforms_[i];
	}


	Size
	has_atom( SafeAtomID const & id ) const
	{
		return ( std::find( atom_ids_.begin(), atom_ids_.end(), id ) != atom_ids_.end() );
	}


	Size
	atom_index( SafeAtomID id ) const;


	Vector const &
	xyz( SafeAtomID const & id ) const
	{
		return xyz_[ atom_index( id ) ];
	}


	Size
	natoms() const
	{
		return xyz_.size();
	}


	////////////
	// insertion
	////////////

	/// @brief  Insert into a conformation (1-component frags only)
	void
	insert( Conformation & conf, Size const offset ) const;


	/// @brief  Insert into a conformation
	void
	insert( Conformation & conf, ComponentOffsets const & offsets ) const;


	///////////////
	// modification
	///////////////


	void
	set_torsion_angle(
										SafeAtomID const & id1,
										SafeAtomID const & id2,
										SafeAtomID const & id3,
										SafeAtomID const & id4,
										Real const setting
										);


	void
	set_bond_angle(
								 SafeAtomID const & id1,
								 SafeAtomID const & id2,
								 SafeAtomID const & id3,
								 Real const setting
								 );


	Size
	n_outgoing() const
	{
		return stub_transforms_.size();
	}


private: // private methods:


	Size
	add_atom( SafeAtomID const & id, Stub const & instub, Vector const & atom_xyz );

	int
	get_seqpos_offset_for_new_component( Size const component_root_seqpos, Size const component_index ) const;


	int
	get_seqpos_for_fragment_insertion( int const frag_seqpos, ComponentOffsets const & offsets ) const;


	int
	get_seqpos_offset_for_torsion_stub(
		TorsionStubID const & id,
		Conformation const & conf,
		SeqposOffsetMap const & seqpos_offset_map
	) const;


	void
	add_stub_atoms(
		TorsionStubID const & torstub_id,
		Stub const & instub,
		Conformation const & conf,
		SeqposOffsetMap const & seqpos_offset_map,
		std::string const & atom_name_prefix
	);

	/// @brief  c-tor helper
	void
	initialize(
		TorsionStubID const & incoming,
		SafeAtomID const & root,
		utility::vector1< TorsionStubID > const & outgoing,
		utility::vector1< SafeBondID > const & new_bonds_in,
		Conformation const & conf
	);


	/// @brief  recursive routine for adding atoms from a conformation, returns atomindex of newly added atom
	Size
	add_frag_atom(
		core::id::AtomID const & id,
		SeqposOffsetMap & seqpos_offset_map,
		core::kinematics::Stub const & instub,
		Conformation const & conf,
		utility::vector1< core::id::BondID > const & cuts,
		utility::vector1< core::id::BondID > const & new_bonds,
		core::id::AtomID_Mask & added
	);


	bool
	skip_bond( AtomID const & id1, AtomID const & id2, utility::vector1< core::id::BondID > const & cuts ) const
	{
		return ( std::find( cuts.begin(), cuts.end(), core::id::BondID( id1, id2 ) ) != cuts.end() ||
						 std::find( cuts.begin(), cuts.end(), core::id::BondID( id2, id1 ) ) != cuts.end() );
	}


	AtomID
	atom_tree_id( SafeAtomID const & id )
	{
		return AtomID( atom_index(id), 1 );
	}


	void
	update_xyz_from_tree();


	void
	setup_atom_tree();


private: // data

	// NOTE: all stubs have had their residue numbers adjusted by seqpos_offset_map[ frag atom in stub ]

	// incoming stub
	TorsionStubID incoming_;
	SafeStubID incoming_stub_id_;

	// outgoing stubs
	utility::vector1< TorsionStubID > outgoing_;
	utility::vector1< SafeStubID > outgoing_stub_ids_;

	// atom positions
	utility::vector1< Vector > xyz_;

	// atom ids
	utility::vector1< SafeAtomID > atom_ids_;

	// transforms between the incoming and outgoing stubs
	utility::vector1< RT > stub_transforms_;

	// subset of bonds used in building the frag -- in principal should form a tree
	utility::vector1< core::chemical::AtomIndices > links_;

	// private atomtree for kinematic updates, empty until we are asked to modify a dof
	core::kinematics::AtomTree tree_;


	Size ncomponents_;
	Size BIG_OFFSET_;

};

}
}

#endif
