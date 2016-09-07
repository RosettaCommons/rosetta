// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
#include <devel/cartesian_frags/CartesianFragment.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/kinematics/util.hh>

#include <core/id/AtomID_Map.hh>

#include <basic/Tracer.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

#include <core/pose/util.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

namespace devel {
namespace cartesian_frags {

using namespace core;
using namespace ObjexxFCL;
using utility::vector1;
using std::string;

/// tracer object:
static THREAD_LOCAL basic::Tracer tt( "devel.cartesian_frags.CartesianFragment", basic::t_trace );


/// @details  Private helper routine for adding an atom, updates required data. Returns atom-index of new atom.

core::Size
CartesianFragment::add_atom( SafeAtomID const & id, Stub const & instub, Vector const & atom_xyz )
{
	xyz_.push_back( instub.global2local( atom_xyz ) );
	atom_ids_.push_back( id );
	links_.push_back( AtomIndices() );
	return xyz_.size();
}


/// @details  Returns the seqpos_offset value for a new connected component of the fragment
/// This is used in converting from AtomID's in the source conformation to our own internal SafeAtomID's.
/// Private.
int
CartesianFragment::get_seqpos_offset_for_new_component(
	Size const component_root_seqpos,
	Size const component_index
) const
{
	return component_root_seqpos - ( component_index - 1 ) * BIG_OFFSET_;
}


/// @details  Maps from a frag-internal sequence number to the sequence position in the
///   conformation into which we are being inserted. Uses the offsets data which tells which
///   sequence positions are the roots for the different connected components of this fragment.

int
CartesianFragment::get_seqpos_for_fragment_insertion( int const frag_seqpos, ComponentOffsets const & offsets ) const
{
	// frag_seqpos = old_conf_seqpos - seqpos_offset[old_conf_id]
	//             = old_conf_seqpos - component_root_seqpos + component_index * BIG_OFFSET_
	//             = component_local_seqpos                  + component_index * BIG_OFFSET_
	//
	// in this function, we want to return:
	//
	// new_conf_seqpos = old_conf_seqpos - component_root_seqpos + offsets[ component ]
	//                 = component_local_seqpos                  + offsets[ component ]
	//
	// so the idea is that we are replacing component_root_seqpos with offsets[ component ]


	int const component_index( 1 + ( ( frag_seqpos + BIG_OFFSET_/2 ) / BIG_OFFSET_ ) );
	int const component_local_seqpos( frag_seqpos - ( component_index - 1 ) * BIG_OFFSET_ );
	assert( component_index >= 1 && component_index <= int(ncomponents_) );
	int const new_conf_seqpos( component_local_seqpos + offsets[ component_index ] );
	assert( new_conf_seqpos - get_seqpos_offset_for_new_component( offsets[ component_index ], component_index ) ==
		frag_seqpos );
	return new_conf_seqpos;
}


/// @details  Given a TorsionStubID corresponding to one of the bonds that bounds a fragment, return
///   the seqpos_offset of the atom in this bond which belongs to the fragment.
///
int
CartesianFragment::get_seqpos_offset_for_torsion_stub(
	TorsionStubID const & id,
	Conformation const & conf,
	SeqposOffsetMap const & seqpos_offset_map
) const
{
	id::BondID const bond( torsion_bond( id, conf ) );
	if ( seqpos_offset_map[ bond.atom1 ] != 0 ) {
		assert( seqpos_offset_map[ bond.atom2 ] == 0 );
		return seqpos_offset_map[ bond.atom1 ];
	} else {
		assert( seqpos_offset_map[ bond.atom2 ] != 0 );
		assert( seqpos_offset_map[ bond.atom1 ] == 0 );
		return seqpos_offset_map[ bond.atom2 ];
	}
}


/// @details  Insert myself into a conformation (this version only works for single-component fragments)
void
CartesianFragment::insert( Conformation & conf, Size const offset ) const
{
	insert( conf, ComponentOffsets( 1, offset ) );
}


/// @details  Insert myself into a conformation, using the offsets to get sequence numbers for my various
///   connected components

void
CartesianFragment::insert( Conformation & conf, ComponentOffsets const & offsets ) const
{
	assert( offsets.size() == ncomponents_ );

	kinematics::AtomTree::FragXYZ frag_xyz;
	kinematics::AtomTree::FragRT frag_rt;


	// frag atom xyz's
	for ( Size i=1; i<= natoms(); ++i ) {
		SafeAtomID id( atom_ids_[i] );
		if ( id.name.substr(0,6) == "instub" || id.name.substr(0,7) == "outstub" ) continue;
		id.seqpos = get_seqpos_for_fragment_insertion( id.seqpos, offsets );
		frag_xyz[ id.id( conf ) ] = xyz_[i];
	}

	tt << "Frag::insert instub_id= " << incoming_stub_id_ <<
		" nout= " << outgoing_stub_ids_.size() << " natoms= " << natoms() << std::endl;


	SafeStubID frag_incoming( incoming_stub_id_ );
	frag_incoming.atom1.seqpos = get_seqpos_for_fragment_insertion( frag_incoming.atom1.seqpos, offsets );
	frag_incoming.atom2.seqpos = get_seqpos_for_fragment_insertion( frag_incoming.atom2.seqpos, offsets );
	frag_incoming.atom3.seqpos = get_seqpos_for_fragment_insertion( frag_incoming.atom3.seqpos, offsets );

	{ // check if the incoming torsionstubid matches
		TorsionStubID tmp_tor( incoming_ );
		tmp_tor.id.rsd() = get_seqpos_for_fragment_insertion( tmp_tor.id.rsd(), offsets );
		SafeStubID tmp_safe( tmp_tor, conf );
		if ( tmp_safe != frag_incoming ) {
			tt << "Frag::insert: warning -- mismatch between incoming torsionstubid and incoming safestubid, " <<
				"using id from torsionstubid " << tmp_safe << ' ' << frag_incoming << std::endl;
			frag_incoming = tmp_safe;
		}
	}

	// stub transforms
	for ( Size i=1; i<= outgoing_.size(); ++i ) {
		tt << "Frag::insert outstub_id= " << outgoing_stub_ids_[i] << std::endl;

		// check if the torsionstubid matches
		SafeStubID frag_outgoing( outgoing_stub_ids_[i] );
		frag_outgoing.atom1.seqpos = get_seqpos_for_fragment_insertion( frag_outgoing.atom1.seqpos, offsets );
		frag_outgoing.atom2.seqpos = get_seqpos_for_fragment_insertion( frag_outgoing.atom2.seqpos, offsets );
		frag_outgoing.atom3.seqpos = get_seqpos_for_fragment_insertion( frag_outgoing.atom3.seqpos, offsets );
		{
			TorsionStubID tmp_tor( outgoing_[i] );
			tmp_tor.id.rsd() = get_seqpos_for_fragment_insertion( tmp_tor.id.rsd(), offsets );
			SafeStubID const tmp_safe( tmp_tor, conf );
			if ( tmp_safe != frag_outgoing ) {
				tt << "Frag::insert: warning -- mismatch between outgoing torsionstubid and safestubid, " <<
					"using id from torsionstubid " << tmp_safe << ' ' << frag_outgoing << std::endl;
				frag_outgoing = tmp_safe;
			}
		}

		id::StubID const id( frag_outgoing.id( conf ) );
		frag_rt[ id ] = stub_transforms_[i];
	}

	id::StubID instub_id( frag_incoming.id( conf ) );
	conf.insert_fragment( instub_id, frag_rt, frag_xyz );
}


/// @details  (private)  Given one of my boundary bond torsions, add the stub atoms that don't already belong to me.
///   Give them names like "instub_atom2" or "outstub3_atom1"
///   I don't recall all the reasons for doing this, but one is so that when we update xyz's using an atomtree
///   after a call to set_bond_angle or set_torsion_angle, we can update the in-out stub transforms by recalculating
///   the outgoing stubs (since we have xyz's for all three atoms of each)
///
void
CartesianFragment::add_stub_atoms(
	TorsionStubID const & torstub_id,
	Stub const & instub,
	Conformation const & conf,
	SeqposOffsetMap const & seqpos_offset_map,
	std::string const & atom_name_prefix
)
{
	int const seqpos_offset( get_seqpos_offset_for_torsion_stub( torstub_id, conf, seqpos_offset_map ) );

	Stub const stub( torsion_stub( torstub_id, conf ) );
	id::BondID const bond( torsion_bond( torstub_id, conf ) );
	SafeAtomID const bond_atom1( bond.atom1, conf, seqpos_offset ); // is stub_atom1 of torsion_stub
	SafeAtomID const bond_atom2( bond.atom2, conf, seqpos_offset ); // is stub_atom2 of torsion_stub

	if ( has_atom( bond_atom1 ) ) {
		assert( !has_atom( bond_atom2 ) );

		// need to add stub_atoms 2 and 3
		SafeAtomID const atom1( bond_atom1 );
		SafeAtomID const atom2( atom_name_prefix + "atom2", 0 );
		SafeAtomID const atom3( atom_name_prefix + "atom3", 0 );

		Size const atom1_index( atom_index( atom1 ) );
		Size const atom2_index( add_atom( atom2, instub, stub.build_fake_xyz(2) ) );
		links_[ atom1_index ].insert( links_[ atom1_index ].begin(), atom2_index );

		Size const atom3_index( add_atom( atom3, instub, stub.build_fake_xyz(3) ) );
		links_[ atom2_index ].push_back( atom3_index );

	} else {
		assert(  has_atom( bond_atom2 ) );
		assert( !has_atom( bond_atom1 ) );

		// only need to add stub_atom 1
		SafeAtomID const atom2( bond_atom2 );
		SafeAtomID const atom1( atom_name_prefix + "atom1", 0 );

		Size const atom2_index( atom_index( atom2 ) );
		Size const atom1_index( add_atom( atom1, instub, stub.build_fake_xyz(1) ) );
		links_[ atom2_index ].insert( links_[ atom2_index ].begin(), atom1_index );
	}
}


/// @details  c-tor helper function that collects code common to the various constructors.
///
void
CartesianFragment::initialize(
	TorsionStubID const & incoming,
	SafeAtomID const & root,
	vector1< TorsionStubID > const & outgoing,
	vector1< SafeBondID > const & new_bonds_in,
	conformation::Conformation const & conf
)
{
	BIG_OFFSET_ = conf.size() * 10;
	ncomponents_ = new_bonds_in.size() + 1;


	// get incoming stub
	Stub const instub( torsion_stub( incoming, conf ) );

	// convert new bonds to BondID's for this conformation
	vector1< id::BondID > new_bonds;
	for ( Size i=1; i<= new_bonds_in.size(); ++i ) {
		new_bonds.push_back( new_bonds_in[i].id( conf ) );
	}

	// calculate the transforms between the incoming stub and the outgoing stubs
	// and setup list of bonds not to cross
	vector1< id::BondID > cuts;
	cuts.push_back( torsion_bond( incoming, conf ) );
	for ( Size i=1; i<= outgoing.size(); ++i ) {
		// store the transform from incoming stub to this outgoing stub
		Stub const outstub( torsion_stub( outgoing[i], conf ) );
		stub_transforms_.push_back( RT( instub, outstub ) );

		// record this bond as one of the boundaries of the fragment
		cuts.push_back( torsion_bond( outgoing[i], conf ) );
	}


	// recursively collect the atoms needed, store in local coord sys defined by incoming stub
	id::AtomID_Mask added;
	core::pose::initialize_atomid_map( added, conf );
	xyz_.clear();
	atom_ids_.clear();

	SeqposOffsetMap seqpos_offset_map( 0 );
	core::pose::initialize_atomid_map( seqpos_offset_map, conf );

	seqpos_offset_map[ root.id( conf ) ] = get_seqpos_offset_for_new_component( root.seqpos, 1 /*1st component*/);
	add_frag_atom( root.id( conf ), seqpos_offset_map, instub, conf, cuts, new_bonds, added );


	//nreal_atoms_ = xyz_.size();
	add_stub_atoms( incoming, instub, conf, seqpos_offset_map, "instub_" );
	for ( Size i=1; i<= outgoing.size(); ++i ) {
		string const prefix( "outstub"+string_of(i)+"_" );
		add_stub_atoms( outgoing[i], instub, conf, seqpos_offset_map, prefix );
	}


	//// store all the miscellaneous stub ids ///////////////////////////////////////////////////////
	int const incoming_seqpos_offset( get_seqpos_offset_for_torsion_stub( incoming, conf, seqpos_offset_map ) );
	incoming_stub_id_ = SafeStubID( incoming, conf, incoming_seqpos_offset );
	incoming_ = incoming;
	incoming_.id.rsd() -= incoming_seqpos_offset;

	outgoing_ = outgoing;
	for ( Size i=1; i<= outgoing.size(); ++i ) {
		int const outgoing_seqpos_offset( get_seqpos_offset_for_torsion_stub( outgoing[i], conf, seqpos_offset_map ) );
		outgoing_[i].id.rsd() -= outgoing_seqpos_offset;
		outgoing_stub_ids_.push_back( SafeStubID( outgoing[i], conf, outgoing_seqpos_offset ) );
	}


} // end of initialize(...)


/// @details  (private)  Recursive routine that adds all the fragment atoms defined by a conformation and a set
///   of bonds to cut and extra connections to jump across.
/// returns atomindex of new atom

Size
CartesianFragment::add_frag_atom(
	id::AtomID const & id,
	SeqposOffsetMap & seqpos_offset_map,
	kinematics::Stub const & instub,
	conformation::Conformation const & conf,
	vector1< id::BondID > const & cuts,
	vector1< id::BondID > const & new_bonds,
	id::AtomID_Mask & added
)
{
	using namespace conformation;
	using namespace chemical;

	assert( !added[ id ] );
	assert( seqpos_offset_map[ id ] != 0 );

	// now mark this id as added
	added[ id ] = true;

	// some local variables
	Size const seqpos( id.rsd() );
	Size const atomno( id.atomno() );
	Residue const & rsd( conf.residue( seqpos ) );
	ResidueType const & rsd_type( rsd.type() );

	// store the coordinates, using a local atomid
	Size const my_atom_index( add_atom( SafeAtomID( id, conf, seqpos_offset_map[id]), instub, rsd.xyz( atomno ) ) );


	//// now recursively add neighbors:

	// first any new bonds, ie jumps for example, or nonpolymeric chemical bonds like disulfides
	//
	for ( Size i=1; i<= new_bonds.size(); ++i ) {
		if ( new_bonds[i].has( id ) ) {
			AtomID const & other_id( new_bonds[i].other_atom( id ) );
			if ( added[ other_id ] ) continue;
			assert( !skip_bond( id, other_id, cuts ) ); // that would be weird
			seqpos_offset_map[ other_id ] = get_seqpos_offset_for_new_component( other_id.rsd(), i+1 /* the new component */);
			links_[ my_atom_index ].push_back( add_frag_atom( other_id, seqpos_offset_map,
				instub, conf, cuts, new_bonds, added ) );
		}
	}


	// first the inter-residue nbrs: //////////////////////////
	{
		Size const ncon( rsd.n_possible_residue_connections() );
		for ( Size connid=1; connid <= ncon; ++connid ) {
			if ( rsd.residue_connect_atom_index( connid ) != atomno ) continue;
			// we don't cross non-polymeric connections, since the sequence numbers on the other
			// side would not necessarily be conserved when we insert. So we require user to add such a connection
			// explicitly in the new_bonds array, and thus to specify the additional sequence offset when inserting
			//
			if ( !rsd_type.residue_connection_is_polymeric( connid ) ) continue;
			ResConnID const & resconnid( rsd.actual_residue_connection( connid ) );
			Size const other_seqpos( resconnid. resid() );
			Size const other_connid( resconnid.connid() );
			Size const other_atomno( conf.residue( other_seqpos ).residue_connect_atom_index( other_connid ) );
			AtomID const other_id( other_atomno, other_seqpos );
			if ( added[ other_id ] ) continue;
			if ( skip_bond( id, other_id, cuts ) ) continue;
			seqpos_offset_map[ other_id ] = seqpos_offset_map[ id ];
			links_[ my_atom_index ].push_back( add_frag_atom( other_id, seqpos_offset_map, instub, conf, cuts,
				new_bonds, added ) );
		}
	} // scope

	// now the intra-residue nbrs: ////////////////////////
	{
		AtomIndices const & nbrs( rsd_type.nbrs( atomno ) );
		for ( Size i=1; i<= nbrs.size(); ++i ) {
			AtomID const nbr_id( nbrs[i], seqpos );
			if ( added[ nbr_id ] ) continue;
			if ( skip_bond( id, nbr_id, cuts ) ) continue;
			seqpos_offset_map[ nbr_id ] = seqpos_offset_map[ id ];
			links_[ my_atom_index ].push_back( add_frag_atom( nbr_id, seqpos_offset_map, instub, conf, cuts,
				new_bonds, added ) );
		}
	} // scope


	return my_atom_index;
}


/// @details  Get the atomindex for the atom with the desired safeatomid
///
Size
CartesianFragment::atom_index( SafeAtomID id ) const
{
	auto it( std::find( atom_ids_.begin(), atom_ids_.end(), id ) );
	if ( it == atom_ids_.end() ) {
		std::string const name( id.name );
		// look for stub
		if ( name.substr(0,6) == "instub" ) {
			// "instub_atom2"
			//  012345678901
			int const atomno = int_of( name.substr( 11, 1 ) );
			id = incoming_stub_id_.atom( atomno );
		} else if ( name.substr(0,7) == "outstub" ) {
			// "outstub3_atom2"
			//  01234567890123
			int const stubno = int_of( name.substr(  7, 1 ) );
			int const atomno = int_of( name.substr( 13, 1 ) );
			id = outgoing_stub_ids_[ stubno ].atom( atomno );
		}
		it = std::find( atom_ids_.begin(), atom_ids_.end(), id );
		if ( it == atom_ids_.end() ) utility_exit_with_message( "unknown frag atom" );
	}

	return it - atom_ids_.begin() + 1;
}


/// @details  Set a torsion angle, used for updating a fragment with new torsions derived from some optimization

void
CartesianFragment::set_torsion_angle(
	SafeAtomID const & id1,
	SafeAtomID const & id2,
	SafeAtomID const & id3,
	SafeAtomID const & id4,
	Real const setting
)
{
	if ( tree_.empty() ) setup_atom_tree();
	tree_.set_torsion_angle( atom_tree_id( id1 ),
		atom_tree_id( id2 ),
		atom_tree_id( id3 ),
		atom_tree_id( id4 ),
		setting
	);
	update_xyz_from_tree(); // silly: refold every time
}

/// @details  Set a bond angle, used for updating a fragment with new angles derived from some optimization

void
CartesianFragment::set_bond_angle(
	SafeAtomID const & id1,
	SafeAtomID const & id2,
	SafeAtomID const & id3,
	Real const setting
)
{
	if ( tree_.empty() ) setup_atom_tree();
	tree_.set_bond_angle( atom_tree_id( id1 ),
		atom_tree_id( id2 ),
		atom_tree_id( id3 ),
		setting
	);
	update_xyz_from_tree(); // silly: refold every time
}

/// @details  Update the xyz's using our atomtree in response to a bond angle/torsion change
void
CartesianFragment::update_xyz_from_tree()
{

	for ( Size i=1; i<= natoms(); ++i ) {
		xyz_[i] = tree_.xyz( AtomID( i,1 ) );
	}

	// now have to update the outgoing stub transforms
	for ( Size i=1; i<= n_outgoing(); ++i ) {
		std::string const prefix( "outstub"+string_of(i)+"_atom" );
		Stub const outstub( xyz( SafeAtomID( prefix+string_of(1), 0 ) ),
			xyz( SafeAtomID( prefix+string_of(2), 0 ) ),
			xyz( SafeAtomID( prefix+string_of(3), 0 ) ) );
		stub_transforms_[i] = RT( Stub(), outstub );
	}
}


/// @details  Setup an atomtree to manage our xyz's in response to dof changes
void
CartesianFragment::setup_atom_tree()
{
	if ( tree_.size() == 1 ) return;
	assert( tree_.size() == 0 );

	kinematics::AtomPointer2D atom_pointer( 1 );
	atom_pointer[1].resize( natoms() );

	int const root_atom(1);
	int const seqpos(1);
	kinematics::add_atom( root_atom, seqpos, links_, atom_pointer[1], true );

	// set coords
	for ( Size i=1; i<= natoms(); ++i ) {
		atom_pointer[1][i]->xyz( xyz_[i] );
	}

	tree_.replace_tree( atom_pointer );

	assert( tree_.size() == 1 );
}

} // namespace cartesian_frags
} // namespace devel

