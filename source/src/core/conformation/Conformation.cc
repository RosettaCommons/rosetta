// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/Conformation.cc
/// @brief  Method definitions for the Conformation class.
/// @author Phil Bradley
/// @author Modified by Vikram K. Mulligan (vmullig@uw.edu) to try to start to remove the assumption that all polymer connections are to the i+1 and i-1 residues.


// Unit headers
#include <core/conformation/Conformation.hh>

// Package headers
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/PseudoBond.hh>
#include <core/conformation/util.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.hh>

// Project headers
#include <core/id/AtomID_Map.hh>
#include <core/id/TorsionID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/constants.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/conformation/membrane/MembraneInfo.fwd.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.fwd.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>
#include <basic/basic.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/broker.OptionKeys.gen.hh>

// Utility headers
#include <utility/assert.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <boost/foreach.hpp>

// C++ headers
#include <algorithm>
#include <cassert>
#include <set>
//#include <stdio.h>


using basic::T;


static basic::Tracer TR("core.conformation.Conformation");


namespace core {
namespace conformation {

using namespace ObjexxFCL;

// Standard class methods ////////////////////////////////////////////////////////////////////////////////////////////

// default destructor
Conformation::~Conformation()
{
	clear();
	notify_connection_obs( ConnectionEvent( this, ConnectionEvent::DISCONNECT ) );
}

Conformation::Conformation() :
	utility::pointer::ReferenceCount(),
	fold_tree_( new FoldTree ),
	atom_tree_( new AtomTree ),
	contains_carbohydrate_residues_( false) ,
	residue_coordinates_need_updating_( false ),
	residue_torsions_need_updating_( false ),
	structure_moved_( true )
{
	atom_tree_->set_weak_pointer_to_self( atom_tree_() );
}

// copy constructor
Conformation::Conformation( Conformation const & src ) :
	utility::pointer::ReferenceCount()
{
	basic::ProfileThis doit( basic::CONFORMATION_COPY );
	// residues
	for ( Size i=1; i<= src.size(); ++i ) {
		residues_.push_back( src.residues_[i]->clone() );
	}

	// kinematics
	fold_tree_ = new FoldTree( *src.fold_tree_ );
	atom_tree_ = new AtomTree( *src.atom_tree_ );
	atom_tree_->set_weak_pointer_to_self( atom_tree_() );

	// carbohydrates?
	contains_carbohydrate_residues_ = src.contains_carbohydrate_residues_;

	// chain info
	chain_endings_ = src.chain_endings_;
	// secstruct
	secstruct_ = src.secstruct_;

	// bookkeeping
	residue_coordinates_need_updating_ = src.residue_coordinates_need_updating_;
	residue_torsions_need_updating_ = src.residue_torsions_need_updating_;

	dof_moved_ = src.dof_moved_;
	xyz_moved_ = src.xyz_moved_;

	structure_moved_ = src.structure_moved_;

}

// equals operator
/// @details If lengths & sequence of source and target are different, will fire a
///  LengthEvent::INVALIDATE signal and then an IdentityEvent::INVALIDATE signal.
Conformation &
Conformation::operator=( Conformation const & src )
{
	if ( &src == this ) return *this;

	basic::ProfileThis doit( basic::CONFORMATION_COPY );

	if ( src.size() == size() && sequence_matches( src ) ) {
		in_place_copy( src );
	} else {

		// delete current data
		clear();

		// residues
		for ( Size i=1; i<= src.size(); ++i ) {
			residues_.push_back( src.residues_[i]->clone() );
		}

		// kinematics
		(*fold_tree_) = (*src.fold_tree_);
		(*atom_tree_) = (*src.atom_tree_);

		// carbohydrates?
		contains_carbohydrate_residues_ = src.contains_carbohydrate_residues_;

		// chain info
		chain_endings_ = src.chain_endings_;
		// secstruct
		secstruct_ = src.secstruct_;

		// bookkeeping
		residue_coordinates_need_updating_ = src.residue_coordinates_need_updating_;
		residue_torsions_need_updating_ = src.residue_torsions_need_updating_;

		dof_moved_ = src.dof_moved_;
		xyz_moved_ = src.xyz_moved_;

		structure_moved_ = src.structure_moved_;

		// length may have radically changed, tell length observers to invalidate their data
		notify_length_obs( LengthEvent( this, LengthEvent::INVALIDATE, 0, 0, NULL ), false );
		// identity may have radically changed, tell length observers to invalidate their data
		notify_identity_obs( IdentityEvent( this, IdentityEvent::INVALIDATE, 0, NULL ), false );
	}

	notify_xyz_obs( XYZEvent( this ) );

	return *this;
}

///@details make a copy of this conformation( allocate actual memory for it )
ConformationOP
Conformation::clone() const
{
  return new Conformation( *this );
}

// clear data
void
Conformation::clear()
{
	pre_nresidue_change(); // nuke old atom tree update-data before destroying coordinates.

	fold_tree_->clear();
	atom_tree_->clear();
	residues_.clear();
	dof_moved_.clear();
	xyz_moved_.clear();
	chain_endings_.clear();
}


// Debugging /////////////////////////////////////////////////////////////////////////////////////////////////////////

/// @details check that the residue torsions are in sync with the residue coords and atomtree coords
void
Conformation::debug_residue_torsions( bool verbose ) const
{
	using namespace ObjexxFCL::format;
	using basic::subtract_degree_angles;

	update_residue_coordinates();
	update_residue_torsions();
	basic::Tracer my_tracer( "core.conformation", basic::t_warning );

	if ( verbose ) {
		int width = 14;
		my_tracer
			// << A( width, "torsion_id"   )
			<< A( width, "rsd_dihedral" )
			<< A( width, "atr_dihedral" )
			<< A( width, "rsd_torsion"  )
			<< A( width, "torsion"      )
			<< A( width, "atr_torsion"  )
			<< std::endl;
	}
	for ( Size i=1; i<= size(); ++i ) {
		Residue const & rsd( residue(i) );

		for ( Size r=1; r<= 2; ++r ) {
			id::TorsionType const type( r == 1 ? id::BB :
																					id::CHI );
			Size const n( r == 1 ? rsd.mainchain_atoms().size() : rsd.nchi() );


			for ( Size j=1; j<= n; ++j ) {
				TorsionID const tor_id( i, type, j );
				AtomID atom1,atom2,atom3,atom4;
				bool const fail
					( get_torsion_angle_atom_ids( tor_id, atom1, atom2, atom3, atom4 ) );
				if ( fail ) {
					if ( ( r == 1 ) &&
							( ( j ==   1 && fold_tree_->is_cutpoint( i-1 ) ) ||
							( j >= n-1 && fold_tree_->is_cutpoint( i ) ) ) ) continue;

					my_tracer << " missed torsion: " << tor_id << std::endl;
					continue;
				}

				Real const rsd_dihedral
					(numeric::dihedral(residue(atom1.rsd()).atom(atom1.atomno()).xyz(),
							residue(atom2.rsd()).atom(atom2.atomno()).xyz(),
							residue(atom3.rsd()).atom(atom3.atomno()).xyz(),
							residue(atom4.rsd()).atom(atom4.atomno()).xyz()));

				Real const atom_tree_dihedral
					( numeric::dihedral( atom_tree_->xyz( atom1 ),
						atom_tree_->xyz( atom2 ),
						atom_tree_->xyz( atom3 ),
						atom_tree_->xyz( atom4 ) ) );

				Real const rsd_torsion
					( r == 1 ? rsd.mainchain_torsion(j) : rsd.chi(j) );

				ASSERT_ONLY(Real const dev
					( std::abs( subtract_degree_angles(rsd_dihedral,atom_tree_dihedral))+
						std::abs( subtract_degree_angles(rsd_dihedral,rsd_torsion))+
						std::abs( subtract_degree_angles(rsd_dihedral,torsion(tor_id)))+
						std::abs( subtract_degree_angles(rsd_dihedral,
								atom_tree_torsion(tor_id))));)

				assert( dev < 1e-3 );

				if ( verbose ) {
					int width = 14;
					int precision = 4;
					my_tracer
						// << F( width, precision, tor_id                    )
						<< F( width, precision, rsd_dihedral              )
						<< F( width, precision, atom_tree_dihedral        )
						<< F( width, precision, rsd_torsion               )
						<< F( width, precision, torsion(tor_id)           )
						<< F( width, precision, atom_tree_torsion(tor_id) )
						<< std::endl;
				}
			} // j = torsionID
		} // r = torsiontype
	} // i = seqpos
}

// Show each residue in the conformation and its connections.
void
Conformation::show_residue_connections() const
{
	show_residue_connections(TR);
}



// Show each residue in the conformation and its connections.
/// @note This method is a rewrite of an earlier version to include an argument
///  for desired output stream.  This is to be more consistent with typical
///  show() methods and allows for simple << operator overloading for use in
///  PyRosetta. ~ Labonte
void
Conformation::show_residue_connections(std::ostream &os) const
{
	for (Size i = 1; i <= size(); ++i) {
		Residue const &res(*(residues_[i]));
		Size const nconn(res.n_residue_connections());
		os << "RESCON: " << i << ' ' << res.name() << " n-conn= " << nconn <<
			" n-poly= " << res.n_polymeric_residue_connections() <<
			" n-nonpoly= " << res.n_non_polymeric_residue_connections();
		for (Size j = 1; j <= nconn; ++j) {
			os << " conn# " << j << ' ' << res.residue_connect_atom_index( j )
			<< ' ' << res.connect_map(j).resid() << ' ' <<
			res.connect_map(j).connid();
		}
		os << std::endl;
	}
}


// Comparisons ///////////////////////////////////////////////////////////////////////////////////////////////////////

// determine the type of the ConformationOP
bool
Conformation::same_type_as_me( Conformation const & other, bool recurse /* = true */ ) const
{
	if ( recurse ) {
		return other.same_type_as_me( *this, false );
	} else {
		return true;
	}
}

/// @details Am I composed of the same sequence of ResidueType objects as other?
/// DANGER! Fails, unexpectedly, if a histidine pair has a different tautomerization.
/// tex, 9/12/08
/// The above comment points out that this method isn't a clearly defined idea, because there
/// are many ways that someone could imagine comparing the a sequence, including:
/// - cutpoint variants
/// - tautomers
/// - disulfides
/// For some applications, it's probably to get strings of the one-letter names of
/// the sequence and compare those instead.
bool
Conformation::sequence_matches( Conformation const & other ) const
{
	if ( size() != other.size() ) return false;
	for ( Size ii = 1, iiend = size(); ii <= iiend; ++ii ) {
		if ( & residue(ii).type() != & other.residue(ii).type() )	  return false;
		if ( residue(ii).chain() != other.residue(ii).chain() ) 	  return false;
		if ( !residue(ii).connections_match( other.residue( ii ) ) ) return false;
	}
	return true;
}

// General Properties ////////////////////////////////////////////////////////////////////////////////////////////////

///@note this is not a good test --Doug
bool
Conformation::is_fullatom() const {
	return is_residue_typeset( core::chemical::FA_STANDARD );
}

///@note this is not a good test --Doug
bool
Conformation::is_centroid() const {
	return is_residue_typeset( core::chemical::CENTROID );
}

///@note this is not a good test --Doug
bool
Conformation::is_residue_typeset( std::string tag ) const {
	// Empty poses aren't any residue typeset (special case prevents segfault)
	if( size() == 0 ) {
		TR.Warning << "WARNING: Attempted to determine the residue type set of an empty pose." << std::endl;
		return false;
	}

	Size const seqpos1( std::max( 1, (int) size() / 2 ) );
	Size const seqpos2( std::min( (int) size(), (int) size() / 2 + 2 ) );

	// if one of them is correct we are happy!
	bool correct = residue_type( seqpos1 ).residue_type_set().name() == tag;
	correct = correct || residue_type( seqpos2 ).residue_type_set().name() == tag;
	return correct;
}


// Chains ////////////////////////////////////////////////////////////////////////////////////////////////////////////

utility::vector1< Size > const &
Conformation::chain_endings() const
{
	return chain_endings_;
}

/// @remarks All positions must be strictly less than the number of
/// residues in the Conformation, otherwise the routine will fail fast.
/// Note that the last residue position is not counted as a chain end.
void
Conformation::chain_endings( utility::vector1< Size > const & endings )
{
	// make sure that all positions in the new endings list < the size of
	// the Conformation
	for ( utility::vector1< Size >::const_iterator i = endings.begin(), ie = endings.end(); i != ie; ++i ) {
		if ( (*i) >= size() ) {
			utility_exit_with_message("new chain endings list contains positions >= Conformation::size()");
		}
	}

	chain_endings_ = endings;
	std::sort( chain_endings_.begin(), chain_endings_.end() );
	rederive_chain_ids();
}

/// @remarks The last residue position is not counted as a chain ending.
/// Also increases the chain ID number by 1 for all residues upstream from seqpos.
void
Conformation::insert_chain_ending( Size const seqpos )
{
	assert( seqpos >= 1 && seqpos < size() ); // dont count last residue as chain ending

	chain_endings_.push_back( seqpos );
	std::sort( chain_endings_.begin(), chain_endings_.end() );

	rederive_chain_ids();
}

/// @remarks The last residue position is not counted as a chain ending.
void
Conformation::delete_chain_ending( Size const seqpos )
{
	utility::vector1< Size >::iterator it( std::find( chain_endings_.begin(), chain_endings_.end(), seqpos ) );
	if ( it == chain_endings_.end() ) {
		utility_exit_with_message("tried to delete nonexistent chain ending");
	}
	chain_endings_.erase( it );
	rederive_chain_ids();
}

void
Conformation::reset_chain_endings() {
	chain_endings_.clear();
	rederive_chain_ids();
}

void
Conformation::chains_from_termini()
{

	chain_endings_.clear();

	for ( Size i=1; i< size(); ++i ) {
		// should there be a chain ending between i and i+1 ?
		if ( residues_[i]->is_polymer() && residues_[i]->is_upper_terminus() ) {
			chain_endings_.push_back( i );
		}
	}

	rederive_chain_ids();
}

// Membranes /////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Given Embedding and Membrane Positions, Setup a Membrane Info Object
void
Conformation::setup_membrane( utility::vector1< std::pair< int, int > > embres_map, int membrane ) {

	// Create new membrane info object
	membrane_info_ = new MembraneInfo( this, embres_map, membrane );
}

/// @brief Returns a Membrane Info Object in the conformation
membrane::MembraneInfoOP
Conformation::membrane() {
	return membrane_info_;
}

// Trees /////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// @details setup atom tree as well from the fold tree
void
Conformation::fold_tree( FoldTree const & fold_tree_in )
{
	if ( size() != Size(fold_tree_in.nres()) ) {
		std::string msg;
		msg += "Conformation: fold_tree nres should match conformation nres. conformation nres: ";
		msg += string_of( size() );
		msg += " fold_tree nres: " + string_of( fold_tree_in.nres() ) ;
		utility_exit_with_message( msg );
	}
	update_residue_coordinates();
	(*fold_tree_) = fold_tree_in;
	setup_atom_tree();
}


// Residues //////////////////////////////////////////////////////////////////////////////////////////////////////////

ResidueCAPs
Conformation::const_residues() const
{
	// setup ResidueCAPs object
	conformation::ResidueCAPs const_rsds;
	for ( Size i=1; i<= size(); ++i ) {
		// pointer conversions are really annoying
		const_rsds.push_back( &( *residues_[i] ) );
	}
	return const_rsds;
}


/// @details add a residue into residues_ container, update its seqpos, chainid as well
/// fold tree and atoms.
/// Fires a LengthEvent::RESIDUE_APPEND signal.
void
Conformation::append_residue_by_jump(
	conformation::Residue const & new_rsd,
	Size const anchor_pos,
	std::string const& anchor_atom, // could be zero
	std::string const& root_atom, // ditto
	bool const start_new_chain // default false
)
{
	pre_nresidue_change();

	// handle first residue: special case
	if ( size() < 1 ) {
		append_residue( new_rsd, false, "", id::BOGUS_NAMED_ATOM_ID, start_new_chain );
		return;
	}

	assert( anchor_pos ); // should be set

	// now call our internal append method
	append_residue( new_rsd, true, root_atom, id::NamedAtomID( anchor_atom, anchor_pos ), start_new_chain );
}

/// @details insert a residue by jump
///  Fires a LengthEvent::RESIDUE_PREPEND signal.
void
Conformation::insert_residue_by_jump(
	Residue const & new_rsd_in,
	Size const seqpos, // desired seqpos of new_rsd
	Size anchor_pos, // in the current sequence numbering, ie before insertion of seqpos
	std::string const& anchor_atom, // could be ""
	std::string const& root_atom, // ditto
	bool new_chain
)
{
	pre_nresidue_change();
	ASSERT_ONLY(Size const old_size( size() );) //, new_size( old_size+1 );
	assert( old_size );
	runtime_assert( fold_tree_->is_cutpoint( seqpos-1 ) );

	// this handles all renumbering internal to the Residues, *_moved arrays
	residues_insert( seqpos, new_rsd_in, false, new_chain );

	Residue const & new_rsd( residue_( seqpos ) );

	assert( new_rsd.seqpos() == seqpos );

	fold_tree_->insert_residue_by_jump( seqpos, anchor_pos /* in the OLD numbering system */, anchor_atom, root_atom );

	insert_residue_into_atom_tree( new_rsd, *fold_tree_, const_residues(), *atom_tree_ );

	residue_torsions_need_updating_ = true;

	assert( atom_tree_->size() == size() && Size(fold_tree_->nres()) == size() );

	notify_length_obs( LengthEvent( this, LengthEvent::RESIDUE_PREPEND, seqpos, 1, &new_rsd ), false );
}

/// @details insert a residue by bond
///  Fires a LengthEvent::RESIDUE_PREPEND signal.
void
Conformation::insert_residue_by_bond(
                                     Residue const & new_rsd_in,
                                     Size const seqpos, // desired seqpos of new_rsd
                                     Size anchor_pos, // in the current sequence numbering, ie before insertion of seqpos
                                     bool const build_ideal_geometry, // = false,
                                     std::string const& anchor_atom, // could be ""
                                     std::string const& root_atom, // ditto
                                     bool new_chain,
                                     bool const lookup_bond_length // default false
                                     )
{
    pre_nresidue_change();
    Size const old_size( size() );
    assert( old_size );
    runtime_assert( fold_tree_->is_cutpoint( seqpos-1 ) );

    //Size const new_size( old_size+1 );
    
	// first renumber things
	utility::vector1< Size > old2new( old_size, 0 );
	for ( Size i=1; i<= old_size; ++i ) {
		if ( i< seqpos ) old2new[i] = i;
		else old2new[i] = i+1;
	}

    
    Residue const & anchor_rsd( residue_( anchor_pos ) ); // no call to residue(anchor_pos)

    ResidueOP ideal_geometry_rsd;
	if ( build_ideal_geometry ) {
		// get the geometries of the new bond
		chemical::ResidueConnection new_rsd_connection;
		chemical::ResidueConnection anchor_rsd_connection;
        
        core::Size anchor_atomid(anchor_rsd.type().atom_index(anchor_atom));
        core::Size anchor_connecting_id = anchor_rsd.type().residue_connection_id_for_atom(anchor_atomid);

        core::Size root_atomid(new_rsd_in.type().atom_index(root_atom));
        core::Size root_connecting_id = new_rsd_in.type().residue_connection_id_for_atom(root_atomid);

        // using a non-polymer residue-residue connection
        new_rsd_connection    = new_rsd_in.residue_connection(   root_connecting_id );
        anchor_rsd_connection = anchor_rsd.residue_connection( anchor_connecting_id );
        
		// this is a little wasteful, creating a new copy, but we need non-const access
		ideal_geometry_rsd = new_rsd_in.clone();
        
		if ( residue_coordinates_need_updating_ ) update_residue_coordinates( anchor_pos ); // safety
		orient_residue_for_ideal_bond( *ideal_geometry_rsd, new_rsd_connection, anchor_rsd, anchor_rsd_connection, *this, lookup_bond_length );

        // this handles all renumbering internal to the Residues, *_moved arrays
        residues_insert( seqpos, *ideal_geometry_rsd, false );

	    residues_[ seqpos ]->residue_connection_partner(        root_connecting_id, old2new[anchor_pos], anchor_connecting_id );
        residues_[ old2new[anchor_pos] ]->residue_connection_partner( anchor_connecting_id, seqpos,      root_connecting_id );
    }
    else {
        // this handles all renumbering internal to the Residues, *_moved arrays
        residues_insert( seqpos, new_rsd_in, false, new_chain );
    }

    Residue const & new_rsd( residue_( seqpos ) );
    assert( new_rsd.seqpos() == seqpos );


    fold_tree_->insert_residue_by_chemical_bond( seqpos, anchor_pos /* in the OLD numbering system */, anchor_atom, root_atom );
    
    insert_residue_into_atom_tree( new_rsd, *fold_tree_, const_residues(), *atom_tree_ );
    
    residue_torsions_need_updating_ = true;
    
    assert( atom_tree_->size() == size() && Size(fold_tree_->nres()) == size() );

    notify_length_obs( LengthEvent( this, LengthEvent::RESIDUE_PREPEND, seqpos, 1, &new_rsd ), false );
}

/// @details The default behavior is to append by a polymeric connection to the preceding residue
/// If we want to connect via a non-polymer connection, we give the connection number, anchor residue
/// and the connection number for the anchor residue. These connection numbers are wrt the connections_
/// arrays in Residue and ResidueType
/// If build_ideal_bond is TRUE it will transform the coordinates of the new residue so that the bond
/// geometry of the new bond is ideal according to the icoor_internal data in the residues.
/// Otherwise the incoming coordinates of new_rsd are preserved.
/// Fires a LengthEvent::RESIDUE_APPEND signal.
void
Conformation::append_residue_by_bond(
	conformation::Residue const & new_rsd,
	bool const build_ideal_geometry, // = false,
	int residue_connection_index, // = 0,
	Size anchor_pos, // = 0,
	int anchor_residue_connection_index, // = 0
	bool const start_new_chain, // default false
	bool const lookup_bond_length // default false
)
{
	if(!build_ideal_geometry) assert(!lookup_bond_length); // lookup only possible if building ideal geometry

	pre_nresidue_change();

	// handle first residue: special case
	if ( size() < 1 ) {
		append_residue( new_rsd, false, "", id::BOGUS_NAMED_ATOM_ID, start_new_chain );
		return;
	}

	Size const seqpos( size() + 1 );
	bool const polymer_connection( residue_connection_index == 0 );

	if ( polymer_connection ) anchor_pos = seqpos - 1; // polymer connection is to the preceding residue
	if ( polymer_connection ) anchor_residue_connection_index  = residue_( anchor_pos ).type().upper_connect_id();
	if ( polymer_connection ) residue_connection_index = new_rsd.type().lower_connect_id();

	Residue const & anchor_rsd( residue_( anchor_pos ) ); // no call to residue(anchor_pos)

	// debug
	if ( polymer_connection ) {
		// otherwise confirm that we are making a valid polymer bond
		//Don't use "is_upper_terminus" here as it fails to catch non-"terminus" residue types that have no upper connect id (e.g. methylamidated reisdues)
		if(!anchor_rsd.is_polymer() || !anchor_rsd.type().upper_connect_id()){
			std::stringstream err;
			err << "Can't create a polymer bond after residue " << anchor_pos
				<< " due to incompatible type: " << anchor_rsd.type().name();
			utility_exit_with_message(err.str());
		}
		else if(!new_rsd.is_polymer() || !new_rsd.type().lower_connect_id()){
			std::stringstream err;
			err << "Can't create a polymer bond to new residue " << seqpos
				<< " due to incompatible type: " << new_rsd.type().name();
			utility_exit_with_message(err.str());
		}
	} else {
		// if using a non-polymer connection, confirm that anchor_pos & anchor_residue_connection_index are set
		assert( anchor_pos && anchor_residue_connection_index );
	}

	//////////////////////////////////////////////////
	ResidueOP ideal_geometry_rsd;
	if ( build_ideal_geometry ) {
		// get the geometries of the new bond
		chemical::ResidueConnection new_rsd_connection;
		chemical::ResidueConnection anchor_rsd_connection;

		if ( polymer_connection ) {
			// polymer connection
			new_rsd_connection    =    new_rsd.lower_connect();
			anchor_rsd_connection = anchor_rsd.upper_connect();
		} else {
			// using a non-polymer residue-residue connection
			new_rsd_connection    =    new_rsd.residue_connection(        residue_connection_index );
			anchor_rsd_connection = anchor_rsd.residue_connection( anchor_residue_connection_index );
		}

		// this is a little wasteful, creating a new copy, but we need non-const access
		ideal_geometry_rsd = new_rsd.clone();

		if ( residue_coordinates_need_updating_ ) update_residue_coordinates( anchor_pos ); // safety
		orient_residue_for_ideal_bond( *ideal_geometry_rsd, new_rsd_connection, anchor_rsd, anchor_rsd_connection, *this, lookup_bond_length );

	}

	//////////////////////////////////////////////////////
	// determine anchor and root atoms
	id::AtomID anchor_id(0,anchor_pos);
	int root_atomno;

	if ( polymer_connection ) {
		// connecting by a polymer bond
		root_atomno        =    new_rsd.lower_connect_atom();
		anchor_id.atomno() = anchor_rsd.upper_connect_atom();
	} else {
		// non-polymer connection
		root_atomno        =    new_rsd.residue_connection(        residue_connection_index ).atomno();
		anchor_id.atomno() = anchor_rsd.residue_connection( anchor_residue_connection_index ).atomno();
	}

	std::string root_atom = "";
	if ( root_atomno ) root_atom = new_rsd.atom_name( root_atomno );

	// Now call our internal append method
	if ( build_ideal_geometry ) {
        if (polymer_connection) {
            append_residue( *ideal_geometry_rsd,
                           false,
                           "",
                           id::BOGUS_NAMED_ATOM_ID,
                           start_new_chain );
        }
        else {
            append_residue(
                *ideal_geometry_rsd,
                false /* not by a jump */,
                root_atom,
                atom_id_to_named_atom_id(anchor_id, anchor_rsd),
                start_new_chain
            );
        }
	}	else {
		append_residue(
			new_rsd,
			false /* not by a jump */,
			root_atom,
			atom_id_to_named_atom_id(anchor_id, anchor_rsd),
			start_new_chain );
	}

	if ( build_ideal_geometry && polymer_connection ) rebuild_polymer_bond_dependent_atoms( seqpos - 1 );

	// update the connection info stored inside the residues
	//if ( !polymer_connection ) {
		residues_[ seqpos     ]->residue_connection_partner(        residue_connection_index, anchor_pos, anchor_residue_connection_index );
		residues_[ anchor_pos ]->residue_connection_partner( anchor_residue_connection_index, seqpos,            residue_connection_index );
	//}

} // append_residue_by_bond

/// @details Fires a LengthEvent::RESIDUE_PREPEND signal at seqpos + 1 due to
///  insert_polymer_residue() call.
void
Conformation::append_polymer_residue_after_seqpos(
	Residue const & new_rsd,
	Size const seqpos,
	bool const build_ideal_geometry // = false
)
{
	pre_nresidue_change();
	assert( !new_rsd.is_lower_terminus() );

	ResidueOP ideal_geometry_rsd;
	if ( build_ideal_geometry ) {
		if ( residue_coordinates_need_updating_ ) update_residue_coordinates( seqpos );
		Residue const & anchor_rsd( residue_( seqpos ) ); // not residue(seqpos)
		assert( !anchor_rsd.is_upper_terminus() );
		// this is a little wasteful, creating a new copy, but we need non-const access
		ideal_geometry_rsd = new_rsd.clone();

		orient_residue_for_ideal_bond(*ideal_geometry_rsd, new_rsd.lower_connect(), anchor_rsd, anchor_rsd.upper_connect(), *this );
	}
	bool const join_lower( true );
	bool const join_upper( !fold_tree_->is_cutpoint( seqpos ) );

	if ( build_ideal_geometry ) insert_polymer_residue( *ideal_geometry_rsd, seqpos+1, join_lower, join_upper );
	else                        insert_polymer_residue(             new_rsd, seqpos+1, join_lower, join_upper );

	if ( build_ideal_geometry ) {
		rebuild_polymer_bond_dependent_atoms( seqpos );
	}

	assert( dof_moved_[seqpos+1].size() == new_rsd.natoms() &&
					xyz_moved_[seqpos+1].size() == new_rsd.natoms() );
}

void
Conformation::safely_append_polymer_residue_after_seqpos(
	Residue const & new_rsd,
	Size const seqpos,
	bool const build_ideal_geometry // = false
)
{
	pre_nresidue_change();
	core::conformation::remove_upper_terminus_type_from_conformation_residue( *this, seqpos );
	append_polymer_residue_after_seqpos( new_rsd, seqpos, build_ideal_geometry );
}

/// @details Fires a LengthEvent::RESIDUE_PREPEND signal.
void
Conformation::prepend_polymer_residue_before_seqpos(
	Residue const & new_rsd,
	Size const seqpos,
	bool const build_ideal_geometry // = false
)
{
	pre_nresidue_change();
	assert( !new_rsd.is_upper_terminus() );

	ResidueOP ideal_geometry_rsd;
	if ( build_ideal_geometry ) {
		if ( residue_coordinates_need_updating_ ) update_residue_coordinates( seqpos );
		Residue const & anchor_rsd( residue_( seqpos ) ); // not residue(seqpos)
		assert( !anchor_rsd.is_lower_terminus() );
		// this is a little wasteful, creating a new copy, but we need non-const access
		ideal_geometry_rsd = new_rsd.clone();

		orient_residue_for_ideal_bond(*ideal_geometry_rsd, new_rsd.upper_connect(), anchor_rsd, anchor_rsd.lower_connect(), *this );
	}

	bool const join_upper( true );
	bool const join_lower( !fold_tree_->is_cutpoint( seqpos-1 ) );

	if ( build_ideal_geometry ) insert_polymer_residue( *ideal_geometry_rsd, seqpos, join_lower, join_upper );
	else                        insert_polymer_residue(             new_rsd, seqpos, join_lower, join_upper );

	if ( build_ideal_geometry ) {
		rebuild_polymer_bond_dependent_atoms( seqpos );
	}

	assert( dof_moved_[seqpos].size() == new_rsd.natoms() &&
					xyz_moved_[seqpos].size() == new_rsd.natoms() );
}

/// @details Fires a LengthEvent::RESIDUE_PREPEND signal.
void
Conformation::safely_prepend_polymer_residue_before_seqpos(
	Residue const & new_rsd,
	Size const seqpos,
	bool const build_ideal_geometry // = false
)
{
	core::conformation::remove_lower_terminus_type_from_conformation_residue( *this, seqpos );
	prepend_polymer_residue_before_seqpos( new_rsd, seqpos, build_ideal_geometry );
}

/// @details replace a residue
///  Fires an IdentityEvent signal.
void
Conformation::replace_residue(
	Size const seqpos,
	Residue const & new_rsd_in,
	bool const orient_backbone
)
{
	// temporarily hold onto the current residue
	ResidueOP old_rsd_ptr( residues_[ seqpos ] );
	Residue const & old_rsd( *old_rsd_ptr );

	// helper function, hides direct access to residues_
	residues_replace( seqpos, new_rsd_in );

	// reference to new residue
	Residue & new_rsd( *residues_[ seqpos ] );

	// transform coordinates to align backbone with current backbone
	// note that this also copies backbone coords from old_rsd
	// ie, it's the same as the rotamer-building residue placement operation
	//
	if ( orient_backbone ) {
		residues_[seqpos]->place( old_rsd, *this );
	}

	replace_residue_in_atom_tree( new_rsd, *fold_tree_, const_residues(), *atom_tree_ );

	residue_torsions_need_updating_ = true;
// 	if ( !residue_torsions_need_updating_ ) update_residue_torsions( seqpos );

	notify_identity_obs( IdentityEvent( this, IdentityEvent::RESIDUE, seqpos, &new_rsd ), false );
}

// function to replace a residue based on superposition on the specified input atom pairs
/// @note  NOTE: at the moment, only superposition on 3 atoms works
/// Fires an IdentityEvent signal.
void
Conformation::replace_residue(
	Size const seqpos,
	Residue const & new_rsd_in,
	utility::vector1< std::pair< std::string, std::string > > const & atom_pairs
)
{
	//first we'll set the correct atom coordinates, then we'll let
	//the other replace residue function deal with everything else

	// temporarily hold onto the current residue
	ResidueOP old_rsd_ptr( residues_[ seqpos ] );
	Residue const & old_rsd( *old_rsd_ptr );

	Residue new_rsd = new_rsd_in;

	new_rsd.orient_onto_residue( old_rsd, atom_pairs );

	//2. call replace res without orienting backbone( no coordinates will be changed ),
	//   only the conformation internal data structures are updated

	replace_residue( seqpos, new_rsd, false);

}

/// @details delete a polymer residues
///  Fires a LengthEvent::RESIDUE_DELETE signal.
void
Conformation::delete_polymer_residue( Size const seqpos )
{
	pre_nresidue_change();
	assert( !fold_tree_->is_jump_point( seqpos ) );

	residues_delete( seqpos ); // handles renumbering of residues, _moved, chains

	// delete from AtomTree, foldtree
	// this atom_tree_ call could be more robust if we passed in the fold_tree_
	// currently it assumes 1 incoming connxn,  <=1 outgoing connxn
	atom_tree_->delete_seqpos( seqpos );
	fold_tree_->delete_seqpos( seqpos );

	residue_torsions_need_updating_ = true; // could reupdate before and after

	assert( atom_tree_->size() == size() && Size(fold_tree_->nres()) == size() );

	notify_length_obs( LengthEvent( this, LengthEvent::RESIDUE_DELETE, seqpos, -1, NULL ), false );
}

/// @details  Delete a residue from the Conformation the slow way -- triggers a rebuild of the atomtree
///  Fires a LengthEvent::RESIDUE_DELETE signal.
/// @note  Could be upstream and/or downstream of a jump or chemical edge, or the root of the tree
/// @note  Not well-tested. Expect funny behavior in new or different situations (email pbradley@fhcrc.org)
///
/// LOGIC: uses fold_tree.delete_seqpos to handle shifting the topology around if necessary, then calls setup_atom_tree
void
Conformation::delete_residue_slow( Size const seqpos )
{
	pre_nresidue_change();
	fold_tree_->delete_seqpos( seqpos );

	residues_delete( seqpos );

	// setup_atom_tree() can fire an XYZEvent - make sure the LengthEvent fires first
	notify_length_obs( LengthEvent( this, LengthEvent::RESIDUE_DELETE, seqpos, -1, NULL ), false );

	setup_atom_tree();

	residue_torsions_need_updating_ = true;

	assert( atom_tree_->size() == size() && Size(fold_tree_->nres()) == size() );

}

/// @details  Like above but only one call to setup_atom_tree
///  Fires a LengthEvent::RESIDUE_DELETE signal.
void
Conformation::delete_residue_range_slow( Size const range_begin, Size const range_end )
{
	pre_nresidue_change();
	Size const range_size( range_end - range_begin + 1 );
	assert( range_size >= 1);
	for ( Size i=1; i<= range_size; ++i ) {
		fold_tree_->delete_seqpos( range_begin );
		residues_delete( range_begin );
	}

	// setup_atom_tree() can fire an XYZEvent - make sure the LengthEvent fires first
	notify_length_obs( LengthEvent( this, LengthEvent::RESIDUE_DELETE, range_begin, -range_size,  NULL ), false );

	setup_atom_tree();

	residue_torsions_need_updating_ = true;

	assert( atom_tree_->size() == size() && Size(fold_tree_->nres()) == size() );

}


/// @details Returns a mask of residues over which scoring is restricted.
/// Only these residues will be used in constructing the neighbor list.
utility::vector1<bool>
Conformation::get_residue_mask() const {
	Size const nres( residues_.size() );
	return 	utility::vector1<bool>(nres, true);
}

/// @details Returns a weight to be used when scoring this residue.
Real
Conformation::get_residue_weight(core::Size , core::Size) const {
	return 1.0;
}


// Bonds, Connections, Atoms, & Stubs ////////////////////////////////////////////////////////////////////////////////

void
Conformation::detect_bonds()
{
	using namespace graph;

	// Determine how many residues have incomplete connections.
	utility::vector1< Size > resid_2_incomp( size(), 0 );  // residue id to "incomplete connection number"
	Size num_incomp( 0 );
	for ( Size ii = 1; ii <= size(); ++ii ) {
		if ( residue(ii).has_incomplete_connection() ) {
			++num_incomp;
			resid_2_incomp[ ii ] = num_incomp;
		}
	}
	if ( num_incomp == 0 ) return;

	TR.Debug << "Looking for connection partners for " << num_incomp << " residues" << std::endl;

	utility::vector1< Size > incomp_2_resid( num_incomp );  // "incomplete connection number" to residue id
	for ( Size ii = 1; ii <= size(); ++ii ) {
		if ( resid_2_incomp[ ii ] != 0 ) {
			incomp_2_resid[ resid_2_incomp[ ii ] ] = ii;
		}
	}

	// Create point graph of nbr_atoms of incomplete residues only.
	// Also, calculate maximum distance that a connected residue could be.
	PointGraphOP pg = new PointGraph;
	pg->set_num_vertices( num_incomp );
	Distance maxrad( 0.0 );
	Distance maxd( 0.0 );
	for ( Size ii = 1; ii <= num_incomp; ++ii ) {
		Residue const & ii_res = residue( incomp_2_resid[ ii ] );
		pg->get_vertex(ii).data().xyz() = ii_res.atoms()[ ii_res.nbr_atom() ].xyz();
		if ( ii_res.nbr_radius() > maxrad ) maxrad = ii_res.nbr_radius();
		for ( Size jj = 1; jj <= ii_res.type().n_residue_connections(); ++jj ) {
			if ( ii_res.connection_incomplete( jj ) ) {
				if ( maxd < ii_res.type().residue_connection(jj).icoor().d() ) {
					maxd = ii_res.type().residue_connection(jj).icoor().d();
				}
			}
		}
	}

	// two Angstrom extra radius for finding bonds... very generous
	maxd += 2.0;
	Distance neighbor_cutoff = maxrad + maxd;
	find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, neighbor_cutoff );

	// Iterate across neighbors of incomplete residues; compare incomplete connection points against each other.
	for ( Size ii = 1; ii <= num_incomp; ++ii ) {
		Size const ii_resid = incomp_2_resid[ ii ];
		Size const ii_n_conn = residue( ii_resid ).type().n_residue_connections();
		Residue const & ii_res( residue( ii_resid ) );
		for ( Size jj = 1; jj <= ii_n_conn; ++jj ) {
			if ( ! ii_res.connection_incomplete( jj ) ) continue;

			// Get the atom index of the jjth ResidueConnection of the iith incomplete residue.
			Size const jjatom = ii_res.residue_connection( jj ).atomno();
			bool multiple_connections_for_jjatom = ii_res.type().n_residue_connections_for_atom( jjatom ) > 1;

			Distance best_match( 0.0 ), best_jj( 0.0 ), best_kk( 0.0 );
			Size best_match_resid( 0 );
			Size best_match_connid( 0 );

			// Search the PointGraph for neighbor residues to the iith incomplete residue.
			for ( PointGraph::UpperEdgeListConstIter ii_iter = pg->get_vertex( ii ).upper_edge_list_begin(),
					ii_end_iter = pg->get_vertex( ii ).upper_edge_list_end();
					ii_iter != ii_end_iter; ++ii_iter ) {
				Size const neighb_id = ii_iter->upper_vertex();
				Size const neighb_resid = incomp_2_resid[ neighb_id ];
				Residue const & neighb( residue( neighb_resid ) );
				Size const neighb_n_conn = neighb.type().n_residue_connections();

				// Get the atom index of the kkth ResidueConnection of the iith incomplete residue's neighbor.
				for ( Size kk = 1; kk <= neighb_n_conn; ++kk ) {
					if ( ! neighb.connection_incomplete( kk ) ) continue;

					Size const kkatom = neighb.residue_connection( kk ).atomno();
					bool multiple_connections_for_kkatom = neighb.type().n_residue_connections_for_atom(kkatom) > 1;

					// Calculate distances between expected atoms and actual atoms for these two connections.
					Distance kk_distance = ii_res.connection_distance( *this,
						jj, neighb.atom( kkatom /*neighb.type().residue_connection( kk ).atomno()*/ ).xyz() );
					Distance jj_distance = neighb.connection_distance( *this,
						kk, ii_res.atom( jjatom /*ii_res.type().residue_connection( jj ).atomno()*/ ).xyz() );

					//std::cout << "kk_distance: " << kk_distance << " jj_distance: " << jj_distance << std::endl;

					if ( multiple_connections_for_jjatom ) {
						Vector jj_expected_coord = ii_res.residue_connection( jj ).icoor().build( ii_res, *this );
						jj_distance += jj_expected_coord.distance( neighb.xyz( kkatom ) );
					}
					if ( multiple_connections_for_kkatom ) {
						Vector kk_expected_coord = neighb.residue_connection( kk ).icoor().build( neighb, *this );
						kk_distance += kk_expected_coord.distance( ii_res.xyz( jjatom ));
					}

					if ( best_match_resid == 0 || best_match > kk_distance + jj_distance ) {
						best_match = kk_distance + jj_distance;
						best_jj = jj_distance; best_kk = kk_distance;
						best_match_resid = neighb_resid;
						best_match_connid = kk;
					}
				}
			}  // Loop while the iith incomplete residue still has neighbors.

			if ( best_match_resid == 0 ) {
				TR.Warning << "Failed to find a residue connection for residue " << ii_resid <<
						" with connection point " << jj << std::endl;
				continue;
			}

			TR.Info << "Connecting residues: " << ii_resid << " ( " << ii_res.name();
			TR.Info << " ) and " << best_match_resid << " ( " << residue( best_match_resid ).name() << " )";
			TR.Info << " at atoms ";
			TR.Info << ii_res.type().atom_name( ii_res.type().residue_connection( jj ).atomno() );
			TR.Info << " and ";
			TR.Info << residue( best_match_resid ).type().atom_name(
					residue( best_match_resid ).type().residue_connection( best_match_connid ).atomno() );
			TR.Info << std::endl;
			TR.Info << " with mutual distances: " << best_jj << " and " << best_kk << std::endl;

			residues_[ ii_resid ]->residue_connection_partner( jj, best_match_resid, best_match_connid );
			residues_[ best_match_resid ]->residue_connection_partner( best_match_connid, ii_resid, jj );
		}  // Loop while the iith incomplete residue still has connections.
	}  // Loop to the next incomplete residue.
}

void
Conformation::detect_pseudobonds()
{
	for ( Size ii = 1; ii <= size(); ++ii ) {
		Residue const & ii_res( residue_( ii ) );
		if ( ii_res.is_coarse() ) continue;
		Size const ii_nresconn = ii_res.type().n_residue_connections();
		for ( Size jj = 1; jj <= ii_nresconn; ++jj ) {
			Size const jj_atid = ii_res.residue_connection( jj ).atomno();
			for ( Size kk = jj + 1; kk <= ii_nresconn; ++kk ) {
				Size const kk_atid = ii_res.residue_connection( kk ).atomno();
				if ( ii_res.path_distance( kk_atid, jj_atid ) < 2 ) {
					// found a through-1 pseudobond
					// add pseudobonds between the residues that connect
					// to connection points kk and jj
					Size const jj_conn_res( ii_res.residue_connection_partner( jj ));
					Size const jj_conn_id(  ii_res.residue_connection_conn_id( jj ));
					Size const kk_conn_res( ii_res.residue_connection_partner( kk ));
					Size const kk_conn_id(  ii_res.residue_connection_conn_id( kk ));

					Size const lr  = jj_conn_res < kk_conn_res ? jj_conn_res : kk_conn_res;
					Size const lri = jj_conn_res < kk_conn_res ? jj_conn_id  : kk_conn_id ;
					Size const ur  = jj_conn_res < kk_conn_res ? kk_conn_res : jj_conn_res;
					Size const uri = jj_conn_res < kk_conn_res ? kk_conn_id  : jj_conn_id ;
					if(! (lr && lri && ur && uri)) continue;
					add_pseudobond( lr, lri, ur, uri, ii_res.path_distance( kk_atid, jj_atid ));

					TR.Info << "Adding PseudoBond between residues " << lr << " " << ur << ", connecting atoms ";
					TR.Info << residue_( lr ).type().atom_name( residue_( lr ).type().residue_connection( lri ).atomno() ) << " and ";
					TR.Info << residue_( ur ).type().atom_name( residue_( ur ).type().residue_connection( uri ).atomno() );
					TR.Info << " through residue " << ii << " at atoms " << ii_res.type().atom_name( jj_atid );
					TR.Info << " & " << ii_res.type().atom_name( kk_atid );
					TR.Info << " which are separated by " << ii_res.path_distance(kk_atid, jj_atid );
					TR.Info << " bond"<< (ii_res.path_distance(kk_atid, jj_atid ) == 1 ? "" : "s") <<  std::endl;
				}

				if ( kk_atid == jj_atid ) {
					// two connections on the same atom. Look for through-2 pseudobonds.
					utility::vector1< chemical::ResConnID > two_neighbors(2);
					two_neighbors[1] = ii_res.actual_residue_connection( jj );
					two_neighbors[2] = ii_res.actual_residue_connection( kk );
					for ( Size ll = 1; ll <= 2; ++ll ) {
						Size const ll_resid = two_neighbors[ ll ].resid();
						if ( ll_resid < ii ) continue; // only add through-2 pseudobonds through upper residues to avoid duplication

						Size const ll_connid = two_neighbors[ ll ].connid();
						Residue const & ll_res = residue( ll_resid );
						Size const ll_atid = ll_res.residue_connection( ll_connid ).atomno();
						Size const ll_nconn = ll_res.type().n_residue_connections();
						for ( Size mm = 1; mm <= ll_nconn; ++mm ) {
							if ( mm == ll_connid ) continue;
							Size const mm_atid = ll_res.type().residue_connection( mm ).atomno();
							if ( mm_atid != ll_atid ) continue;

							// arriving here, we've now identified that kk_atid == jj_atid and that
							// mm_atid == ll_atid; in other words, we have two atoms that have two
							// connections each, and these two atoms are bound to each other;  this pair of atoms
							// forms a through-2 connection.

							// add a psuedobond between the residue two_neighbors[ ll == 1 ? : 2, 1 ].resid()
							// and ll_res.residue_conections( mm ).resid()

							Size const ll_conn_res( ii_res.residue_connection_partner( ll == 1 ? kk : jj ));
							Size const ll_conn_id(  ii_res.residue_connection_conn_id( ll == 1 ? kk : jj ));
							Size const mm_conn_res( ll_res.residue_connection_partner( mm ));
							Size const mm_conn_id(  ll_res.residue_connection_conn_id( mm ));

							Size const lr  = ll_conn_res < mm_conn_res ? ll_conn_res : mm_conn_res;
							Size const lri = ll_conn_res < mm_conn_res ? ll_conn_id  : mm_conn_id ;
							Size const ur  = ll_conn_res < mm_conn_res ? mm_conn_res : ll_conn_res;
							Size const uri = ll_conn_res < mm_conn_res ? mm_conn_id  : ll_conn_id ;

							add_pseudobond( lr, lri, ur, uri, 3 /*nbonds*/ );

							TR.Info << "Adding PseudoBond between residues " << lr << " " << ur << ", connecting atoms ";
							TR.Info << residue_( lr ).type().atom_name( residue_( lr ).type().residue_connection( lri ).atomno() ) << " and ";
							TR.Info << residue_( ur ).type().atom_name( residue_( ur ).type().residue_connection( uri ).atomno() );
							TR.Info << " through two residues: " << ii << " and " << ll_resid << " at atoms " << ii_res.type().atom_name( jj_atid );
							TR.Info << " & " << ll_res.type().atom_name( mm_atid );
							TR.Info << " which are separated by 1 bond" <<  std::endl;
						}
					}
				}
			}
		}
	}
}


// Declare that a chemical bond exists between two residues
void
Conformation::declare_chemical_bond(
	Size const seqpos1,
	std::string const & atom_name1,
	Size const seqpos2,
	std::string const & atom_name2
)
{

	Residue & rsd1( *residues_[ seqpos1 ] );
	Residue & rsd2( *residues_[ seqpos2 ] );

	// find the connection ids
	Size const atom1( rsd1.atom_index( atom_name1 ) );
	Size const atom2( rsd2.atom_index( atom_name2 ) );
	Size connid1(0);
	for ( Size connid=1; connid<= rsd1.n_residue_connections(); ++connid ) {
		if ( Size(rsd1.residue_connection( connid ).atomno()) == atom1 ) {
			connid1 = connid;
		}
	}
	Size connid2(0);
	for ( Size connid=1; connid<= rsd2.n_residue_connections(); ++connid ) {
		if ( Size(rsd2.residue_connection( connid ).atomno()) == atom2 ) {
			connid2 = connid;
		}
	}
	if ( !connid1 ) utility_exit_with_message( rsd1.name()+" doesnt have connection at "+atom_name1 );
	if ( !connid2 ) utility_exit_with_message( rsd2.name()+" doesnt have connection at "+atom_name2 );
	rsd1.residue_connection_partner( connid1, seqpos2, connid2 );
	rsd2.residue_connection_partner( connid2, seqpos1, connid1 );

	//A new chemical bond might result in new torsion angles that need to be updated.  Update now:
	residue_torsions_need_updating_ = true;
	update_residue_torsions();
}

/// @brief Rebuilds the atoms that are depenent on polymer bonds for the specified residue only.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
Conformation::rebuild_polymer_bond_dependent_atoms_this_residue_only ( Size const seqpos )
{
	rebuild_polymer_bond_dependent_atoms( seqpos, 1 );
	rebuild_polymer_bond_dependent_atoms( seqpos, -1 );
	return;
}

/// @details rebuilds the atoms that are dependent on the bond between seqpos and seqpos+1 for their torsion offset
void
Conformation::rebuild_polymer_bond_dependent_atoms( Size const seqpos )
{
	rebuild_polymer_bond_dependent_atoms( seqpos  ,  1 );
	rebuild_polymer_bond_dependent_atoms( seqpos+1, -1 );
}

/// @note  This could be rewritten to avoid a refold using set_bond_angle, etc might be safer since
/// by-xyz ignores possibility of propagating dependencies from the moving atoms...
void
Conformation::rebuild_residue_connection_dependent_atoms( Size const seqpos, Size const connid )
{
	// assumes no interdependencies among the atoms being built
	update_residue_coordinates();
	Residue const & rsd( residue_(seqpos) );
	for ( Size i=1, ie=rsd.natoms(); i<= ie; ++i ) {
		if ( rsd.icoor(i).depends_on_residue_connection( connid ) ) {
			set_xyz( AtomID(i,seqpos), rsd.icoor(i).build( rsd, *this ) );
		}
	}
}


id::AtomID
Conformation::inter_residue_connection_partner(
	Size const seqpos,
	int const connection_index
) const
{
	Residue const & rsd( *residues_[ seqpos ] );
	return id::AtomID( rsd.residue_connection        ( connection_index ).atomno(),
										 rsd.residue_connection_partner( connection_index ) );
}

/// @details virtual atoms are excluded by default
utility::vector1<id::AtomID>
Conformation::bonded_neighbor_all_res(
	id::AtomID atomid,
	bool virt // = false
) const
{
	conformation::Residue const & primary_residue(residue(atomid.rsd()));

	utility::vector1<id::AtomID> neighbors;

	// add all the atoms in the same residue
	chemical::AtomIndices const & intrares_atomnos(primary_residue.bonded_neighbor(atomid.atomno()));
	for (Size i = 1; i <= intrares_atomnos.size(); ++i) {
		if (virt || ! primary_residue.is_virtual(intrares_atomnos[i]) ) {
			neighbors.push_back(id::AtomID(intrares_atomnos[i], atomid.rsd()));
		}
	}

	// add all the atoms in other residues
	Size const num_connections(primary_residue.n_residue_connections());
	for (Size i = 1; i <= num_connections; ++i) {
		if (primary_residue.residue_connect_atom_index(i) == atomid.atomno() &&
				!primary_residue.connection_incomplete(i)) {
			chemical::ResConnID resconid(primary_residue.actual_residue_connection(i));
			Size connected_atomno(residue(resconid.resid()).residue_connect_atom_index(resconid.connid()));
			if (virt || ! residue(resconid.resid()).is_virtual(connected_atomno)) {
				neighbors.push_back(id::AtomID(connected_atomno, resconid.resid()));
			}
		}
	}
	return neighbors;
}

void
Conformation::fill_missing_atoms(
	id::AtomID_Mask missing // make local copy
)
{
	using namespace basic;

	update_residue_coordinates(); //since we will be using residue_(seqpos) below

	for ( Size i=1; i<= size(); ++i ) { //Loop through all residues in this conformation
		Residue const & rsd( residue_(i) ); // prevent many many calls to update_residue_torsions()
		Size const natoms( rsd.natoms() );
		Size tries(0);
		while ( true ) {
			++tries;
			if ( tries > 10000 ) {
				utility_exit_with_message("too many tries in fill_missing_atoms!");
			}
			bool any_missing( false );
			Size num_present = 0;
			for ( Size j=1; j<= natoms; ++j ) {
				if ( !missing[ AtomID( j, i ) ] ) num_present += 1;
			}
			for ( Size j=1; j<= natoms; ++j ) { //Loop through all atoms in this residue
				AtomID const id( j, i );
				if ( missing[ id ] ) {
					// Virtual atoms don't get written to the PDB file, so we shouldn't expect them in input.
					if ( rsd.atom_type(j).is_heavyatom() && ! rsd.is_virtual(j)) {
						TR.Warning << "[ WARNING ] missing heavyatom: " << rsd.atom_name(j) <<
							" on residue " << rsd.name() << ' ' << i << std::endl;
					}
					any_missing = true;

					if(rsd.icoor(j).depends_on_residue_connection()) {
						//if the coordinates of this atom depend on a residue connection
						//NOTE: we assume that the coordinates of the root atom are never connection-dependent (which is true).
						bool needed_atom_is_missing(false);
						for(core::Size ii=1; ii<=3; ++ii) { //Check whether needed residues are missing, ignoring connections.
							if(!rsd.icoor(j).stub_atom(ii).is_connect() && missing[rsd.icoor(j).stub_atom(ii).atom_id(i, *this)] ) {
								needed_atom_is_missing = true;
								break;
							}
						}
						if(!needed_atom_is_missing) { //If no needed atoms are missing
								set_xyz( id , rsd.build_atom_ideal( j, *this ) );
								missing[ id ] = false;
								num_present += 1;
						}
					} else { //if the coordinates of this atom do not depend on a residue connection
						// check if our stub atoms are all present
						AtomID const
							stub_atom1( rsd.icoor( j ).stub_atom1().atom_id( i, *this ) ),
							stub_atom2( rsd.icoor( j ).stub_atom2().atom_id( i, *this ) ),
							stub_atom3( rsd.icoor( j ).stub_atom3().atom_id( i, *this ) );

						if ( num_present < 3 && !missing[ stub_atom1 ] ) {
							// with < 3 atoms, we can't build any stubs, so we're stuck forever unless we punt:
							//if ( rsd.natoms() == 3 && !missing[ stub_atom1 ] ) {
							// special case, eg HOH water, O is present H's missing
							using numeric::random::uniform;
							Vector xyz1( xyz( stub_atom1 ) ), xyz2( xyz( stub_atom2 ) ), xyz3( xyz( stub_atom3 ) );
							if ( missing[ stub_atom2 ] ) xyz2 = Vector( uniform(), uniform(), uniform() );
							if ( missing[ stub_atom3 ] ) xyz3 = Vector( uniform(), uniform(), uniform() );
							kinematics::Stub const stub( xyz1, xyz2, xyz3 );
							set_xyz( id, stub.spherical( rsd.icoor(j).phi(), rsd.icoor(j).theta(), rsd.icoor(j).d() ) );
							missing[id] = false;
							num_present += 1;

						} else if ( id == stub_atom1 || id == stub_atom2 || id == stub_atom3 ) {
							// the root atom of the default residue tree or one of it's stub_atoms...
							// special case requires careful handling:
							// build new residue trees for the ideal-coordinates residue, look for one
							// that doesn't have this problem...

							ResidueOP tmp_rsd( new Residue( rsd.type(), false /*dummy arg*/ ) );
							tmp_rsd->seqpos( i );

							for ( Size root_atomno=1; root_atomno<= natoms; ++root_atomno ) {
								if ( root_atomno == j ) continue;

								kinematics::AtomPointer2D atom_pointer( i );
								build_residue_tree( root_atomno, *tmp_rsd, atom_pointer[i], true/*root is jump*/ );
								AtomTree rsd_tree( atom_pointer );

								AtomID const
									new_stub_atom1( rsd_tree.atom( id ).input_stub_atom1_id() ),
									new_stub_atom2( rsd_tree.atom( id ).input_stub_atom2_id() ),
									new_stub_atom3( rsd_tree.atom( id ).input_stub_atom3_id() );

								if ( !missing[ new_stub_atom1 ] && !missing[ new_stub_atom2 ] && !missing[ new_stub_atom3 ] ) {
									TR.Warning << "[ WARNING ] Building missing atom (" << rsd.atom_name( j ) <<
										") at root of residue tree, using stubs: " <<
										rsd.atom_name( new_stub_atom1.atomno() ) << ' ' <<
										rsd.atom_name( new_stub_atom2.atomno() ) << ' ' <<
										rsd.atom_name( new_stub_atom3.atomno() ) <<
										"\nThis probably means that a torsion angle is being taken from the ideal residue and"
										"\nshould be further optimized..." << std::endl;
									kinematics::Stub stub
										( rsd.xyz( new_stub_atom1.atomno() ),
											rsd.xyz( new_stub_atom2.atomno() ),
											rsd.xyz( new_stub_atom3.atomno() ) );
									set_xyz( id, stub.spherical( rsd_tree.dof( DOF_ID( id, id::PHI ) ),
										rsd_tree.dof( DOF_ID( id, id::THETA ) ),
										rsd_tree.dof( DOF_ID( id, id::D ) ) ) );
									missing[ id ] = false;
									num_present += 1;
									break;
								}
							} // root_atomno = 1,natoms
						} else {
							// typical case
							if ( !missing[ stub_atom1 ] && !missing[ stub_atom2 ] && !missing[ stub_atom3 ] ) {
								set_xyz( id , rsd.build_atom_ideal( j, *this ) );
								missing[ id ] = false;
								num_present += 1;
							}
						}
					}
				} //if(rsd.icoor(j).depends_on_residue_connection())
			}
			if ( !any_missing ) break;
		}
	} // i = 1..size()
}

/// @details There is a version of this in Residue.hh.
/// The problem with that function is that accessing a residue in the pose triggers a refold
/// which is really slow. We do not need a correctly folded residue to makethis check
bool
Conformation::atom_is_backbone_norefold( Size const pos, Size const atomno ) const
{
	return ( ( *residues_[ pos ] ).atom_is_backbone( atomno ) );
}


void
Conformation::set_stub_transform(
  id::StubID const & stub_id1,
  id::StubID const & stub_id2,
  kinematics::RT const & target_rt
)
{
  set_dof_moved( atom_tree_->set_stub_transform( stub_id1, stub_id2, target_rt ) );
}

/// @brief  get the transform between two stubs
kinematics::RT
Conformation::get_stub_transform(
	id::StubID const & stub_id1,
	id::StubID const & stub_id2
) const
{
	return atom_tree_->get_stub_transform( stub_id1, stub_id2 );
}


void
Conformation::set_jump_atom_stub_id( id::StubID const& id )
{
  atom_tree_->set_jump_atom_stub_id( id );
}

kinematics::Stub
Conformation::stub_from_id( id::StubID const& id ) const {
	return atom_tree_->stub_from_id( id );
}


// The upstream and downstream Stubs are the coordinate frames between which this jump is transforming
kinematics::Stub
Conformation::upstream_jump_stub( int const jump_number ) const
{
	return atom_tree_->atom( jump_atom_id( jump_number ) ).get_input_stub();
}

// The upstream and downstream Stubs are the coordinate frames between which this jump is transforming
kinematics::Stub
Conformation::downstream_jump_stub( int const jump_number ) const
{
	return atom_tree_->atom( jump_atom_id( jump_number ) ).get_stub();
}

void
Conformation::update_polymeric_connection(
	Size const lower_seqpos,
	bool const update_connection_dep_atoms
) {
	if ( lower_seqpos < 1 || lower_seqpos >= size() ) return;

	Residue const & lower( *residues_[ lower_seqpos ] ), upper( *residues_[ lower_seqpos + 1 ] );

	bool const disconnected(
		lower.chain() != upper.chain() || lower.is_upper_terminus() ||
		upper.is_lower_terminus() || !lower.is_polymer() || !upper.is_polymer() ||
		lower.has_variant_type( "C_METHYLAMIDATION" ) ||
		upper.has_variant_type( "N_ACETYLATION" ) );


	if ( !disconnected ) set_polymeric_connection( lower_seqpos, lower_seqpos+1 );
	if (update_connection_dep_atoms) { //Update connection-dependent atoms, EVEN if there's no connection between lower_seqpos and lower_seqpos+1
		//TR << "Updating connection-dependent atoms." << std::endl; TR.flush();
		rebuild_polymer_bond_dependent_atoms( lower_seqpos    ,  1 );
		rebuild_polymer_bond_dependent_atoms( lower_seqpos + 1, -1 );
	} //else {
		//TR << "NOT updating connection-dependent atoms." << std::endl; TR.flush();
	//}
}

void
Conformation::update_noncanonical_connection(
                                             Size const lower_seqpos,
                                             Size const lr_conn_id,
                                             Size const upper_seqpos,
                                             Size const ur_conn_id)
{
    if ( lower_seqpos < 1 || lower_seqpos > size() ) return;
    if ( upper_seqpos < 1 || upper_seqpos > size() ) return;
    
    set_noncanonical_connection( lower_seqpos, lr_conn_id, upper_seqpos, ur_conn_id );
}

// identify polymeric connections
/// @details The lower residue connects to the upper residue through its upper_connect connection
/// The upper residue connects to the lower residue through its lower_connect connection
void
Conformation::set_polymeric_connection(
	Size res_id_lower,
	Size res_id_upper
)
{
	Size const lr_conn_id( residue_( res_id_lower ).type().upper_connect_id());
	Size const ur_conn_id( residue_( res_id_upper ).type().lower_connect_id());
	residues_[ res_id_lower ]->residue_connection_partner( lr_conn_id, res_id_upper, ur_conn_id );
	residues_[ res_id_upper ]->residue_connection_partner( ur_conn_id, res_id_lower, lr_conn_id );
}

void
Conformation::set_noncanonical_connection(
                                       Size res_id_lower,
                                       Size lr_conn_id,
                                       Size res_id_upper,
                                       Size ur_conn_id
                                       )
{
    residues_[ res_id_lower ]->residue_connection_partner( lr_conn_id, res_id_upper, ur_conn_id );
    residues_[ res_id_upper ]->residue_connection_partner( ur_conn_id, res_id_lower, lr_conn_id );
}



// Assigns disulfide bonds based on a pre-determined list
/// @note works in centroid and full-atom modes
void
Conformation::fix_disulfides(utility::vector1< std::pair<Size, Size> > disulf_bonds) {
	typedef std::pair<Size,Size> SizePair;
	BOOST_FOREACH(SizePair disulfide_bond, disulf_bonds){
		using utility::vector1;

		Size l_index = (disulfide_bond).first; //Lower residue
		Size u_index = (disulfide_bond).second; //Upper residue, usually l<u

		// Check that the residues exist
		if(l_index > size() ) {
			TR.Error << "[ERROR] Residue " << l_index << " is out of range." << std::endl;
			utility_exit();
		}
		if(u_index > size() ) {
			TR.Error << "[ERROR] Residue " << u_index << " is out of range." << std::endl;
			utility_exit();
		}

		//Swap the CYS for CYD
		bool replaced = conformation::change_cys_state(l_index, "CYD", *this);
		replaced = replaced && conformation::change_cys_state(u_index, "CYD", *this);
		if(! replaced) {
			TR.Error << "Failed to introduce CYD for disulfide ("
				<< l_index <<", "<< u_index << ")." << std::endl;
			continue;
		}

		//Next, form a bond between the two residues.
		//This is a little messy since we don't know the residue types of l & u.
		Residue const& l_res( residue( l_index ));
		Residue const& u_res( residue( u_index ));

		// Determine which atom forms the disulfide bond
		// Prefer SG to SG (fullatom) or CEN to CEN (centroid) bonds
		// Allow SG to CEN bonds, but give a warning
		// If neither atom is found (neither fullatom nor centroid) throw an error.
		Size l_bond_atom, u_bond_atom;

		bool l_has_sg  = l_res.type().has("SG");
		bool u_has_sg  = u_res.type().has("SG");
		bool l_has_cen = l_res.type().has("CEN");
		bool u_has_cen = u_res.type().has("CEN");

		if( l_has_sg ) {
			l_bond_atom = l_res.atom_index( "SG" );
			if(! u_has_sg)
				TR.Warning << "Bonding SG of residue " << l_index
					<< " to non-SG of residue "<< u_index << std::endl;
		}
		else if(l_has_cen) {
			l_bond_atom = l_res.atom_index( "CEN" );
		}
		else {
			TR.Error << "Cannot form disulfide bond with residue "<<l_index<< std::endl;
			continue;
		}

		if( u_has_sg ) {
			u_bond_atom = u_res.atom_index( "SG" );
			if(! l_has_sg)
				TR.Warning << "Bonding SG of residue " << u_index
					<< " to non-SG of residue "<< l_index << std::endl;
		}
		else if(u_has_cen) {
			u_bond_atom = u_res.atom_index( "CEN" );
			//			TR.Debug << "return index before res-type finalize " << u_bond_atom << std::endl;
			//residues_[ u_index ]->type().finalize();
			//			u_bond_atom = u_res.atom_index( "CEN" );
			//			TR.Debug << "return index after res-type finalize " << u_bond_atom << std::endl;
		}
		else {
			TR.Error << "Cannot form disulfide bond with residue "<<l_index<< std::endl;
			continue;
		}

		//Now have the correct atoms, so bond them
		Size l_connid = l_res.type().residue_connection_id_for_atom( l_bond_atom );
		Size u_connid = u_res.type().residue_connection_id_for_atom( u_bond_atom );

		residues_[ l_index ]->residue_connection_partner( l_connid, u_index, u_connid);
		residues_[ u_index ]->residue_connection_partner( u_connid, l_index, l_connid);

		assert( residue(l_index).has_variant_type( chemical::DISULFIDE ));
		assert( residue(u_index).has_variant_type( chemical::DISULFIDE ));

	} //done with this disulfide
}

// Detect existing disulfides from the protein structure.
/// @details For full atom confomations, looks at SG-SG distance. If the SG-SG
/// are about 2.02 A apart, calls it a disulfide bond. For centroid and other
/// conformations, the less accurate CB-CB distance is used instead. In this
/// case a CB-CB distance of 3.72 A is optimal.
/// @note Assumes full atom
void
Conformation::detect_disulfides()
{
	basic::ProfileThis doit( basic::CONFORMATION_DETECT_DISULF );
	using namespace graph;
	using namespace basic::options;

	// gather all cys, construct mapping from resid to cys index
	utility::vector1< Size > resid_2_cysid( size(), 0 );
	Size num_cys( 0 );
	for ( Size ii = 1; ii <= size(); ++ii ) {
		if ( residue(ii).aa() == chemical::aa_cys ) {
			++num_cys;
			resid_2_cysid[ ii ] = num_cys;
		}
	}
	if ( num_cys == 0 ) return;

	// construct reverse mapping from cys index to resid
	utility::vector1< Size > cysid_2_resid( num_cys );
	for ( Size ii = 1; ii <= size(); ++ii ) {
		if ( resid_2_cysid[ ii ] != 0 ) {
			cysid_2_resid[ resid_2_cysid[ ii ]] = ii;
		}
	}

	// If all the cys are fullatom, use stricter criteria
	bool fullatom(true);
	for( Size ii = 1; ii <= num_cys; ++ii) {
		if( residue_type(cysid_2_resid[ii]).residue_type_set().name()
				!= core::chemical::FA_STANDARD )
		{
			fullatom = false;
			break;
		}
	}
	// SG-SG distance for fullatom, CB-CB distance otherwise
	Real const typical_disulfide_distance = fullatom? 2.02 : 3.72;
	Real const tolerance = option[OptionKeys::in::detect_disulf_tolerance].user()
	                     ? option[OptionKeys::in::detect_disulf_tolerance]()
	                     : ( fullatom? 0.5 : 1.0 );
	std::string const distance_atom = fullatom? "SG" : "CB";

	// Create point graph
	PointGraphOP pg = new PointGraph;
	pg->set_num_vertices( num_cys );
	Distance maxrad( 0.0 );
	Distance maxd( typical_disulfide_distance + tolerance );
	for ( Size ii = 1; ii <= num_cys; ++ii ) {
		Residue const & ii_res = residue( cysid_2_resid[ ii ] );
		pg->get_vertex(ii).data().xyz() = ii_res.atoms()[ ii_res.nbr_atom() ].xyz();
		if ( ii_res.nbr_radius() > maxrad ) maxrad = ii_res.nbr_radius();
	}
	Distance neighbor_cutoff = maxrad + maxd;
	find_neighbors( pg, neighbor_cutoff );

	// Iterate across neighbors of cys residues; examine SG-SG distances.
	// Note that since graph only stores upper neighbors iterating through
	// the (upper) edge list will automatically prevent double counting.
	std::set< Size > processed_cys; // track cys that have already been processed
	for ( Size ii = 1; ii <= num_cys; ++ii ) {
		Size const ii_resid = cysid_2_resid[ ii ];
		//Size const ii_n_conn = residue( ii_resid ).type().n_residue_connections();
		Residue const & ii_res( residue( ii_resid ) );

		//if ii already processed, continue
		if( processed_cys.find( ii_resid) != processed_cys.end() ) {
			continue;
		}

		//Determine which atom makes the disulfide bond
		Size ii_sg_atomno(0);
		if(ii_res.type().has("SG")) {
			ii_sg_atomno = residue( cysid_2_resid[ ii ] ).atom_index( "SG" );
		} else if(ii_res.type().has("CEN")) {
			ii_sg_atomno = residue( cysid_2_resid[ ii ] ).atom_index( "CEN" );
		} else {
			TR.Error << "Error: Can't find an atom to disulfide bond from at residue "<< ii_resid <<std::endl;
			utility_exit();
		}

		Distance best_match( 0.0 );
		Size best_neighbor( 0 );
		//Size best_neighbor_cysid( 0 );

		for ( PointGraph::UpperEdgeListConstIter
				ii_iter = pg->get_vertex( ii ).upper_edge_list_begin(),
				ii_end_iter = pg->get_vertex( ii ).upper_edge_list_end();
				ii_iter != ii_end_iter; ++ii_iter ) {
			Size const jj = ii_iter->upper_vertex();

			Size const jj_resid = cysid_2_resid[ jj ];

			Residue const & jj_res( residue( jj_resid ) );

			//if jj already processed, continue
			if( processed_cys.find( jj_resid) != processed_cys.end() ) {
				continue;
			}

			Distance dist = ii_res.atom( ii_res.atom_index( distance_atom ) ).xyz().distance(
				jj_res.atom( jj_res.atom_index( distance_atom )).xyz() );
			if ( best_neighbor == 0 || dist < best_match ) {
				best_neighbor = jj_resid;
				//best_neighbor_cysid = jj;  // set but never used ~Labonte
				best_match = dist;
			}
		}

		if ( best_neighbor == 0 || best_match >= typical_disulfide_distance + tolerance ) {

			// handle case where old disulfide doesn't exist anymore and
			// needs to be cleared
			if ( processed_cys.find( ii_resid ) == processed_cys.end() && ii_res.has_variant_type( chemical::DISULFIDE ) ) {
				TR << "Reverting out-of-date disulfide CYD to CYS at resid " << ii_resid << std::endl;

				bool const successful_revert = conformation::change_cys_state( ii_resid, "CYS", *this ) && !residues_[ ii_resid ]->has_variant_type( chemical::DISULFIDE );
				if ( !successful_revert ) {
					TR.Error << "ERROR: unable to revert CYD to CYS for removal of disulfide at resid " << ii_resid << std::endl;
				}
			}

			// mark cys as processed
			processed_cys.insert( ii_resid );

			continue;

		} else { // found disulfide bond

			// Create disulfide bond using CYD residues.  Note that this will
			// end up doing a dummy replace for already existing disulfides,
			// but it doesn't necessarily hurt just in case something weird
			// happened.
			TR << "Found "<< (fullatom?"":"CEN ") << "disulfide between residues " << ii_resid << " " << best_neighbor << std::endl;
			TR << "current variant for " << ii_resid << " " << ( residues_[ ii_resid ]->has_variant_type( chemical::DISULFIDE ) ? "CYD" : "CYS" ) << std::endl;
			TR << "current variant for " << best_neighbor << " " << ( residues_[ best_neighbor ]->has_variant_type( chemical::DISULFIDE ) ? "CYD" : "CYS" ) << std::endl;

			bool const success_at_ii = conformation::change_cys_state( ii_resid, "CYD", *this )
				&& residues_[ ii_resid ]->has_variant_type( chemical::DISULFIDE );
			bool const success_at_best_neighbor = conformation::change_cys_state( best_neighbor, "CYD", *this )
				&& residues_[ best_neighbor ]->has_variant_type( chemical::DISULFIDE );

			if ( !success_at_ii ) {
				TR.Error << "ERROR: unable to create residue type CYD for disulfide at resid " << ii_resid << std::endl;
			}

			if ( !success_at_best_neighbor ) {
				TR.Error << "ERROR: unable to create residue type CYD for disulfide at resid " << best_neighbor << std::endl;
			}

			TR << "current variant for " << ii_resid << " " << ( residues_[ ii_resid ]->has_variant_type( chemical::DISULFIDE ) ? "CYD" : "CYS" ) << std::endl;
			TR << "current variant for " << best_neighbor << " " << ( residues_[ best_neighbor ]->has_variant_type( chemical::DISULFIDE ) ? "CYD" : "CYS" ) << std::endl;

			// Record SG-SG connections.
			if ( success_at_ii && success_at_best_neighbor ) {
				Residue const & ii_new_res( residue( ii_resid ) );
				// ASSUMPTION Disulfide forming cystein SG atoms for exactly one inter-residue chemical bond.
				Size ii_connid = ii_new_res.type().residue_connection_id_for_atom( ii_sg_atomno );
				Size jj_resid  = best_neighbor;
				Residue const & jj_res = residue( jj_resid );
				Size jj_sg_atomno(0);
				if(jj_res.type().has("SG")) {
					jj_sg_atomno = jj_res.atom_index( "SG" );
				} else if(jj_res.type().has("CEN")) {
					jj_sg_atomno = jj_res.atom_index( "CEN" );
				} else {
					TR.Error << "Error: Can't find an atom to disulfide bond from at residue "<< jj_resid <<std::endl;
					utility_exit();
				}

				Size jj_connid = jj_res.type().residue_connection_id_for_atom( jj_sg_atomno );

				residues_[ ii_resid ]->residue_connection_partner( ii_connid, jj_resid, jj_connid );
				residues_[ jj_resid ]->residue_connection_partner( jj_connid, ii_resid, ii_connid );
			}

			// mark both cys as processed
			processed_cys.insert( ii_resid );
			processed_cys.insert( best_neighbor );

		}
	}
}


// Conformation Cutting/Pasting //////////////////////////////////////////////////////////////////////////////////////

/// @details  Insert one conformation into another. Some tricky issues:
/// (1) residue connections: assume all residue connections within conf carry over, after renumbering residues positions
/// (2) jump numbers: see FoldTree::insert_fold_tree_by_jump
/// (3) chains: chain endings are inserted before insert_seqpos, and after insert_seqpos+conf.size()-1
///     ie. at the beginning and ending of the inserted conformation; all internal chain endings from conf are used
/// (4) atom_tree: setup_atom_tree is called to rebuild from scratch using the fold_tree
void
Conformation::insert_conformation_by_jump(
	Conformation const & conf,
	Size const insert_seqpos,         // rsd 1 in conf goes here
	Size const insert_jumppos,        // jump#1 in conf goes here
	Size const anchor_pos,            // in the current sequence numbering, ie before insertion of conf
	Size const anchor_jump_number,    // the desired jump number of the anchoring jump, default=0
	std::string const & anchor_atom,  // could be ""
	std::string const & root_atom     // ditto
)
{
	// ensure that residue data is OK
	pre_nresidue_change();

	// save some info
	Size const old_size( size() );
	Size const insert_size( conf.size() );
	Size const new_size( old_size + insert_size );
	assert( old_size );

	// sanity checks
	bool const fold_tree_polymer_bond( !fold_tree_->is_cutpoint( insert_seqpos-1 ) );
	bool const residues_polymer_bond( insert_seqpos > 1 && insert_seqpos <= old_size &&
																		residue( insert_seqpos-1 ).is_polymer_bonded( insert_seqpos ) );
	if ( fold_tree_polymer_bond || residues_polymer_bond ) {
		utility_exit_with_message("cant insert 'by_jump' into a polymer stretch");
	}

	// sequence numbering of existing residues, *_moved --> do this before changing residues_ array
	{
		utility::vector1< Size > old2new;
		for ( Size i=1; i<= old_size; ++i ) old2new.push_back( ( i >= insert_seqpos ) ? i+insert_size : i );
		update_sequence_numbering( new_size, old2new );
	}

	// insert the residues from conf
	{
		utility::vector1< Size > old2new;
		for ( Size i=1; i<= insert_size; ++i ) old2new.push_back( i+insert_seqpos-1 );

		utility::vector1< ResidueOP > new_residues;
		for ( Size i=1; i<= insert_size; ++i ) {
			ResidueOP new_rsd( conf.residue(i).clone() );
			new_rsd->update_sequence_numbering( old2new );
			Size const new_seqpos( old2new[i] );
			residues_.insert( residues_.begin() + new_seqpos-1, new_rsd );
			xyz_moved_[ new_seqpos ].clear(); xyz_moved_[ new_seqpos ].resize( new_rsd->natoms(), true );
			dof_moved_[ new_seqpos ].clear(); dof_moved_[ new_seqpos ].resize( new_rsd->natoms(), true );
			secstruct_.insert( secstruct_.begin() + ( new_seqpos - 1 ),  secstruct_[i] );
		}
	}


	// chains
	{
		utility::vector1< Size > new_chain_endings;
		for ( utility::vector1< Size >::const_iterator ch=chain_endings_.begin(); ch != chain_endings_.end(); ++ch ) {
			new_chain_endings.push_back( (*ch<insert_seqpos ? *ch : *ch + insert_size ) );
		}
		for ( Size i=1; i< conf.num_chains(); ++i ) new_chain_endings.push_back( conf.chain_end(i) + insert_seqpos-1 );
		if ( insert_seqpos > 1 &&
		     std::find(new_chain_endings.begin(),new_chain_endings.end(),insert_seqpos-1) == new_chain_endings.end())
			new_chain_endings.push_back( insert_seqpos-1 );
		if ( insert_seqpos <= old_size &&
		     std::find(new_chain_endings.begin(),new_chain_endings.end(),insert_seqpos+insert_size-1) == new_chain_endings.end())
			new_chain_endings.push_back( insert_seqpos+insert_size-1 );
		std::sort( new_chain_endings.begin(), new_chain_endings.end() );
		chain_endings_ = new_chain_endings;
		rederive_chain_ids();
	}

	// FoldTree
	fold_tree_->insert_fold_tree_by_jump( conf.fold_tree(), insert_seqpos, insert_jumppos, anchor_pos, anchor_jump_number,
																			 anchor_atom, root_atom );

	// AtomTree
	setup_atom_tree();

	// not sure if this is necessary, perhaps around the insertion? can't hurt though...
	residue_torsions_need_updating_ = true;

	// Notify length observers
	for (core::Size i = insert_seqpos; i < insert_size + insert_seqpos; i++)
	{
		notify_length_obs( LengthEvent( this, LengthEvent::RESIDUE_APPEND, i - 1, 1, &residue(i)), false );
	}
}


/// @details copy a stretch of coordinates/torsions from another pose
///  Fires IdentityEvent signals as residues are replaced.
void
Conformation::copy_segment(
	Size const size,
	Conformation const & src,
	Size const begin,
	Size const src_begin
)
{
	for ( Size i=0; i< size; ++i ) {
		Size const seqpos    (     begin + i );
		Size const seqpos_src( src_begin + i );
		replace_residue( seqpos, src.residue(seqpos_src), false );
	}
}

void
Conformation::insert_fragment(
	id::StubID const & instub_id,
	FragRT const & outstub_transforms,
	FragXYZ const & frag_xyz
)
{
	utility::vector1< AtomID > moved_atoms;

	atom_tree_->insert_fragment( instub_id, outstub_transforms, frag_xyz, moved_atoms );

	for ( Size i=1; i<= moved_atoms.size(); ++i ) {
		set_dof_moved( moved_atoms[i] );
	}
}


// DoFs/xyzs //////////////////////////////////////////////////////////////////////////////////////////////////////////

// Returns the AtomTree degree of freedom (DOF)  <id>
Real
Conformation::dof( DOF_ID const & id ) const
{
	return atom_tree_->dof( id );
}

// Sets the AtomTree degree of freedom (DOF)  <id>  to  <setting>
void
Conformation::set_dof( DOF_ID const & id, Real const setting )
{
	set_dof_moved( id );
	residue_torsions_need_updating_ = true; // might have been a torsion angle
	atom_tree_->set_dof( id, setting );
}


Real
Conformation::torsion( TorsionID const & tor_id ) const
{
	using numeric::conversions::degrees;

	if ( tor_id.type() == id::JUMP ) {
		// jump rigid-body offset
		return atom_tree_->dof( dof_id_from_torsion_id( tor_id ) );

	} else {
		// get torsion angle from the Residue
		//
		// note the use of the residue(...) access method
		// this will ensure that the torsions in the residue are in sync
		// with the atomtree torsions. Use of residues_[...]-> would not
		// be safe.
		//
		if ( tor_id.type() == id::BB ) {
			return residue( tor_id.rsd() ).mainchain_torsion( tor_id.torsion() );
		} else {
			return residue( tor_id.rsd() ).chi( tor_id.torsion() );
		}
	}
}

void
Conformation::set_torsion(
		TorsionID const & tor_id,
		Real const setting
)
{
	using numeric::conversions::radians;

	if ( tor_id.type() == id::JUMP ) {
		// jump rigid-body offset degree of freedom
		DOF_ID const dof_id( dof_id_from_torsion_id( tor_id ) );
		atom_tree_->set_dof( dof_id, setting );

		// update book-keeping to reflect that this torsion has changed
		set_dof_moved( dof_id );

	} else /* bb, chi, or nu */ {
		// update residue torsions
		switch (tor_id.type()) {
			case id::BB:
				residues_[tor_id.rsd()]->mainchain_torsions()[tor_id.torsion()] = setting;
				break;
			case id::CHI:
				residues_[tor_id.rsd()]->chi()[tor_id.torsion()] = setting;
				break;
			case id::NU:
				residues_[tor_id.rsd()]->nus()[tor_id.torsion()] = setting;
				break;
			default:
				TR.Error << "Unknown torsion type; no torsion set." << std::endl;
		}

		// find out what are the four atoms that define this torsion angle
		AtomID id1, id2, id3, id4;
		bool const fail( get_torsion_angle_atom_ids( tor_id, id1, id2, id3, id4 ) );

		if ( fail ) 	return;

		// AtomTree works in radians
		DOF_ID const dof_id( atom_tree_->set_torsion_angle( id1, id2, id3, id4, radians(setting) ) );

		if ( !dof_id.valid() ) {
			TR.Warning << "Unable to find torsion angle in atom_tree: " <<
				tor_id << std::endl;
			return;
		}

		// update book-keeping to reflect that this torsion has changed
		set_dof_moved( dof_id );
	}
}


// Returns the torsion angle defined by  <atom[1-4]>
Real
Conformation::torsion_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		AtomID const & atom4
) const
{
	return atom_tree_->torsion_angle( atom1, atom2, atom3, atom4 );
}

void
Conformation::set_torsion_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3,
	AtomID const & atom4,
	Real const setting,
	bool const quiet
)
{
	residue_torsions_need_updating_ = true;
	DOF_ID const dof_id( atom_tree_->set_torsion_angle( atom1, atom2, atom3, atom4, setting, quiet ) );
	if ( dof_id.valid() ) {
		set_dof_moved( dof_id );
	} else if ( !quiet ) {
		TR << "set_torsion_angle failed, unable to find dof_id: " << atom1 << ' ' << atom2 << ' ' << atom3 << ' ' <<
			atom4 << std::endl;
	}
}


// Returns the bond angle defined by  <atom[1-3]>
Real
Conformation::bond_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3
) const
{
	return atom_tree_->bond_angle( atom1, atom2, atom3 );
}

/// @note in radians!!
void
Conformation::set_bond_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3,
	Real const setting
)
{
	DOF_ID const dof_id( atom_tree_->set_bond_angle( atom1, atom2, atom3, setting ) );
	if ( dof_id.valid() ) {
		set_dof_moved( dof_id );
	}
}


// Returns the bond length between  <atom1>  and  <atom2>
Real
Conformation::bond_length(
		AtomID const & atom1,
		AtomID const & atom2
) const
{
	return atom_tree_->bond_length( atom1, atom2 );
}

void
Conformation::set_bond_length(
	AtomID const & atom1,
	AtomID const & atom2,
	Real const setting
)
{
	DOF_ID const dof_id( atom_tree_->set_bond_length( atom1, atom2, setting ) );
	if ( dof_id.valid() ) {
		set_dof_moved( dof_id );
	}
}


Conformation::Jump const &
Conformation::jump( int const jump_number ) const
{
	return atom_tree_->jump( jump_atom_id( jump_number ) );
}

// access a jump
Conformation::Jump const &
Conformation::jump( AtomID const & id ) const
{
	return atom_tree_->jump( id );
}

// Sets the jump  <jump_number>  to  <new_jump>
void
Conformation::set_jump(
		int const jump_number,
		Jump const & new_jump
)
{
	assert( new_jump.ortho_check() );
	AtomID const id( jump_atom_id( jump_number ) );
	atom_tree_->set_jump( id, new_jump );
	set_dof_moved( id );
}


void
Conformation::set_jump(
		AtomID const & id,
		Jump const & new_jump
)
{
	assert( new_jump.ortho_check() );
	atom_tree_->set_jump( id, new_jump );
	set_dof_moved( id );
}


// access xyz coordinates of an atom
PointPosition const &
Conformation::xyz( AtomID const & id ) const
{
	return atom_tree_->xyz( id );
}

void
Conformation::set_xyz(
	AtomID const & id,
	PointPosition const & position
)
{
	// update atomtree coords
	if ( !atom_tree_->empty() ) atom_tree_->set_xyz( id, position );

	// update residue coords
	residues_[ id.rsd() ]->set_xyz( id.atomno(), position );

	// notify scoring
	set_xyz_moved( id );
}

void
Conformation::batch_set_xyz(
	utility::vector1<AtomID> const & ids,
	utility::vector1<PointPosition> const & positions
)
{
	runtime_assert( ids.size() == positions.size() );

	// update atomtree coords
	if ( !atom_tree_->empty() ) atom_tree_->batch_set_xyz( ids, positions );

	// update residue coords
	for (core::Size i=1; i<=ids.size(); ++i)
		residues_[ ids[i].rsd() ]->set_xyz( ids[i].atomno(), positions[i] );

	// notify scoring
	set_xyz_moved( ids );
}

void
Conformation::batch_get_xyz(
	utility::vector1<AtomID> const & ids,
	utility::vector1<PointPosition> & positions
) const
{
	positions.resize( ids.size() );

	// update residue coords
	for (core::Size i=1; i<=ids.size(); ++i)
		positions[i] = xyz( ids[i] );
}


/// @details Set two bond angles and a bond length. DOES NOT DO ANY DIHEDRALS -- NOT EVEN OMEGA IF ITS A PROTEIN
void
Conformation::insert_ideal_geometry_at_polymer_bond( Size const seqpos )
{
	// use residuetype's rather than residues to avoid triggering a refold
	chemical::ResidueType const & lower_rsd_type( residue_type( seqpos   ) );
	chemical::ResidueType const & upper_rsd_type( residue_type( seqpos+1 ) );

	// what are the four atoms that define the bonds/angles?
	Size const nbb( lower_rsd_type.mainchain_atoms().size() );
	AtomID const atom1( lower_rsd_type.mainchain_atom( nbb - 1 ), seqpos   );
	AtomID const atom2( lower_rsd_type.mainchain_atom( nbb     ), seqpos   );
	AtomID const atom3( upper_rsd_type.mainchain_atom( 1       ), seqpos+1 );
	AtomID const atom4( upper_rsd_type.mainchain_atom( 2       ), seqpos+1 );

	// what is the ideal geometry?
	chemical::ResidueConnection const & connect1( lower_rsd_type.upper_connect() );
	chemical::ResidueConnection const & connect2( upper_rsd_type.lower_connect() );

	Real const bond_distance( connect1.icoor().d() );
	Real const bond_angle1( numeric::constants::d::pi - connect1.icoor().theta() );
	Real const bond_angle2( numeric::constants::d::pi - connect2.icoor().theta() );

	assert( ( connect1.icoor().stub_atom2().atomno() == Size( atom1.atomno() ) ) &&
					( connect1.icoor().stub_atom1().atomno() == Size( atom2.atomno() ) ) &&
					( connect2.icoor().stub_atom1().atomno() == Size( atom3.atomno() ) ) &&
					( connect2.icoor().stub_atom2().atomno() == Size( atom4.atomno() ) ) );

	set_bond_angle( atom1, atom2, atom3, bond_angle1 );
	set_bond_angle( atom2, atom3, atom4, bond_angle2 );
	set_bond_length( atom2, atom3, bond_distance );

	// unfortunately this appears to trigger a refold, in that it works in cartesian space to update the atom
	// positions.
	//
	// could be rewritten to communicate directly with the atomtree...
	rebuild_polymer_bond_dependent_atoms( seqpos );
}

/// @details Set two bond angles and a bond length. DOES NOT DO ANY DIHEDRALS -- NOT EVEN OMEGA IF ITS A PROTEIN
void
Conformation::insert_ideal_geometry_at_residue_connection( Size const pos1, Size const connid1 )
{
	// we use residues_[ xx ] rather than residue(xx) to avoid triggering refold/angle update
	//
	// determine what the other residue for this connection is:
	ResidueCOP rsd1( residues_[ pos1 ]() );
	Size const    pos2( rsd1->connect_map( connid1 ).resid() );
	Size const connid2( rsd1->connect_map( connid1 ).connid() );
	ResidueCOP rsd2( residues_[ pos2 ]() );

	// what is the ideal geometry?
	chemical::ResidueConnection const & connect1( rsd1->residue_connection( connid1 ) );
	chemical::ResidueConnection const & connect2( rsd2->residue_connection( connid2 ) );

	Real const bond_distance( connect1.icoor().d() );
	Real const bond_angle1( numeric::constants::d::pi - connect1.icoor().theta() );
	Real const bond_angle2( numeric::constants::d::pi - connect2.icoor().theta() );

	AtomID const atom1( connect1.icoor().stub_atom2().atomno(), pos1 );
	AtomID const atom2( connect1.icoor().stub_atom1().atomno(), pos1 );
	AtomID const atom3( connect2.icoor().stub_atom1().atomno(), pos2 );
	AtomID const atom4( connect2.icoor().stub_atom2().atomno(), pos2 );

	set_bond_angle( atom1, atom2, atom3, bond_angle1 );
	set_bond_angle( atom2, atom3, atom4, bond_angle2 );
	set_bond_length( atom2, atom3, bond_distance );

	// need to do something similar here:
	rebuild_residue_connection_dependent_atoms( pos1, connid1 );
	rebuild_residue_connection_dependent_atoms( pos2, connid2 );
}


void
Conformation::update_actcoords()
{
	for ( Size ii = 1; ii <= size(); ++ii )
	{
		if ( residues_[ ii ]->requires_actcoord() )
		{
			residues_[ ii ]->update_actcoord();
		}
	}
}

void
Conformation::update_actcoord( Size resid )
{
	if ( residues_[ resid ]->requires_actcoord() ) {
		residues_[ resid ]->update_actcoord();
	}
}

void
Conformation::update_orbital_coords( Size resid )
{
		residues_[ resid ]->update_orbital_coords();
}

void
Conformation::update_orbital_coords( Residue & rsd) const{
	BOOST_FOREACH(core::Size atom_with_orbitals, rsd.atoms_with_orb_index()){
		utility::vector1<core::Size> const & orbital_indices(rsd.bonded_orbitals(atom_with_orbitals));
		BOOST_FOREACH(core::Size orbital_index, orbital_indices){
			Vector orb_xyz(rsd.build_orbital_xyz(orbital_index));
			rsd.set_orbital_xyz(orbital_index, orb_xyz );
		}
	}
}


// ID access and conversions /////////////////////////////////////////////////////////////////////////////////////////

id::DOF_ID
Conformation::dof_id_from_torsion_id( TorsionID const & tor_id ) const
{
	if ( tor_id.type() == id::JUMP ) {
		// jump rigid-body offset degree of freedom
		int const rb_no( tor_id.torsion() );
		assert( rb_no >= 1 && rb_no <= 6 );
		int const jump_number( tor_id.rsd() );
		AtomID const id( jump_atom_id( jump_number ) );
		return DOF_ID( id, id::get_rb_type( rb_no ) );
	} else {
		// bb or chi
		// find out what are the four atoms that define this torsion angle
		AtomID id1, id2, id3, id4;
		bool const fail
			( get_torsion_angle_atom_ids( tor_id, id1, id2, id3, id4 ) );

		if ( fail ) {
			// probably a backbone torsion undefined b/c of a cutpoint
			return id::BOGUS_DOF_ID;
		}
		return atom_tree_->torsion_angle_dof_id( id1, id2, id3, id4 );
	}
}

id::AtomID
Conformation::jump_atom_id( int const jump_number ) const
{
	// the fold_tree may or may not contain the information necessary
	// to determine the exact AtomID of the corresponding JumpAtom

	kinematics::Edge const & edge( fold_tree_->jump_edge( jump_number ) );

	Size const seqpos( edge.stop() );

	if ( edge.has_atom_info() ) {
		return id::AtomID( residues_[ seqpos ]->atom_index( edge.downstream_atom() ), seqpos );
	} else {
		// get the default atomno used for Jump attachment
		int const atomno( get_root_atomno( residue_( seqpos ),
			kinematics::dir_jump));
		return AtomID( atomno, seqpos );
	}
}

///@note  Returns TRUE to signal FAILURE.
bool
Conformation::get_torsion_angle_atom_ids(
		TorsionID const & tor_id,
		AtomID & id1,
		AtomID & id2,
		AtomID & id3,
		AtomID & id4
) const
{
	uint const torsion( tor_id.torsion() );
	uint const seqpos( tor_id.rsd() );

	bool fail( false );

	switch (tor_id.type()) {
		case id::BB:
			// backbone torsion, eg if protein: 1==phi, 2==psi, 3==omega
			// may fail if we are at a chainbreak
			fail = backbone_torsion_angle_atoms(tor_id, id1, id2, id3, id4);
			break;
		case id::CHI:
			id1.rsd() = id2.rsd() = id3.rsd() = id4.rsd() = seqpos;
			id1.atomno() = residue_(seqpos).chi_atoms(torsion)[1];
			id2.atomno() = residue_(seqpos).chi_atoms(torsion)[2];
			id3.atomno() = residue_(seqpos).chi_atoms(torsion)[3];
			id4.atomno() = residue_(seqpos).chi_atoms(torsion)[4];
			break;
		case id::NU:
			id1.rsd() = id2.rsd() = id3.rsd() = id4.rsd() = seqpos;
			id1.atomno() = residue_(seqpos).nu_atoms(torsion)[1];
			id2.atomno() = residue_(seqpos).nu_atoms(torsion)[2];
			id3.atomno() = residue_(seqpos).nu_atoms(torsion)[3];
			id4.atomno() = residue_(seqpos).nu_atoms(torsion)[4];
			break;
		case id::JUMP:
			TR.Error << "Conformation::get_torsion_angle_atom_ids: " << tor_id <<
				" Jump 'torsions' are not described by four atoms!" << std::endl;
			utility_exit();
			break;
		default:
			TR.Error << "Invalid torsion type provided." << std::endl;
			fail = true;
	}

	return fail;
}


///@brief get two atoms connect by jump
bool
Conformation::get_jump_atom_ids(
  core::Size const jump_number,
	AtomID & upstream_id,
	AtomID & downstream_id ) const {
	downstream_id = jump_atom_id( jump_number );
	upstream_id   = atom_tree_->atom( downstream_id ).input_stub_atom1()->id();
	return false;
}


// for tracking changes to the structure /////////////////////////////////////////////////////////////////////////////

void
Conformation::update_domain_map( DomainMap & domain_map ) const
{
	domain_map.dimension( size() );
	domain_map = 0;
	atom_tree_->update_domain_map( domain_map, dof_moved_, xyz_moved_ );
}

/// @details called after domain map information is transferred to the Energies class
void
Conformation::reset_move_data()
{
	structure_moved_ = false;
	dof_moved_.fill_with( false );
	xyz_moved_.fill_with( false );
}


// observer management ///////////////////////////////////////////////////////////////////////////////////////////////

// clear all observers
/// @remarks ConnectionEvent::DISCONNECT will be sent to all observers
void Conformation::clear_observers() {
	// send the signal first
	notify_connection_obs( ConnectionEvent( this, ConnectionEvent::DISCONNECT ) );

	// clear all hubs
	xyz_obs_hub_.clear();
	length_obs_hub_.clear();
	identity_obs_hub_.clear();
	general_obs_hub_.clear();
	connection_obs_hub_.clear();
}

// fire a ConnectionEvent::TRANSFER to transfer observers from some source Conformation
/// @param src Take observers from this source Conformation.
/// @remarks Only observers that properly honor the TRANSFER event by
///  re-attaching themselves to 'this' Conformation will be transferred.
void Conformation::receive_observers_from( Conformation const & src ) {
	// tell the observers in source Conformation that a transfer is
	// occurring and they need to re-register themselves with 'this'
	// conformation
	if ( this == &src ) return; // do not transfer observers from yourself
	src.notify_connection_obs( ConnectionEvent( this, ConnectionEvent::TRANSFER ) );
}


// additional observer behavior //////////////////////////////////////////////////////////////////////////////////////

// wait for stdin after sending a GeneralEvent signal
void
Conformation::debug_pause( bool const flag ) const
{
	if ( flag ) {
		general_obs_hub_.pause();
	} else {
		general_obs_hub_.unpause();
	}
}

// waiting for stdin after sending a GeneralEvent signal?
bool
Conformation::debug_pause() const
{
	return general_obs_hub_.pausing();
}


// signal management /////////////////////////////////////////////////////////////////////////////////////////////////

// block signals from being sent and buffer them to be sent after unblocking
void
Conformation::buffer_signals()
{
	general_obs_hub_.buffer();
	identity_obs_hub_.buffer();
	length_obs_hub_.buffer();
	xyz_obs_hub_.buffer();
}


// block signals from being sent
/// @warning for safety, ConnectionEvents are never blocked
void
Conformation::block_signals()
{
	general_obs_hub_.block();
	identity_obs_hub_.block();
	length_obs_hub_.block();
	xyz_obs_hub_.block();
}


// allow signals to be sent
/// @details If unblocking after buffering, buffered/held signals will be sent.
void
Conformation::unblock_signals()
{
	general_obs_hub_.unblock();
	identity_obs_hub_.unblock();
	length_obs_hub_.unblock();
	xyz_obs_hub_.unblock();
}


// are signals being blocked and buffered?
bool
Conformation::buffering_signals() const
{
	return general_obs_hub_.buffering();
}


// are signals being blocked?
bool
Conformation::blocking_signals() const
{
	return general_obs_hub_.blocked();
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// private
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// @details private
/// need to update:
/// 1. seqpos of Residues, residue_connection_partners_ info in Residues
/// 2. xyz_moved, dof_moved
void
Conformation::update_sequence_numbering(
	Size const new_size,
	utility::vector1< Size > const & old2new
)
{
	// residues:
	for ( Size i=1; i<= size(); ++i ) {
		// seqpos and residue_connection_partners
		residues_[i]->update_sequence_numbering( old2new );
	}

	// chain_endings not done here -- use rederive_chain_endings after the residues_ array is OK

	// xyz_moved, dof_moved
	xyz_moved_.update_sequence_numbering( new_size, old2new );
	dof_moved_.update_sequence_numbering( new_size, old2new );
}


void
Conformation::rebuild_polymer_bond_dependent_atoms( Size const seqpos, int const upper_lower )
{
	// assumes no interpendencies among the atoms being built
	update_residue_coordinates();
	Residue const & rsd( residue_(seqpos) );

 	for ( Size i=1, ie=rsd.natoms(); i<= ie; ++i ) {
		// OP1 and OP2 could in principle be figured out to 'depend on polymer lower'; they
		// currently depend on atoms which themselves can depend on LOWER. Following has hacks
		// until a general fix is tested and checked in. --rhiju.
 		if ( ( upper_lower == -1 && rsd.icoor(i).depends_on_polymer_lower() ) ||
				 ( upper_lower == -1 && rsd.atom_name( i ) == " OP1" ) || // hack1
				 ( upper_lower == -1 && rsd.atom_name( i ) == " OP2" ) || // hack2
 				 ( upper_lower ==  1 && rsd.icoor(i).depends_on_polymer_upper() ) ) {
 			set_xyz( AtomID(i,seqpos), rsd.icoor(i).build( rsd, *this ) );
 		}
	}
}


/// @details insert a polymer residue
///  Fires a LengthEvent::RESIDUE_PREPEND signal.
void
Conformation::insert_polymer_residue(
	Residue const & new_rsd_in,
	Size const seqpos, // desired seqpos of new_rsd
	bool const join_lower,
	bool const join_upper
)
{
	pre_nresidue_change();
	// debug termini status
	if ( join_lower ) assert( !new_rsd_in.is_lower_terminus() );
	if ( join_upper ) assert( !new_rsd_in.is_upper_terminus() );

	// this handles all renumbering internal to the Residues, *_moved arrays
	residues_insert( seqpos, new_rsd_in, !join_upper );

	Residue const & new_rsd( residue_( seqpos ) );

	fold_tree_->insert_polymer_residue( seqpos, join_lower, join_upper );

	insert_residue_into_atom_tree( new_rsd, *fold_tree_, const_residues(), *atom_tree_ );

	residue_torsions_need_updating_ = true; // could reupdate before and after

	assert( atom_tree_->size() == size() && Size(fold_tree_->nres()) == size() );

	notify_length_obs( LengthEvent( this, LengthEvent::RESIDUE_PREPEND, seqpos, 1, &new_rsd ), false );
}


/// @details add a residue into residues_ container, update its seqpos, chainid as well
/// fold tree and atoms.
/// private now -- this is our internal routine with everything determined ahead of time
/// root_atomno and anchor_id.atomno() may not be filled in yet
/// Fires a LengthEvent::RESIDUE_APPEND signal.
void
Conformation::append_residue(
	conformation::Residue const & new_rsd_in,
	bool const attach_by_jump,
	std::string const& root_atom,
	id::NamedAtomID anchor_id,
	bool const start_new_chain
)
{
	pre_nresidue_change();

	Size const seqpos( size() + 1 );
	assert( seqpos == fold_tree_->nres() + 1 );

	// is this the first residue?
	bool const first_residue( seqpos == 1 );

	// append to residues
	residues_append( new_rsd_in, start_new_chain, attach_by_jump, root_atom, anchor_id );

	// get reference to new rsd
	Residue const & new_rsd( *residues_[ seqpos ] );

	int anchor_atomno( 0 );

	// If there is an anchor atom specified, find its index number in the residue.
	if ( anchor_id.atom().size() ) anchor_atomno = residues_[ anchor_id.rsd() ]->atom_index( anchor_id.atom() );
	// update the fold_tree
	if ( first_residue ) {
		fold_tree_->simple_tree( 1 );
	} else {
		int root_atomno = 0;
		if ( root_atom.size() ) root_atomno = new_rsd.atom_index( root_atom );

		// try to detect a non-polymer chemical bond:
		if ( !attach_by_jump && root_atomno != 0 && ( new_rsd.is_ligand() || // needed because new_rsd.is_lower_terminus() asserts this is a polymer
				new_rsd.is_lower_terminus() || anchor_id.rsd() != seqpos-1 ||
					residues_[ seqpos-1 ]->is_upper_terminus() ||      // if we got here we're connect to seqpos-1...
					root_atomno != int(new_rsd_in.lower_connect_atom()) ||  // ...by our lower-connect-atom
					anchor_atomno != int(residues_[ seqpos-1 ]->upper_connect_atom()) ) ) {
			// must be a chemical bond since the criteria for a polymer connection are not met
			TR << "appending residue by a chemical bond in the foldtree: " << seqpos << ' ' <<
				new_rsd.name() << " anchor: " << anchor_id << " root: " << root_atom <<  std::endl;
            
			fold_tree_->append_residue_by_chemical_bond( anchor_id.rsd(), anchor_id.atom(), root_atom );
		} else {
			fold_tree_->append_residue( attach_by_jump, anchor_id.rsd(), anchor_id.atom(), root_atom );
		}
	}
	TR.Trace << "CURRENT_" << *fold_tree_ << std::endl;

	// update the atom_tree.
	// alternatively we could call setup_atom_tree, would be more expensive.
	// note that, inside these routines, new_rsd's sequence number is being used.
	// good thing we set it already.

	if ( first_residue ) {
		assert( atom_tree_->empty() );
		setup_atom_tree(); // just this once

	} else {
		insert_residue_into_atom_tree( new_rsd, *fold_tree_, const_residues(), *atom_tree_ );
	}

	residue_torsions_need_updating_ = true;
	//if ( !residue_torsions_need_updating_ ) update_residue_torsions( seqpos );

	notify_length_obs( LengthEvent( this, LengthEvent::RESIDUE_APPEND, seqpos - 1, 1, &new_rsd ), false );
} // append_residue


/// @details PRIVATE: wrap direct access to the Residues container for replacement
///  also handles redimensioning and setting of the _moved arrays
/// @warning Since this is the innermost routine when replacing a residue, it does not
///  fire any signals to ensure that observer access to the residues container does not
///  trigger an atom tree refold prematurely.
void
Conformation::residues_replace(
	Size const seqpos,
	Residue const & new_rsd
)
{

	int const old_chain = residues_[ seqpos ]->chain();
	ResidueOP old_residue = residues_[ seqpos ];
	residues_[ seqpos ] = new_rsd.clone();
	residues_[ seqpos ]->seqpos( seqpos );
	residues_[ seqpos ]->chain( old_chain );

	residues_[ seqpos ]->copy_residue_connections( *old_residue );
	//Loop through all the connections of the new residue and ensure that the residues connected to it have their
	//connect_map_ updated appropriately:
	for(core::Size i=1, imax=residues_[seqpos]->type().n_residue_connections(); i<=imax; ++i)
	{
		if(!residues_[seqpos]->connection_incomplete(i)) { //If we're connected to something,
			core::Size const res_to_update = residues_[seqpos]->connected_residue_at_resconn(i); //Get the index of the residue to update
			residues_[res_to_update]->update_connections_to_other_residue(*(residues_[seqpos])); //Update the connections
		}
	}

	// mark this position as having changed, resize the _moved arrays
	structure_moved_ = true;

	utility::vector1< bool > & xyz_m( xyz_moved_[seqpos] );
	utility::vector1< bool > & dof_m( dof_moved_[seqpos] );
	xyz_m.clear(); // guarantee that they're all set to true
	xyz_m.resize( new_rsd.natoms(), true );
	{ // this is a little tricky:
		// on the one hand, we dont want to obliterate any existing dof_moved info
		// on the other, if we are just replacing a residue, setting dof_moved to TRUE for this position
		// would likely lead to a massive score recalculation. So here's a compromise
		bool any_dof_moved( false );
		for ( Size i=1; i<= dof_m.size(); ++i ) {
			if ( dof_m[ i ] ) {
				any_dof_moved = true;
				break;
			}
		}
		dof_m.clear();
		dof_m.resize( new_rsd.natoms(), any_dof_moved );
	}
}


/// @details PRIVATE: wrap direct access to the Residues container for replacement
/// also handles redimensioning and setting of the _moved arrays and updating of sequence numbers
/// @warning Since this is the innermost routine when replacing a residue, it does not
///  fire any signals to ensure that observer access to the residues container does not
///  trigger an atom tree refold prematurely.
void
Conformation::residues_insert(
	Size const seqpos,
	Residue const & new_rsd,
	bool const use_lower_chain, // = false,
	bool const new_chain // = false
)
{
	assert( ! new_chain || residues_[ seqpos-1 ]->chain() != residues_[ seqpos ]->chain() );

	Size const old_size( residues_.size() ), new_size( old_size+1 );

	// first renumber things
	utility::vector1< Size > old2new( old_size, 0 );
	for ( Size i=1; i<= old_size; ++i ) {
		if ( i< seqpos ) old2new[i] = i;
		else old2new[i] = i+1;
	}
	update_sequence_numbering( new_size, old2new ); // numbering in residues_ and *_moved

	Size const old_chain( ( new_chain || use_lower_chain || Size(seqpos) == new_size ) ?
		residues_[ seqpos-1 ]->chain() :
		residues_[ seqpos ]->chain() );

	Size const newrsd_chain( new_chain ? old_chain + 1 :old_chain );

	residues_.insert( residues_.begin() + (seqpos-1), new_rsd.clone() );
	assert( residues_[seqpos]->name() == new_rsd.name() );
	residues_[ seqpos ]->seqpos( seqpos );
	residues_[ seqpos ]->chain( newrsd_chain );
	if ( new_chain ) {
		for ( Size ii = seqpos+1; ii <= new_size; ++ii ) {
			residues_[ ii ]->chain( residues_[ ii ]->chain() + 1 );
		}
	}
	// wipe residue connection data that may have been cloned from the original residue
	residues_[ seqpos ]->clear_residue_connections();

	// rederive polymeric connection status from chain id's or termini status
	update_polymeric_connection( seqpos-1 );
	update_polymeric_connection( seqpos );

	// recompute the chain_endings_ array from residues_[]->chain()
	rederive_chain_endings();

	// mark this position as having changed, resize the _moved arrays
	structure_moved_ = true;

	utility::vector1< bool > & xyz_m( xyz_moved_[seqpos] );
	utility::vector1< bool > & dof_m( dof_moved_[seqpos] );
	assert( xyz_m.empty() && dof_m.empty() );
	xyz_m.resize( new_rsd.natoms(), true  );
	dof_m.resize( new_rsd.natoms(), false );

	secstruct_.insert( secstruct_.begin() + ( seqpos - 1 ), 'L' );
}


/// @details PRIVATE: wrap direct access to the Residues container for replacement
///  also handles redimensioning and setting of the _moved arrays and updating of sequence numbers
/// @warning Since this is the innermost routine when replacing a residue, it does not
///  fire any signals to ensure that observer access to the residues container does not
///  trigger an atom tree refold prematurely.
void
Conformation::residues_delete(
	Size const seqpos
)
{
	Size const old_size( residues_.size() ), new_size( old_size-1 );

	// delete from residues_
	residues_.erase( residues_.begin()+seqpos-1 );

	// now renumber things
	utility::vector1< Size > old2new( old_size, 0 );
	for ( Size i=1; i<= old_size; ++i ) {
		if ( i < seqpos ) old2new[i] = i;
		else if ( i > seqpos ) old2new[i] = i-1;
	}
	update_sequence_numbering( new_size, old2new ); // numbering in residues_ and *_moved

	// rederive polymeric connection status from chain id's or termini status
	if ( seqpos>1 ) update_polymeric_connection( seqpos-1 );

	// recompute the chain_endings_ array from residues_[]->chain()
	rederive_chain_endings();
	rederive_chain_ids(); // necessary if we deleted an entire, single-residue chain

	// update arrays
	secstruct_.erase( secstruct_.begin() + ( seqpos - 1 ) );
}


/// @details PRIVATE: wrap direct access to the Residues container for appending
/// @warning Since this is the innermost routine when replacing a residue, it does not
///  fire any signals to ensure that observer access to the residues container does not
///  trigger an atom tree refold prematurely.
void
Conformation::residues_append( Residue const & new_rsd, bool const start_new_chain, bool const by_jump, std::string const & root_atom, id::NamedAtomID anchor_id)
{
	residues_.push_back( new_rsd.clone() );
	// ensure that the residue number is set
	Size const nres( residues_.size() );
	residues_[ nres ]->seqpos( nres );

	// apl If this is the first residue overall, set the chain id to 1.
	// Otherwise, set the chain id based on whether this a new chain or not.
	if ( nres > 1 ) {
		if ( start_new_chain ) {
			residues_[ nres ]->chain( residues_[ nres - 1 ]->chain() + 1 );
			TR.Debug << "Starting a new chain: chain " << residues_[nres - 1]->chain() << std::endl;
		} else {
			residues_[ nres ]->chain( residues_[ nres - 1 ]->chain() );
		}
		rederive_chain_endings();
		rederive_chain_ids();
	} else {
		residues_[ nres ]->chain( 1 );
	}
	// wipe residue connection data that may have been cloned from the original residue
	residues_[ nres ]->clear_residue_connections();

	// update polymeric connection status from chain id's and termini status (but only if it's possible they're chemically connected)
	if ( nres > 1 && ! by_jump ) {
        if (anchor_id == id::BOGUS_NAMED_ATOM_ID) {
            update_polymeric_connection(nres - 1);
        }
        else {
            Size root_atomno = residues_[nres]->type().atom_index(root_atom);
            Size lr_conn_id = residues_[nres]->type().residue_connection_id_for_atom(root_atomno);
            Size anchor_atomno = residues_[anchor_id.rsd()]->atom_index( anchor_id.atom() );
            Size ur_conn_id = residues_[anchor_id.rsd()]->type().residue_connection_id_for_atom(anchor_atomno);
            if ( anchor_id.rsd() == nres - 1
                && lr_conn_id == residues_[nres]->type().lower_connect_id()
                && ur_conn_id == residues_[anchor_id.rsd()]->type().upper_connect_id()
                )
            {
                update_polymeric_connection(nres - 1);
            }
            else {
                update_noncanonical_connection(nres, lr_conn_id, anchor_id.rsd(), ur_conn_id);
            }
        }
    }

	// have to calculate scores with this guy
	structure_moved_ = true;
	xyz_moved_.resize( nres );
	dof_moved_.resize( nres );
	xyz_moved_.resize( nres, new_rsd.natoms(), true  );
	dof_moved_.resize( nres, new_rsd.natoms(), false );

	secstruct_.resize( nres, 'L' );
}


/// @details used when we want to copy torsions from the atomtree to the residues
/// used when updating the torsion angles stashed in the Residues
///	other torsion access calls look inside the residues since this
///	will be faster than getting from the atomtree, at least in the
///	current implementation of this routine
Real
Conformation::atom_tree_torsion( TorsionID const & tor_id ) const
{
	using numeric::conversions::degrees;

	if ( tor_id.type() == id::JUMP ) {
		// jump rigid-body offset
		return atom_tree_->dof( dof_id_from_torsion_id( tor_id ) );

	} else {
		// bb or chi
		// find out what are the four atoms that define this torsion angle
		AtomID id1, id2, id3, id4;
		bool const fail
			( get_torsion_angle_atom_ids( tor_id, id1, id2, id3, id4 ) );

		if ( fail ) {
			return 0.0;
		}

		//Check whether the four atoms are a valid DOF id:
		if(atom_tree_->torsion_angle_dof_id(id1,id2,id3,id4, true).valid()) {
			//printf("Torsion angle %lu-%s %lu-%s %lu-%s %lu-%s is valid.\n",
			//	id1.rsd(), residue_(id1.rsd()).atom_name(id1.atomno()).c_str(), 
			//	id2.rsd(), residue_(id2.rsd()).atom_name(id2.atomno()).c_str(),
			//	id3.rsd(), residue_(id3.rsd()).atom_name(id3.atomno()).c_str(),
			//	id4.rsd(), residue_(id4.rsd()).atom_name(id4.atomno()).c_str() ); fflush(stdout); //DELETE ME
			// atomtree works in radians
			return degrees( atom_tree_->torsion_angle( id1, id2, id3, id4 ) );
		} else { //If the torsion is not in the atom tree, it might still be a valid torsion angle (e.g. an C-terminus-to-N-terminus peptide bond).
			//printf("Torsion angle %lu-%s %lu-%s %lu-%s %lu-%s is NOT valid.\n",
			//	id1.rsd(), residue_(id1.rsd()).atom_name(id1.atomno()).c_str(), 
			//	id2.rsd(), residue_(id2.rsd()).atom_name(id2.atomno()).c_str(),
			//	id3.rsd(), residue_(id3.rsd()).atom_name(id3.atomno()).c_str(),
			//	id4.rsd(), residue_(id4.rsd()).atom_name(id4.atomno()).c_str() ); fflush(stdout); //DELETE ME
			//Determine whether the atoms are bonded in a way that defines a dihedral angle:
			bool const bonded = atoms_are_bonded(id1, id2) && atoms_are_bonded(id2, id3) && atoms_are_bonded(id3, id4);
			if(bonded) {
				//Calculate the dihedral angle between these atoms and return it.
				core::Real this_torsion(0.0);
				numeric::dihedral_degrees(
					residue_(id1.rsd()).xyz(id1.atomno()),
					residue_(id2.rsd()).xyz(id2.atomno()),
					residue_(id3.rsd()).xyz(id3.atomno()),
					residue_(id4.rsd()).xyz(id4.atomno()),
					this_torsion
				);
				return this_torsion;
			}
		}
	}
	return 0.0; //DEFAULT CASE, if the angle isn't in the atom tree AND the four atoms aren't bonded.
}


void
Conformation::setup_atom_tree()
{
	// this owns the tree
	kinematics::AtomPointer2D atom_pointer;
	build_tree( *fold_tree_, const_residues(), atom_pointer );

	// replace the current data in the atom_tree
	atom_tree_->replace_tree( atom_pointer );

	residue_coordinates_need_updating_ = false;
	residue_torsions_need_updating_ = true;
	update_residue_torsions();
}


// messy... could perhaps keep this cached if speed is an issue
// but then have to keep sync'ed with changes in the residues
//
// I don't think this is going to be performance critical
//
// consider that in current rosetta we routinely re-evaluate all the
// backbone and sidechain torsion angles from the xyz coordinates...
//
// This was hacky to begin with, but I feel that I have made it more hacky now.
// I am currently working on a more elegent method to handle residue types that
// have a a non-standard number of backbone torsions (ie. beta-peptides, beta-peptoids,
// and the dipeptides ACE-X-NME I use for making rotamers) --Doug
//
// Dear 2010 Doug,
// Sorry, I have completely forgoton what we were talking about in your comment. Very Fermat
// of you. I have a few vauge ideas about what this miraculous solution could possibly be
// but we probably had more brain cells then so who knows. Basically defining the atoms
// that comprise the backbone torsions in the params files themselves and storing that in
// the ResidueType, like is currently done for the chi angles. I guess it is just the cacheing
// idea in the comments before mineThat way, pathces could modify which atoms make up backbone
// torsions and we would not need special if checks to look at varient types if we want to
// have 4 backbone torsions. Until then however I am just going to make this function more
// hacky, sorry.
// Love, 2012 Doug
//
/// @note returns TRUE for FAILURE
bool
Conformation::backbone_torsion_angle_atoms(
	TorsionID const & id,
	AtomID & id1,
	AtomID & id2,
	AtomID & id3,
	AtomID & id4
) const
{

	//std::cout << "here" << std::endl;
	using chemical::AtomIndices;

	bool fail( true ); // return value

	Size const seqpos( id.rsd() );
	Size const torsion( id.torsion() );

	Residue const & rsd( *residues_[ seqpos ] );

	AtomIndices const & mainchain( rsd.mainchain_atoms() );

	Size const ntorsions( mainchain.size() );
	assert( torsion >= 1 && torsion <= ntorsions );

	// this is hacky
	// the ACETYLATED_NTERMINUS and METHYLATED_CTERMINUS prepend and append additional backbone atoms which is why the numbers may seem off
	if( rsd.has_variant_type( chemical::ACETYLATED_NTERMINUS ) && rsd.has_variant_type( chemical::METHYLATED_CTERMINUS ) ) {
		// set all id rsds to seqpos since they are all in the same residue
		id1.rsd() = id2.rsd() = id3.rsd() = id4.rsd() = seqpos;

		// torsion 1; epsilon
		if( torsion == 1 ){
			id1.atomno() = rsd.atom_index( "CP2" );
			id2.atomno() = mainchain[1];
			id3.atomno() = mainchain[2];
			id4.atomno() = mainchain[3];
		}
		// torsion 2; phi
		if( torsion == 2 ){
			id1.atomno() = mainchain[1];
			id2.atomno() = mainchain[2];
			id3.atomno() = mainchain[3];
			id4.atomno() = mainchain[4];
		}
		// torsion 3; psi
		if( torsion == 3 ){
			id1.atomno() = mainchain[2];
			id2.atomno() = mainchain[3];
			id3.atomno() = mainchain[4];
			id4.atomno() = rsd.atom_index( "NM" );
		}
		// torsion 4; omega
		if( torsion == 4 ){
			id1.atomno() = mainchain[3];
			id2.atomno() = mainchain[4];
			id3.atomno() = rsd.atom_index( "NM" );
			id4.atomno() = rsd.atom_index( "CN" );
		}

		// the ACETYLATED_NTERMINUS prepends an additional backbone atom which is why the numbers are increased by one
	} else if ( rsd.has_variant_type( chemical::ACETYLATED_NTERMINUS ) ) {
		// torsion 1; phi
		if( torsion == 1 ) {
			//TR << "HI MOM, ACET_NTERM 1" << std::endl;
			id1.rsd() = seqpos; id1.atomno() = mainchain[1];
			id2.rsd() = seqpos; id2.atomno() = mainchain[2];
			id3.rsd() = seqpos; id3.atomno() = mainchain[3];
			id4.rsd() = seqpos; id4.atomno() = mainchain[4];
		}
		// torsion 2; psi
		else if( torsion == 2 ) {
			if ( seqpos+1 <= residues_.size() ) { // should this be the chain end?
				AtomIndices const & next_mainchain( residue_( seqpos+1 ).mainchain_atoms() );
				id1.rsd() = seqpos;   id1.atomno() = mainchain[2];
				id2.rsd() = seqpos;   id2.atomno() = mainchain[3];
				id3.rsd() = seqpos;   id3.atomno() = mainchain[4];
				id4.rsd() = seqpos+1; id4.atomno() = next_mainchain[1];
			} else {
				return true; // torsion not well defined
			}
		}
		// torsion 3; omg
		else if( torsion == 3 ) {
			if ( seqpos+1 <= residues_.size() ) { // should this be the chain end?
				AtomIndices const & next_mainchain( residue_( seqpos+1 ).mainchain_atoms() );
				id1.rsd() = seqpos;   id1.atomno() = mainchain[3];
				id2.rsd() = seqpos;   id2.atomno() = mainchain[4];
				id3.rsd() = seqpos+1; id3.atomno() = next_mainchain[1];
				id4.rsd() = seqpos+1; id4.atomno() = next_mainchain[2];
			} else {
				return true; // torsion not well defined
			}
		}
		// catch other stuff
		else {
			return true; // torsion not well defined
		}
	}
	else {

	///////////////////////////////////////////
	// first atom -- may be in seqpos-1
	if ( torsion > 1 ) {
		id1.rsd()    = seqpos;
		id1.atomno() = mainchain[ torsion-1 ];
	} else {
		if ( rsd.has_variant_type( chemical::NTERM_CONNECT ) && residues_[ chain_end( rsd.chain() ) ]->has_variant_type( chemical::CTERM_CONNECT ) ) {
			// assume that they are connected as a cyclic polymer. phi is defined even though it is first residue, only used for torsion 1
			//TR << "HI MOM, NTERM_CONN 1" << std::endl;
			AtomIndices const & cyclic_partner_mainchain( residue_( chain_end( rsd.chain() ) ).mainchain_atoms() );
			id1.rsd() = chain_end( rsd.chain() ); // last residue in chain
			id1.atomno() = cyclic_partner_mainchain[ cyclic_partner_mainchain.size() ]; // last mainchain atom in last residue in chain
		} else if ( fold_tree_->is_cutpoint( seqpos-1 ) ) { // seems like this should be a bug if seqpos==1
			if ( rsd.has_variant_type( chemical::CUTPOINT_UPPER ) ) {
				id1.rsd() = seqpos;
				id1.atomno() = rsd.atom_index( "OVU1" );
			} else if ( rsd.has_variant_type( "N_ACETYLATION" ) ) {
				id1.rsd() = seqpos;
				id1.atomno() = rsd.atom_index( " CP " );
			} else if ( rsd.has_variant_type( "FIVE_PRIME_PHOSPHATE" ) ){
				id1.rsd() = seqpos;
				id1.atomno() = rsd.atom_index( "XO3'" );
			} else if ( seqpos==1 /*<-- only necessary for this case?*/ && rsd.has_lower_connect() && !rsd.connection_incomplete( rsd.type().lower_connect_id() ) ) {
				id1.rsd() = rsd.residue_connection_partner( rsd.type().lower_connect_id() ); //Get the residue index of the residue connected to this residue at this residue's lower connection.
				if(residue_(id1.rsd()).connect_map_size() < rsd.residue_connection_conn_id( rsd.type().lower_connect_id() ) ) {
					return true; //FAIL if this residue is not connected properly.
				}
				id1.atomno() =  residue_(id1.rsd()).residue_connect_atom_index( rsd.residue_connection_conn_id( rsd.type().lower_connect_id() ) ); //Get the atom index in the connected residue of the atom that's making a connection to THIS residue's connection #1.  (Convoluted, I know.)
			} else {
				// first bb-torsion is not well-defined
				return true; // FAILURE
			}
		} else {
			//AtomIndices const & prev_mainchain
			//	( residue_( seqpos-1 ).mainchain_atoms() );
			//id1.rsd()    = seqpos-1;
			//id1.atomno() = prev_mainchain[ prev_mainchain.size() ];

			//Altered by VKM, 12 June 2014: we want to fish out whatever atom the residue is connected to.
			if(!rsd.has_lower_connect() || rsd.connection_incomplete( rsd.type().lower_connect_id() )) {
				return true; //FAIL if this residue is not connected to anything.
			}
			id1.rsd() = rsd.residue_connection_partner( rsd.type().lower_connect_id() ); //Get the residue index of the residue connected to this residue at this residue's lower connection.
			if(residue_(id1.rsd()).connect_map_size() < rsd.residue_connection_conn_id( rsd.type().lower_connect_id() ) ) {
				return true; //FAIL if this residue is not connected properly.
			}
			id1.atomno() =  residue_(id1.rsd()).residue_connect_atom_index( rsd.residue_connection_conn_id( rsd.type().lower_connect_id() ) ); //Get the atom index in the connected residue of the atom that's making a connection to THIS residue's connection #1.  (Convoluted, I know.)
		}
	}

	///////////////////////////////////////////
	// second atom -- for sure in seqpos
	id2.rsd()    = seqpos;
	id2.atomno() = mainchain[ torsion ];

	///////////////////////////////////////////
	// third and fourth atoms, may be in seqpos+1
	if ( torsion+2 <= ntorsions ) { //In this case, all of the remaining atoms are within the residue (e.g. theta in beta-amino acids, phi in alpha-amino acids).
		id3.rsd()    = seqpos;
		id3.atomno() = mainchain[ torsion+1 ];
		id4.rsd()    = seqpos;
		id4.atomno() = mainchain[ torsion+2 ];

	} else {
		if ( rsd.has_variant_type( chemical::CTERM_CONNECT ) && residues_[ chain_begin( rsd.chain() ) ]->has_variant_type( chemical::NTERM_CONNECT ) ) {
			// assume that they are connected as a cyclic polymer. psi and omg are defined eventhough it is last residue
			AtomIndices const & cyclic_partner_mainchain( residue_( chain_begin( rsd.chain() ) ).mainchain_atoms() );
			if ( torsion == 2 ) {
				//TR << "HI MOM, CTERM, 2" << std::endl;
				id3.rsd() = seqpos;                     id3.atomno() = mainchain[3];
				id4.rsd() = chain_begin( rsd.chain() ); id4.atomno() = cyclic_partner_mainchain[1];
			} else if ( torsion == 3 ) {
				//TR << "HI MOM, CTERM, 3" << std::endl;
				id3.rsd() = chain_begin( rsd.chain() ); id3.atomno() = cyclic_partner_mainchain[1];
				id4.rsd() = chain_begin( rsd.chain() ); id4.atomno() = cyclic_partner_mainchain[2];
			}

	} else if ( fold_tree_->is_cutpoint( seqpos ) &&
				!( seqpos==residues_.size() && rsd.has_upper_connect() && !rsd.connection_incomplete( rsd.type().upper_connect_id() ) /*special case -- last residue is connected to something at its upper connection*/)
		) {
			if ( rsd.has_variant_type( chemical::CUTPOINT_LOWER ) ) {
				if ( torsion+1 == ntorsions ) {
					id3.rsd()    = seqpos;
					id3.atomno() = mainchain[ torsion+1 ];
					id4.rsd()    = seqpos;
					id4.atomno() = rsd.atom_index( "OVL1" );
				} else {
					assert( torsion == ntorsions );
					id3.rsd()    = seqpos;
					id3.atomno() = rsd.atom_index( "OVL1" );
					id4.rsd()    = seqpos;
					id4.atomno() = rsd.atom_index( "OVL2" );
				}
			} else if ( rsd.has_variant_type( "C_METHYLAMIDATION" ) ) {
				//ugly.
				id3.rsd() = seqpos;
				id4.rsd() = seqpos;
				if ( torsion == 2 ) /*psi*/ {
					id3.atomno() = rsd.atom_index( " C  ");
					id4.atomno() = rsd.atom_index( " NR ");
				} else {
					assert( torsion == 3 );
					id3.atomno() = rsd.atom_index( " NR ");
					id4.atomno() = rsd.atom_index( " CS ");
				}
			} else if ( rsd.has_variant_type( "THREE_PRIME_PHOSPHATE" ) ) {
				//ugly again. -- rhiju.
				id3.rsd() = seqpos;
				id4.rsd() = seqpos;
				if ( torsion+1 == ntorsions ) /*epsilon*/ {
					id3.atomno() = mainchain[ torsion+1 ];
					id4.atomno() = rsd.atom_index( "YP  ");
				} else {
					assert( torsion == ntorsions ); /*zeta*/
					id3.atomno() = rsd.atom_index( "YP  ");
					id4.atomno() = rsd.atom_index( "YO5'");
				}
			} else {
				// last two bb-torsions not well-defined
				return true; // FAILURE
			}
		} else if ( rsd.has_variant_type( chemical::METHYLATED_CTERMINUS ) ) {
			if ( torsion+1 == ntorsions ) {
				id3.rsd() = seqpos;
				id3.atomno() = mainchain[ torsion+1 ];
				id4.rsd() = seqpos;
				id4.atomno() = rsd.atom_index( "NM" );
			} else {
				assert( torsion == ntorsions );
				id3.rsd()    = seqpos;
				id3.atomno() = rsd.atom_index( "NM" );
				id4.rsd()    = seqpos;
				id4.atomno() = rsd.atom_index( "CN" );
			}
		} else	{ //If this is NOT a cutpoint
			if(!rsd.has_upper_connect() || rsd.connection_incomplete( rsd.type().upper_connect_id() )) return true; //FAIL if this residue is not connected to anything.
			AtomIndices const & next_mainchain ( residue_( rsd.residue_connection_partner( rsd.type().upper_connect_id() ) ).mainchain_atoms() );

			if ( torsion+1 == ntorsions ) { //If this is the second-to-last torsion angle (e.g. psi, in alpha-amino acids)
				id3.rsd()    = seqpos;
				id3.atomno() = mainchain[ torsion+1 ];
				//id4.rsd()    = seqpos+1;
				//id4.atomno() = next_mainchain[ 1 ];

				//Altered by VKM, 12 June 2014: we want to fish out whatever atom the residue is connected to.
				id4.rsd() = rsd.residue_connection_partner(rsd.type().upper_connect_id()); //Get the residue index of the residue connected to this residue at this residue's connection #2.
				if(residue_(id4.rsd()).connect_map_size() < rsd.residue_connection_conn_id( rsd.type().upper_connect_id() )) {
					return true; //FAIL if the residue connected at upper is connected improperly.
				}
				id4.atomno() =  residue_(id4.rsd()).residue_connect_atom_index( rsd.residue_connection_conn_id(rsd.type().upper_connect_id()) ); //Get the atom index in the connected residue of the atom that's making a connection to THIS residue's connection #2.

			} else { //If this is the last torsion angle (e.g. omega, in alpha- or beta-amino acids).
				assert( torsion == ntorsions );
				//id3.rsd()    = seqpos+1;
				//id3.atomno() = next_mainchain[ 1 ];

				//Altered by VKM, 12 June 2014: we want to fish out whatever atom the residue is connected to.
				id3.rsd() = rsd.residue_connection_partner(rsd.type().upper_connect_id()); //Get the residue index of the residue connected to this residue at this residue's connection #2.
				if(residue_(id3.rsd()).connect_map_size() < rsd.residue_connection_conn_id( rsd.type().upper_connect_id() )) {
					return true; //FAIL if the residue connected at upper is connected improperly.
				}
				id3.atomno() =  residue_(id3.rsd()).residue_connect_atom_index( rsd.residue_connection_conn_id(rsd.type().upper_connect_id()) ); //Get the atom index in the connected residue of the atom that's making a connection to THIS residue's connection #2.

				if ( next_mainchain.size() >= 2 ) {
					id4.rsd()    = id3.rsd();
					if(id3.atomno() == next_mainchain[1]) {
						id4.atomno() = next_mainchain[ 2 ];
					} else {
						id4.atomno() = residue_(id3.rsd()).type().icoor(id3.atomno()).stub_atom1().atomno(); //Let the fourth atom index be the parent atom of the third if it is not the first mainchain atom.
					}
				} else {
					// tricky... a single-mainchain-atom polymer residue.
					// TODO -- remove the seqpos + 1 assumption here.
					if ( fold_tree_->is_cutpoint( seqpos+1 ) ) {
						if ( residue_( seqpos+1 ).has_variant_type( chemical::CUTPOINT_LOWER ) ) {
							id4.rsd()    = seqpos+1;
							id4.atomno() = residue_( seqpos+1 ).atom_index( "OVL1" );
						} else {
							return true; // failure
						}
					} else {
						id4.rsd() = seqpos+2;
						id4.atomno() = residue_( seqpos+2 ).mainchain_atom(1);
					}
				} // next_mainchain has size>=2 ?
			} // torsion+1 == ntorsions?
		} // seqpos is a cutpoint
	} // torsion+2 <= ntorsions

	}

	//std::cout << "2. checking backbone torsion angle atoms for torsionID " << torsion << " in residue " << //DELETE ME -- for debugging only.
	//		 seqpos << " and atoms " << rsd.atom_name(id1.atomno()) << " " << rsd.atom_name(id2.atomno()) //DELETE ME -- for debugging only.
	//		 << " " << rsd.atom_name(id3.atomno()) << " " << rsd.atom_name(id4.atomno()) << std::endl; //DELETE ME -- for debugging only.

	fail = false;

	return fail;
}

/// @brief Helper function to determine whether two atoms have a chemical bond linking them.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool Conformation::atoms_are_bonded(
	AtomID const &id1,
	AtomID const &id2
) const {
	if(id1.rsd()==id2.rsd()) { //If the atoms are in the same residue
		if(residue_(id1.rsd()).type().atoms_are_bonded( id1.atomno(), id2.atomno() ) ) return true;
		else return false; //Should be redundant.
	} else { //If the atoms are in different residues
		if(!residue_(id1.rsd()).is_bonded( residue_(id2.rsd()) ) ) return false;
		//Get the list of connection ids in residue 1 that connect to residue 2:
		utility::vector1< core::Size > connlist = residue_(id1.rsd()).connections_to_residue(id2.rsd());

		//Loop through these, and check each for a connection between the atom pair in question. 
		for(core::Size i=1, imax=connlist.size(); i<=imax; ++i) {
			//Does the current connection id in residue 1 correspond to the correct atom index in residue 1?:
			if(residue_(id1.rsd()).residue_connect_atom_index(connlist[i]) != id1.atomno() ) continue;
			//Redundant check: the current connection id in residue 1 should connect to residue 2:
			assert(residue_(id1.rsd()).connected_residue_at_resconn(connlist[i]) == id2.rsd() );
			//Does the current connection id in residue 1 connect to a residue id in residue 2 with the correct atom id?:
			if(   residue_(id2.rsd()).residue_connect_atom_index( residue_(id1.rsd()).residue_connection_conn_id( connlist[i]) ) == id2.atomno() ) return true;
		} 
	}
	return false;
}


// setting the moved data ////////////////////////////////////////////////////////////////////////////////////////////

/// @details The AtomTree is responsible for tracking the set of residues whose
/// coordinates need updating, and informs the Conformation object of these residues.
/// The Conformation only updates coordinates for this subset of residues.
/// @note Always safe to call. Nothing will happen unless coords_need_updating_ is true.
void
Conformation::update_residue_coordinates() const
{
	if ( !residue_coordinates_need_updating_ ) return;
	residue_coordinates_need_updating_ = false;

	// fill in the new xyz values in the pose's Residues
	// this could be made way faster...
	//
	// options: give the atomtree atoms links to the corresponding
	// atoms in the Residue
	//

	PROF_START( basic::UPDATE_RESIDUE_COORDINATES );
	for ( kinematics::ResidueListIterator iter = atom_tree_->residue_xyz_change_list_begin(),
		iter_end = atom_tree_->residue_xyz_change_list_end(); iter != iter_end; ++iter ) {
		update_residue_coordinates( *iter, false );
	}
	atom_tree_->note_coordinate_change_registered();
	PROF_STOP( basic::UPDATE_RESIDUE_COORDINATES );

	notify_xyz_obs( XYZEvent( this ) );
}


void
Conformation::update_residue_coordinates( Size const seqpos, bool const fire_signal ) const
{
	Residue & rsd( *residues_[ seqpos ] );
	for ( Size j=1, j_end = rsd.natoms(); j<= j_end; ++j ) {
		rsd.set_xyz( j, atom_tree_->xyz( AtomID(j,seqpos) ) );
	}
	rsd.update_actcoord();

	//update orbital coords!
	this->update_orbital_coords(rsd);

	if ( fire_signal ) {
		notify_xyz_obs( XYZEvent( this ) );
	}
}


void
Conformation::rederive_chain_ids()
{
	Size chain(1);
	for ( Size i=1; i<= size(); ++i ) {
		residues_[i]->chain( chain );
		if ( chain <= chain_endings_.size() && (Size) chain_endings_[chain] == i ) {
			++chain;
		}
	}

}

void
Conformation::rederive_chain_endings()
{
	chain_endings_.clear();
	for ( Size i=1, ie=size()-1; i<=ie; ++i ) {
		if ( residues_[i+1]->chain() != residues_[i]->chain() ) {
			assert( residues_[i+1]->chain() > residues_[i]->chain() );
			// 			assert( residues_[i+1]->chain() == residues_[i]->chain() + 1 );
			chain_endings_.push_back( i );
		}
	}
}


/// @note Always safe to call. Nothing will happen unless torsions_need_updating_ is true.
void
Conformation::update_residue_torsions() const
{
	if ( !residue_torsions_need_updating_ ) return;
	residue_torsions_need_updating_ = false;

	// fill in the new torsion angles in the pose's Residues
	// this could be made way faster...
	//
	// eg by keeping track of the subset of torsions that could have changed
	PROF_START( basic::UPDATE_RESIDUE_TORSIONS );
	for ( Size i=1, i_end = size(); i<= i_end; ++i ) {
		update_residue_torsions( i, false );
		residues_[ i ]->update_actcoord();
	}
	PROF_STOP( basic::UPDATE_RESIDUE_TORSIONS );

	notify_xyz_obs( XYZEvent( this ) );
}

void
Conformation::update_residue_torsions( Size const seqpos, bool const fire_signal ) const
{
	using id::BB;
	using id::CHI;

	Residue & rsd( *residues_[seqpos] );

	// mainchain
	for ( Size j=1, j_end = rsd.mainchain_torsions().size(); j<= j_end; ++j ) {
		// these hacks are for cyclization
		// phi of first residue in cycle
		if ( rsd.has_variant_type( chemical::NTERM_CONNECT ) && residues_[ chain_end( rsd.chain() ) ]->has_variant_type( chemical::CTERM_CONNECT ) && j == 1 ) {
			// calc torsion from xyz since those DOFs are not in atom tree
			// not sure how well this will work since the xyz might not be up to date

			TorsionID tid(seqpos,BB,j);
			AtomID id1, id2, id3, id4;
			get_torsion_angle_atom_ids( tid, id1, id2, id3, id4 );

			rsd.mainchain_torsions()[ j ] = numeric::dihedral_degrees(
				residues_[ id1.rsd() ]->xyz( id1.atomno() ),
				residues_[ id2.rsd() ]->xyz( id2.atomno() ),
				residues_[ id3.rsd() ]->xyz( id3.atomno() ),
				residues_[ id4.rsd() ]->xyz( id4.atomno() ) );
		}
		// psi and omg of last residue in cycle
		else if ( rsd.has_variant_type( chemical::CTERM_CONNECT ) && residues_[ chain_begin( rsd.chain() ) ]->has_variant_type( chemical::NTERM_CONNECT ) && ( j == 2 || j == 3 ) ) {
			// calc torsion from xyz since those DOFs are not in atom tree
			// not sure how well this will work since the xyz might not be up to date
			TorsionID tid(seqpos,BB,j);
			AtomID id1, id2, id3, id4;
			get_torsion_angle_atom_ids( tid, id1, id2, id3, id4 );

			rsd.mainchain_torsions()[ j ] = numeric::dihedral_degrees(
				residues_[ id1.rsd() ]->xyz( id1.atomno() ),
				residues_[ id2.rsd() ]->xyz( id2.atomno() ),
				residues_[ id3.rsd() ]->xyz( id3.atomno() ),
				residues_[ id4.rsd() ]->xyz( id4.atomno() ) );
		}
		// normal residue
		else {
			rsd.mainchain_torsions()[ j ] = atom_tree_torsion( TorsionID(seqpos,BB,j) );
		}
	}

	// chi
	for ( Size j=1, j_end = rsd.nchi(); j<= j_end; ++j ) {
		rsd.chi()[ j ] = atom_tree_torsion( TorsionID( seqpos, CHI, j ) );
	}

	//update orbital coords!
	this->update_orbital_coords(rsd);

	if ( fire_signal ) {
		notify_xyz_obs( XYZEvent( this ) );
	}
}


void
Conformation::add_pseudobond(
	Size lr,
	Size lr_connid,
	Size ur,
	Size ur_connid,
	Size nbonds
)
{
	PseudoBondCollectionOP new_pbs;
	assert(residues_.size()>0);
	assert(lr>0);
	assert(lr_connid>0);//?
	assert(ur>0);
	assert(ur_connid>0);//?

	if ( residues_[ lr ]->is_pseudo_bonded( *residues_[ ur ] ) ) {
		PseudoBondCollectionCOP existing_pbs = residues_[ lr ]->get_pseudobonds_to_residue( ur );
		new_pbs = new PseudoBondCollection( *existing_pbs );
	} else {
		new_pbs = new PseudoBondCollection();
	}
	PseudoBond new_pb;
	new_pb.lr( lr ); new_pb.lr_conn_id( lr_connid );
	new_pb.ur( ur ); new_pb.ur_conn_id( ur_connid );
	new_pb.nbonds( 2 + nbonds );
	new_pbs->push_back( new_pb );
	residues_[ lr ]->set_pseudobonds_to_residue( ur, new_pbs );
	residues_[ ur ]->set_pseudobonds_to_residue( lr, new_pbs );
}


void
Conformation::in_place_copy(
	Conformation const & src
)
{
	assert( src.size() == size() );

	/// BEGIN IN PLACE OPTIMIZATION

	// Do not allocate new residue objects, just reuse the old ones
	for ( Size ii = 1; ii <= size(); ++ii ) {
		runtime_assert( & residues_[ ii ]->type() ==  & src.residues_[ ii ]->type() );
		runtime_assert( residues_[ ii ]->seqpos() == src.residues_[ ii ]->seqpos() );
		runtime_assert( residues_[ ii ]->chain() == src.residues_[ ii ]->chain() );
		runtime_assert( residues_[ ii ]->connections_match( *src.residues_[ ii ] ));
		for ( Size jj = 1; jj <= residues_[ ii ]->natoms(); ++jj ) {
			residues_[ ii ]->set_xyz( jj, src.residues_[ ii ]->xyz( jj ) );
		}
		residues_[ ii ]->mainchain_torsions( src.residues_[ ii ]->mainchain_torsions() );
		residues_[ ii ]->chi( src.residues_[ ii ]->chi() );
		residues_[ ii ]->actcoord() = src.residues_[ ii ]->actcoord();
	}

	/// END IN PLACE OPTIMIZATION

	// kinematics
	(*fold_tree_) = (*src.fold_tree_);
	(*atom_tree_) = (*src.atom_tree_); // internally performs an in-place copy if it can

	// chain info
	chain_endings_ = src.chain_endings_;
	// secstruct
	secstruct_ = src.secstruct_;

	// bookkeeping
	residue_coordinates_need_updating_ = src.residue_coordinates_need_updating_;
	residue_torsions_need_updating_    = src.residue_torsions_need_updating_;

	dof_moved_ = src.dof_moved_;
	xyz_moved_ = src.xyz_moved_;

	structure_moved_ = src.structure_moved_;
}


/// @details Copy any un-registered coordinate or DOF changes into the existing residues.
/// For now, the AtomTree only tracks which Residues need external coordinate changes,
/// and not internal coordinate changes.  When internal coordinates go out-of-date in
/// the atom tree, the "update_internal_coordinates" recursion begins at the root.  Moreover,
/// the Conformation updates internal coordinates for all residues when it updates. However,
/// when this changes and the AtomTree starts tracking which residues need to have their
/// internal coordinates updated, then it will be necessary to add a call to
/// update_residue_torsions() here.
void
Conformation::pre_nresidue_change()
{
	update_residue_coordinates();
}


/////////////////////////////////////////////////////
/** work in progress
void
insert_chemical_chainbreak_between_polymer_residues( Size const lower_seqpos )
{
	ResidueOP lower_rsd( residues_[ lower_seqpos ] );

	// the type of the desired variant residue
	ResidueTypeSet const & rsd_set( lower_rsd->residue_type_set() );

	ResidueType const & new_lower_rsd_type
		( rsd_set.get_residue_type_with_variant_removed( lower_rsd.type(), UPPER_TERMINUS ) );
	ResidueOP new_lower_rsd( ResidueFactory::create_residue( new_lower_rsd_type ) );
}
**/
/////////////////////////////////////////////////////


// observer notifications ////////////////////////////////////////////////////////////////////////////////////////////

// notify ConnectionEvent observers
/// @remarks called upon a change in the state of connection between
///  the Conformation and the observer (e.g. destruction of Conformation
///  or transfer of connection)
void
Conformation::notify_connection_obs( ConnectionEvent const & e ) const {
	connection_obs_hub_( e );
}


// notify GeneralEvent observers
/// @remarks should only be called when there are no other suitable event types
///  since specific event notifications will automatically fire a GeneralEvent signal
void
Conformation::notify_general_obs( GeneralEvent const & e ) const {
	general_obs_hub_( e );
}


// notify IdentityEvent observers
/// @param e the event
/// @param fire_general fire a GeneralEvent afterwards? default true
void
Conformation::notify_identity_obs( IdentityEvent const & e, bool const fire_general ) const {
	identity_obs_hub_( e );
	if ( fire_general ) {
		notify_general_obs( e );
	}
}


// notify LengthEvent observers
/// @param e the event
/// @param fire_general fire a GeneralEvent afterwards? default true
void
Conformation::notify_length_obs( LengthEvent const & e, bool const fire_general ) const {
	length_obs_hub_( e );
	if ( fire_general ) {
		notify_general_obs( e );
	}
}


// notify XYZEvent observers
/// @param e the event
/// @param fire_general fire a GeneralEvent afterwards? default true
void
Conformation::notify_xyz_obs( XYZEvent const & e, bool const fire_general ) const {
	xyz_obs_hub_( e );
	if ( fire_general ) {
		notify_general_obs( e );
	}
}


std::ostream &operator<< (std::ostream &os, Conformation const &conf)
{
 	conf.show_residue_connections(os);
    return os;
}

} // namespace kinematics
} // namespace core
