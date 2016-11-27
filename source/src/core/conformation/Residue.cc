// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/conformation/Residue.cc
/// @brief  Method definitions for the Residue class
/// @author Phil Bradley

// Unit header
#include <core/conformation/Residue.hh>

// Package headers
#include <core/conformation/PseudoBond.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/residue_datacache.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/util.hh>

// Project headers
#include <core/kinematics/Stub.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/rings/RingConformer.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomICoor.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// Basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// Utility Headers
#include <utility/assert.hh>
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <iostream>

#ifdef    SERIALIZATION
// Package serialization headers
#include <core/chemical/ResidueType.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace conformation {

static THREAD_LOCAL basic::Tracer TR( "core.conformation.Residue" );

/// @brief This function enforces the fact that ResidueTypes must
/// be constructed with non-null-pointer ResidueTypeCOPs.
/// @details This must be a function as the initialization of the
/// Residue's rsd_type_ reference must occur in the constructor.
chemical::ResidueType const &
reference_from_restype_ptr( chemical::ResidueTypeCOP rsd_type )
{
	debug_assert( rsd_type );
	return *rsd_type;
}

/// @details Constructor from ResidueTypeCOP; sets coords to ideal values
/// create a residue of type residue_type_in.
/// @note Dummmy arg to prevent secret type conversions from ResidueTypeCOP to Residue
Residue::Residue( ResidueTypeCOP rsd_type_in, bool const /*dummy_arg*/ ):
	utility::pointer::ReferenceCount(),
	rsd_type_ptr_( rsd_type_in ),
	rsd_type_( reference_from_restype_ptr( rsd_type_in )),
	seqpos_( 0 ),
	mirrored_relative_to_type_(false),
	chain_( 0 ),
	chi_( rsd_type_.nchi(), 0.0 ), // uninit
	nus_( rsd_type_.n_nus(), 0.0 ),
	mainchain_torsions_( rsd_type_.mainchain_atoms().size(), 0.0 ),
	actcoord_( 0.0 ),
	data_cache_( nullptr ),
	nonstandard_polymer_( false ),
	connect_map_( rsd_type_.n_possible_residue_connections() )
{
	// Assign atoms.
	for ( Size i=1; i<= rsd_type_.natoms(); ++i ) {
		atoms_.push_back( Atom( rsd_type_.atom(i).ideal_xyz(), rsd_type_.atom(i).atom_type_index(),
			rsd_type_.atom(i).mm_atom_type_index() ) );
	}

	update_nus();
	assign_orbitals();

}

/// @details Constructor from residue type; sets coords to ideal values
/// create a residue of type residue_type_in.
/// @note Dummmy arg to prevent secret type conversions from ResidueType to Residue
Residue::Residue( ResidueType const & rsd_type_in, bool const /*dummy_arg*/ ):
	utility::pointer::ReferenceCount(),
	rsd_type_ptr_( rsd_type_in.get_self_ptr() ),
	rsd_type_( rsd_type_in ),
	seqpos_( 0 ),
	mirrored_relative_to_type_(false),
	chain_( 0 ),
	chi_( rsd_type_.nchi(), 0.0 ), // uninit
	nus_( rsd_type_.n_nus(), 0.0 ),
	mainchain_torsions_( rsd_type_.mainchain_atoms().size(), 0.0 ),
	actcoord_( 0.0 ),
	data_cache_( nullptr ),
	nonstandard_polymer_( false ),
	connect_map_( rsd_type_.n_possible_residue_connections() )
{
	// Assign atoms.
	for ( Size i=1; i<= rsd_type_.natoms(); ++i ) {
		atoms_.push_back( Atom( rsd_type_.atom(i).ideal_xyz(), rsd_type_.atom(i).atom_type_index(),
			rsd_type_.atom(i).mm_atom_type_index() ) );
	}

	update_nus();
	assign_orbitals();

}

/// @details Create a residue/rotamer of type rsd_type_in placed at the position occupied by current_rsd
/// Used primarily in rotamer building. The newly created Residue has the same sequence position, chain id
/// and mainchain torsion angles as current_rsd. It has a ResidueType as defined by rsd_type_in. Its side-chain
/// chi angles are uninitialized as all 0.0 and sidechain atom coords are from ideal coords. Its backbone is aligned
/// with that of current_rsd.
/// Its residue connections and its pseudobonds must be initialized from the original residue.
Residue::Residue(
	ResidueType const & rsd_type_in,
	Residue const & current_rsd,
	Conformation const & conformation,
	bool preserve_c_beta
):
	utility::pointer::ReferenceCount(),
	rsd_type_ptr_( rsd_type_in.get_self_ptr() ),
	rsd_type_( rsd_type_in ),
	seqpos_( current_rsd.seqpos() ),
	mirrored_relative_to_type_(current_rsd.mirrored_relative_to_type()),
	chain_( current_rsd.chain() ),
	chi_( rsd_type_.nchi(), 0.0 ), // uninit
	nus_( current_rsd.nus() ),
	mainchain_torsions_( current_rsd.mainchain_torsions() ),
	actcoord_( 0.0 ),
	data_cache_( nullptr ),
	nonstandard_polymer_( current_rsd.nonstandard_polymer_ ),
	connect_map_( current_rsd.connect_map_ ),
	connections_to_residues_( current_rsd.connections_to_residues_ ),
	pseudobonds_( current_rsd.pseudobonds_ )
{
	// Assign atoms.
	for ( Size i=1; i<= rsd_type_.natoms(); ++i ) {
		atoms_.push_back( Atom( rsd_type_.atom(i).ideal_xyz(), rsd_type_.atom(i).atom_type_index(),
			rsd_type_.atom(i).mm_atom_type_index() ));
	}

	debug_assert( current_rsd.mainchain_torsions().size() == rsd_type_.mainchain_atoms().size() );

	// Now orient the residue.
	place( current_rsd, conformation, preserve_c_beta );

	// Assumption: if two residue types have the same number of residue connections,
	// then their residue connections are "the same" residue connections.
	// This assumption works perfectly for all amino acids, except CYD.
	//
	// THIS REALLY NEEDS TO BE FIXED AT SOME POINT

	if ( rsd_type_in.n_possible_residue_connections() != current_rsd.type().n_possible_residue_connections() ) {
		if ( ! current_rsd.pseudobonds_.empty() ) {
			TR.Error << "Unable to handle change in the number of residue connections in the presence of pseudobonds!";
			TR.Error << std::endl;
			utility_exit();
		}

		copy_residue_connections( current_rsd );
	}

	// This seems a little silly, but the update of chis doesn't seem to occur automatically in
	// any of the functions above.
	for ( Size chino = 1; chino <= rsd_type_.nchi(); ++chino ) {
		AtomIndices const & chi_atoms( rsd_type_.chi_atoms( chino ) );

		// get the current chi angle
		Real const current_chi(numeric::dihedral_degrees(
			atom( chi_atoms[1] ).xyz(),
			atom( chi_atoms[2] ).xyz(),
			atom( chi_atoms[3] ).xyz(),
			atom( chi_atoms[4] ).xyz()));
		chi_[ chino ] = current_chi;
	}

	update_nus();
	assign_orbitals();

}

/// @brief Copy constructor.
///
Residue::Residue( Residue const & src ) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< Residue >(),
	rsd_type_ptr_( src.rsd_type_ptr_ ),
	rsd_type_(src.rsd_type_)
{
	init_residue_from_other( src );
}

/// @brief Copy constructor that preserves everything EXCEPT flips ResidueType to its mirror counterpart.
/// @details The flip_chirality option governs whether this residue keeps the old type or gets the mirror-image type.  All other data
/// (including xyz coordinates) are preserved, so this ONLY switches the type, not the geometry.
/// @author Vikram K. Mulligan.
Residue::Residue( Residue const & src, bool const flip_chirality ):
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< Residue >(),
	rsd_type_ptr_( flip_chirality ? src.residue_type_set()->get_mirrored_type( src.rsd_type_ptr_ ) : src.rsd_type_ptr_ ),
	rsd_type_( flip_chirality ? * (src.residue_type_set()->get_mirrored_type( src.rsd_type_.get_self_ptr() ) ) : src.rsd_type_ ) /*I don't assume that src.rsd_type_ and src.rsd_type_ptr_ point to the same type, though they should.*/
{
	init_residue_from_other( src );
	if ( type().is_achiral_backbone() ) set_mirrored_relative_to_type( !src.mirrored_relative_to_type() ); //For achiral residues, we need to record whether the geometry is mirrored relative to the params file.  (1H and 2H may be chemically indistinguishable, but Rosetta distinguishes them.)
}

Residue::~Residue() = default;

/// @brief Function called by both copy constructors, to avoid code duplication.
/// @details As private member variables are added, add them to this to copy them.
/// @author Vikram K. Mulligan.
void
Residue::init_residue_from_other(
	Residue const &src
) {
	atoms_ = src.atoms_;
	orbitals_ = src.orbitals_;
	seqpos_ = src.seqpos_;
	mirrored_relative_to_type_ = src.mirrored_relative_to_type_;
	chain_ = src.chain_;
	chi_ = src.chi_;
	nus_ = src.nus_;
	mainchain_torsions_ = src.mainchain_torsions_;
	actcoord_ = src.actcoord_;
	if ( src.data_cache_ != nullptr ) {
		if ( data_cache_ != nullptr ) ( *data_cache_) = (*src.data_cache_);
		else data_cache_ = basic::datacache::BasicDataCacheOP( new basic::datacache::BasicDataCache( *src.data_cache_) );
	}
	nonstandard_polymer_ = src.nonstandard_polymer_;
	connect_map_ = src.connect_map_;
	connections_to_residues_ = src.connections_to_residues_;
	pseudobonds_ = src.pseudobonds_;
}


/// @details make a copy of this residue( allocate actual memory for it )
ResidueOP
Residue::clone() const
{
	return ResidueOP( new Residue( *this ) );
}

/// @brief Copy this residue( allocate actual memory for it ), keeping everything the same EXCEPT the type.
/// @details Switches the ResidueType to the mirror type (D->L or L->D).  Preserves it for achiral residues.
/// @note This function is the best way to convert a D-residue to its L-counterpart, or an L-residue to its D-counterpart.
/// It assumes that you've already mirrored all of the coordinates, and just allows you to generate a replacement residue
/// of the mirror type that preserves all other Residue information (connections, seqpos, xyz coordinates of all atoms,
/// variant types, etc.).
/// @author Vikram K. Mulligan (vmullig@uw.edu)
ResidueOP
Residue::clone_flipping_chirality() const
{
	return ResidueOP( new Residue( *this, true /*flip the chirality*/ ) );
}



/// @author Labonte <JWLabonte@jhu.edu>
void
Residue::show( std::ostream & output, bool output_atomic_details ) const
{
	using namespace std;
	using namespace chemical::rings;

	output << "Residue " << seqpos_ << ": ";
	rsd_type_.show( output, output_atomic_details );
	if ( rsd_type_.is_cyclic() ) {
		Size const n_rings( rsd_type_.n_rings() );
		for ( core::uint i( 1 ); i <= n_rings; ++i ) {
			output << "Ring Conformer: " << ring_conformer( i ) << endl;

			AtomIndices const ring_atoms( rsd_type_.ring_atoms( i ) );
			Size const n_ring_atoms( ring_atoms.size() );
			for ( core::uint j( 1 ); j <= n_ring_atoms; ++j ) {
				AtomIndices const potential_substituents( get_adjacent_heavy_atoms( ring_atoms[ j ] ) );
				Size const n_potential_substituents( potential_substituents.size() );
				for ( core::uint k( 1 ); k <= n_potential_substituents; ++k ) {
					if ( ( ! ring_atoms.contains( potential_substituents[ k ] ) ) &&
							( ! is_virtual( potential_substituents[ k ] ) ) ) {
						// This atom must be exocyclic.
						// Is it axial or equatorial?
						output << ' ' << atom_name( potential_substituents[ k ] ) << ": ";
						switch ( is_atom_axial_or_equatorial_to_ring( *this, potential_substituents[ k ], ring_atoms ) ) {
						case AXIAL :
							output << "axial" << endl;
							break;
						case EQUATORIAL :
							output << "equatorial" << endl;
							break;
						case NEITHER :
							output << "neither axial nor equatorial" << endl;
							break;
						}
						break;
					}
				}
			}
		}
	}
	output << "Atom Coordinates:" << endl;
	Size n_atoms = natoms();
	for ( core::uint i = 1; i <= n_atoms; ++i ) {
		conformation::Atom const & atom_coords( atoms_[ i ] );
		output << "  " << atom_name( i ) << ": ";
		atom_coords.show( output );
		if ( is_virtual( i ) ) {
			output << " (virtual)";
		}
		output << endl;
	}
	output << "Mirrored relative to coordinates in ResidueType: " << (mirrored_relative_to_type() ? "TRUE" : "FALSE") << std::endl;
}


Size
Residue::atom_type_index( Size const atomno ) const
{
	return atoms_[ atomno ].type();
}

Real
Residue::atomic_charge( int const atomno ) const
{
	return rsd_type_.atom( atomno ).charge();
}

Vector const &
Residue::xyz( Size const atm_index ) const
{
	return atoms_[ atm_index ].xyz();
}

Vector const &
Residue::xyz( std::string const & atm_name ) const
{
	return atom( atm_name ).xyz();
}

void
Residue::set_xyz( core::Size const atm_index, Vector const & xyz_in )
{
	atoms_[ atm_index ].xyz( xyz_in );
}

void
Residue::set_xyz( std::string const & atm_name, Vector const & xyz_in )
{
	atom( atm_name ).xyz( xyz_in );

}


// Return the CarbohydrateInfo object containing sugar-specific properties for this residue.
core::chemical::carbohydrates::CarbohydrateInfoCOP
Residue::carbohydrate_info() const
{
	debug_assert(rsd_type_.is_carbohydrate());
	PyAssert(rsd_type_.is_carbohydrate(), "Residue::carbohydrate_info(): This residue is not a carbohydrate!");

	return rsd_type_.carbohydrate_info();
}


// Return the current RingConformer of this residue's nth ring.
chemical::rings::RingConformer const &
Residue::ring_conformer( core::uint const ring_num, core::Real limit/*=90*/ ) const
{
	PyAssert( ( ring_num <= rsd_type_.n_rings() ),
		"Residue::ring_conformer(core::uint const ring_num): variable ring_num is out of range!" );
	debug_assert( rsd_type_.is_cyclic() );

	// First, figure out which nus belong to this ring.
	Size n_nus_on_previous_rings( 0 );
	for ( uint previous_ring_num( ring_num -1 ); previous_ring_num > 0; --previous_ring_num ) {
		// ( Each ring has one fewer nus associated with it than the size of the ring. )
		n_nus_on_previous_rings += rsd_type_.ring_atoms( previous_ring_num ).size() - 1;
	}

	Size const n_nus( rsd_type_.ring_atoms( ring_num ).size() - 1 );
	utility::vector1< Angle > nus;
	for ( uint i( 1 ); i <= n_nus; ++i ) {
		nus.push_back( nus_[ n_nus_on_previous_rings + i ] );
	}

	return rsd_type_.ring_conformer_set( ring_num )->get_ideal_conformer_from_nus( nus_, limit );
}

core::Size
Residue::n_current_residue_connections() const
{
	//for connection in range( 1, rsd.n_residue_connections() + 1 ):
	//    if rsd.connected_residue_at_resconn( connection ) != 0:
	//        i+=1
	//return i
	core::Size connections = 0;
	for ( core::Size connect_id = 1; connect_id <= n_possible_residue_connections(); ++connect_id ) {
		if ( connected_residue_at_resconn( connect_id ) != 0 ) {
			connections+=1;
		}
	}
	return connections;

}

bool Residue::connections_match( Residue const & other ) const
{
	if ( connect_map_.size() != other.connect_map_.size() ) return false;
	//if ( connections_to_residues_.size() != other.connections_to_residues_.size()) return false; // duplicate data
	if ( pseudobonds_.size() != other.pseudobonds_.size() ) return false;

	for ( Size ii = 1; ii <= connect_map_.size(); ++ii ) {
		if ( connect_map_[ ii ] != other.connect_map_[ ii ] ) return false;
	}
	for ( auto
			iter = pseudobonds_.begin(), iter_end = pseudobonds_.end(),
			other_iter_end = other.pseudobonds_.end();
			iter != iter_end; ++iter ) {
		auto other_iter = other.pseudobonds_.find( iter->first );
		if ( other_iter == other_iter_end ) return false;
		if ( iter->second != other_iter->second ) return false; // pointer comparison
		//if ( ! (*(iter->second) == *(other_iter->second) ) ) return false;
	}
	return true;
}

bool
Residue::is_similar_rotamer( Residue const & other ) const
{
	utility::vector1< Real > this_chi = chi_;
	utility::vector1< Real > other_chi = other.chi();
	bool match = true;
	if ( chi_.size() != other_chi.size() || rsd_type_.aa() != other.aa() || rsd_type_.name3() != other.name3() ) {
		return false;
	} else {
		for ( Size i = 1; i<= chi_.size(); ++i ) {
			if ( std::abs( this_chi[i] - other_chi[i]) >= 5 ) {
				match = false;
			}
		}
	}
	return match;
}

/// @brief Search through the other residue for connections to this residue, and
/// ensure that this residue's connect_map is up to date with that residue's
/// connection indices (and residue number).
/// @details Throws an error if this residue doesn't have a connection id indicated
/// by the other residue.  Overwrites connection ids otherwise, with no consideration
/// of whether the original connection was to other_rsd.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
Residue::update_connections_to_other_residue( Residue const &other_rsd)
{
	for ( core::Size ic=1, ic_max=other_rsd.type().n_possible_residue_connections(); ic<=ic_max; ++ic ) {
		if ( other_rsd.connected_residue_at_resconn(ic)==seqpos() ) { //If the other's connection lists this residue's sequence position
			core::Size const this_conn_id = other_rsd.connect_map(ic).connid();
			//TR << "this_res=" << seqpos() << "other_res=" << other_rsd.seqpos() << " this_conn_id=" << this_conn_id << std::endl; //DELETE ME
			runtime_assert_string_msg(connect_map_size() >= this_conn_id, "Residue::update_connections_to_other_residue() error:  Connection id reported by other residue doesn't exist in current residue!");
			if ( connected_residue_at_resconn(this_conn_id)!=other_rsd.seqpos() ) {
				TR.Warning << "Warning!  While updating residue " << seqpos() << "'s connections to residue " << other_rsd.seqpos() << ", a connection to residue " << connected_residue_at_resconn(this_conn_id) << " was overwritten!" << std::endl;
			}
			residue_connection_partner( this_conn_id, other_rsd.seqpos(), ic ); //Set this residue's connection appropriately for the other residue's connection indices.
		}
	}
	return;
}

void
Residue::copy_residue_connections( Residue const & src_rsd )
{

	/// ASSUMPTION: if two residue types have the same number of residue connections,
	// then their residue connections are "the same" residue connections.
	// This assumption works perfectly for typical polymeric residues, but the following
	// assignment would produce unexpected behavior: take a CYD residue i that's disulfide
	// partner is residue j and replace it with a catalytic glutamate GLC that's bound
	// to ligand residue k.  Connection #3 for cyd will be confused with connection #3 on
	// GLC.  Such weird cases will have to be explicitly detected and handled.
	//
	// THIS REALLY NEEDS TO BE FIXED AT SOME POINT

	// AMW: coming back here, planning to add logic like second case to the first case.
	// let's assume at least pseudobonds have to be the same.

	if ( type().n_possible_residue_connections() == src_rsd.type().n_possible_residue_connections() ) {

		connect_map_ = src_rsd.connect_map_;
		connections_to_residues_ = src_rsd.connections_to_residues_;
		pseudobonds_ = src_rsd.pseudobonds_;

	} else {

		if ( ! src_rsd.pseudobonds_.empty() ) {
			std::cerr << "Unable to handle change in the number of residue connections in the presence of pseudobonds!" <<
				std::endl;
			utility_exit();
		}

		connect_map_.clear();
		connections_to_residues_.clear();
		pseudobonds_.clear();

		connect_map_.resize( type().n_possible_residue_connections() );

		// Find correspondence between src_rsd's connection atoms and atoms on *this.
		for ( Size ii = 1; ii <= src_rsd.type().n_possible_residue_connections(); ++ii ) {
			Size const ii_connatom = src_rsd.type().residue_connection( ii ).atomno();
			if ( has( src_rsd.atom_name( ii_connatom ) ) ) {

				Size const this_connatom = atom_index( src_rsd.atom_name( ii_connatom ));

				/// Simple case: this atom on both residues is connected to only a single other residue.
				if ( type().n_residue_connections_for_atom( this_connatom ) == 1 &&
						src_rsd.type().n_residue_connections_for_atom( ii_connatom ) == 1 ) {
					Size const this_connid = type().residue_connection_id_for_atom( this_connatom );

					// note that we might have the same atom name, but in src_rsd it's connected to something
					// and in our rsd it's not. So check for that now:
					if ( this_connid ) {
						residue_connection_partner(
							this_connid,
							src_rsd.connect_map( ii ).resid(),
							src_rsd.connect_map( ii ).connid() );
						if ( this_connid != ii ) {
							TR.Debug << "WARNING: Residue connection id changed when creating a new residue at seqpos " << seqpos() <<
								std::endl;
							TR.Debug << "WARNING: ResConnID info stored on the connected residue (residue " << src_rsd.connect_map( ii ).resid();
							TR.Debug << ") is now out of date!" << std::endl;
							TR.Debug << "Connection atom name (in src): " << src_rsd.atom_name( ii_connatom ) << std::endl;
						}
					}
				} else {
					/// Preserve residue connections in their input order.
					/// Figure out which residue connection for this atom on the source residue we're looking at.
					/// This *could* lead to weird behavior if you were to remove the middle of three residue connections
					/// for a single atom; e.g. if you had a Zn coordinated to three residues and wanted to replace it
					/// with a Zn coordinated to two residues -- then the question should be, which two residue connections
					/// from the original Zn should you copy.  At that point, in fact, you would have a weird situation
					/// where the residues coordinating the Zn would have out-of-date information about which residue connection
					/// on Zn they're coordinated to.
					/// The logic for altering residue connections in any way besides first building up a molecule
					/// is terribly incomplete.
					/// Fortunately, once a molecule is built and its residue connection topology is finalized, then all
					/// downstream operations are a sinch.
					Size which_connection_on_this_atom( 0 );
					for ( Size jj = 1; jj <= src_rsd.type().residue_connections_for_atom( ii_connatom ).size(); ++jj ) {
						if ( src_rsd.type().residue_connections_for_atom( ii_connatom )[ jj ] == ii ) {
							which_connection_on_this_atom = jj;
							break;
						}
					}
					if ( which_connection_on_this_atom == 0 ) {
						utility_exit_with_message("CATASTROPHIC ERROR in Residue::copy_residue_connections.  ResidueType connection map integrity error");
					}
					if ( which_connection_on_this_atom <= type().residue_connections_for_atom( this_connatom ).size() ) {
						residue_connection_partner(
							type().residue_connections_for_atom( this_connatom )[ which_connection_on_this_atom ],
							src_rsd.connect_map( ii ).resid(),
							src_rsd.connect_map( ii ).connid() );
					} else {
						/// Warn, we've just dropped a residue connection.  Was that intentional?
						/// Actually -- common occurrence when converting a mid-residue to a terminal residue.
						/// std::cerr << "WARNING: Not copying residue connection " << ii << " from " << src_rsd.name()
						/// << " to " << name() << " at position " << seqpos() << std::endl;
					}
				}
			}
		}
	}
}


/// @details loop over all actcoord atoms for this ResidueType,
/// average their actual positions in this residue.
void
Residue::update_actcoord()
{
	rsd_type_.update_actcoord( *this );
}

void Residue::select_orient_atoms(Size & center, Size & nbr1, Size & nbr2) const
{
	rsd_type_.select_orient_atoms(center, nbr1, nbr2);
}

/// @details  Helper function: selects atoms to orient on and transforms all of my atoms to
/// orient onto another residue. Used by place(). Need to think a bit more about the
/// restrictions on src...
void
Residue::orient_onto_residue( Residue const & src )
{
	using kinematics::Stub;

	Size center, nbr1, nbr2;
	select_orient_atoms( center, nbr1, nbr2 );
	// std::cout << " CENTER " << atom_name( center ) << "   NBR1 " << atom_name( nbr1 ) << "    NBR2 " << atom_name( nbr2 ) << std::endl;

	debug_assert( center && nbr1 && nbr2 );
	// this will fail if src doesnt have these atoms -- think more about this!
	//
	orient_onto_residue(
		src,
		center,
		nbr1,
		nbr2,
		src.atom_index( rsd_type_.atom_name( center )),
		src.atom_index( rsd_type_.atom_name( nbr1 )),
		src.atom_index( rsd_type_.atom_name( nbr2 )));

} // orient_onto_residue( Residue const & src)


void
Residue::orient_onto_residue(
	Residue const & src,
	utility::vector1< std::pair< std::string, std::string > > const & atom_pairs
)
{
	using kinematics::Stub;

	// Verify that three atom pairs have been provided
	if ( atom_pairs.size() != 3 ) {
		utility_exit_with_message( "Three atom pairs must be provided in Residue::orient_onto_residue.");
	}

	orient_onto_residue(
		src,
		atom_index( atom_pairs[1].second ),
		atom_index( atom_pairs[2].second ),
		atom_index( atom_pairs[3].second ),
		src.atom_index( atom_pairs[1].first ),
		src.atom_index( atom_pairs[2].first ),
		src.atom_index( atom_pairs[3].first ));
} //orient_onto_residue( Residue src, atom_pairs )

void Residue::orient_onto_residue(
	Residue const & src,
	Size center, Size nbr1, Size nbr2,
	Size src_center, Size src_nbr1, Size src_nbr2)
{
	using kinematics::Stub;

	//NOTE: the implementation of this function might change in the future
	//from strictly superimposing on three atoms to superposition along the lines
	//of what is in numeric::model_quality::findUU()

	// explanation for taking the midpoint...?
	Vector const
		rot_midpoint ( 0.5 * (     atom(     nbr1 ).xyz() +     atom(     nbr2 ).xyz() ) ),
		src_midpoint ( 0.5 * ( src.atom( src_nbr1 ).xyz() + src.atom( src_nbr2 ).xyz() ) );

	Stub rot_stub( atom( center ).xyz(),
		rot_midpoint,
		atom( nbr1 ).xyz() );

	Stub src_stub( src.atom( src_center ).xyz(),
		src_midpoint,
		src.atom( src_nbr1 ).xyz() );

	// this could be made faster by getting the composite rotation and translation

	for ( Size i=1; i<= rsd_type_.natoms(); ++i ) {
		Vector const old_xyz( atoms()[i].xyz() );
		Vector const new_xyz( src_stub.local2global( rot_stub.global2local( old_xyz ) ) );
		atoms()[i].xyz( new_xyz );
	}
}


/// @details oritent onto residue for peptoids
/// Not to proud of this as I think there is probably a more general way but at the moment it is alluding me
/// To align a peptoid the N is the center and the nbrs are the lower connect and the CA. The lower connect is
/// not stored in the atom index so we need to get the xyz coord a different way.
void
Residue::orient_onto_residue_peptoid (
	Residue const & src,
	Conformation const & conformation
)
{
	using kinematics::Stub;
	using namespace core;
	using namespace conformation;
	using namespace chemical;

	// // DEBUG DEBUG DEBUG
	//  using namespace pose;
	//  Pose before_pose;
	//  before_pose.append_residue_by_jump( *this, 1 );
	//  std::string before_filename("before.pdb");
	//  before_pose.dump_pdb( before_filename );

	// static strings to cut down on object creation penelty
	static std::string const bb_ca( " CA " );
	static std::string const bb_co( " C  " );
	static std::string const bb_nx( " N  " );
	static std::string const bb_ac( " CO " );
	static std::string const sc_ca( " CA1" );

	// xyz Vectors
	Vector cent_xyz, nbr1_xyz, nbr2_xyz, src_cent_xyz, src_nbr1_xyz, src_nbr2_xyz;

	// The position of the first sidechain atom in the first residue of a chain of a peptoid is not defined as it is for peptides
	// unless it is part of a macrocycle or it is acetylated as both of those (just like residues in the middle and upper terminus
	// of chains) provide a preceding atom to determin the position of the first sidechain atom.

	// the order that these are called in is important

	//std::cout << "DEBUG DEBUG DEBUG " << src.type().name() << " " << src.seqpos() << " " << src.type().is_lower_terminus() << " "
	//          << ( src.has_variant_type( chemical::NTERM_CONNECT ) && conformation.residue( conformation.chain_end( src.chain() ) ).has_variant_type( chemical::CTERM_CONNECT ) )
	//          << std::endl;

	// AMW: I am recombining these conditionals in a relatively byzantine way to attempt to maintain their order
	// I don't know the last time they should have worked, but various checks in here break in the current state
	// first residue of chain
	//TR << " orienting " << src.type().name() << " onto " << type().name() << std::endl;
	if ( src.type().is_lower_terminus() || src.has_variant_type( chemical::NTERM_CONNECT ) ) {
		//std::cout << "DEBUG LT" << std::endl;
		// NtermConnect (cyclic)
		//std::cout << "Num chains is " << conformation.num_chains() << " src.chain() is " << src.chain() << " conformation.size() is " << conformation.size() << " this chain end is " << conformation.chain_end( src.chain() )  << std::endl;
		Residue const cterm_residue = ( conformation.num_chains() == 2 || src.chain() == conformation.num_chains() - 1 ) ?
			conformation.residue( conformation.size() ) :
			conformation.residue( conformation.chain_end( src.chain() ) );
		if ( /*!(*this).type().is_lower_terminus() && */src.has_variant_type( chemical::NTERM_CONNECT ) && cterm_residue.has_variant_type( chemical::CTERM_CONNECT ) ) {
			//std::cout << "DEBUG LT CYCLIC" << std::endl;
			// cent: peptoid N, src N
			cent_xyz = atom( atom_index( bb_nx ) ).xyz();
			src_cent_xyz = src.atom( src.atom_index( bb_nx ) ).xyz();

			// nbr2; peptoid CA, src CA
			nbr2_xyz = atom( atom_index( bb_ca ) ).xyz();
			src_nbr2_xyz = src.atom( src.atom_index( bb_ca ) ).xyz();

			// nbr1: peptoid lower connect, C of the cyclic paterner of src
			//Residue const cyclic_src( conformation.residue( conformation.chain_end( src.chain() ) ) );
			src_nbr1_xyz = cterm_residue.atom( cterm_residue.atom_index( bb_co ) ).xyz();

			Stub stub( atom( atom_index( bb_nx ) ).xyz(), atom( atom_index( bb_ca ) ).xyz(), atom( atom_index( bb_co ) ).xyz() );
			//AtomICoor icoor ( (*this).type().lower_connect().icoor() );
			// Put it where the lower connect WOULD be put...
			// we can't actually summon lower connect because it's a lower terminal type!
			// AMW: nonsense! we have cterm_residue for that!
			nbr1_xyz = cterm_residue.atom( cterm_residue.atom_index( bb_co ) ).xyz(); //stub.spherical( 167.209, 52.954, 1.367000 );

			// AcetylatedPeptoidNterm
		} else if ( src.has_variant_type( chemical::ACETYLATED_NTERMINUS_VARIANT ) ) {
			//std::cout << "DEBUG LT ACE" << std::endl;
			// cent: peptoid N, src N
			cent_xyz = atom( atom_index( bb_nx ) ).xyz();
			src_cent_xyz = src.atom( src.atom_index( bb_nx ) ).xyz();

			// nbr2; peptoid CA, src CA
			nbr2_xyz = atom( atom_index( bb_ca ) ).xyz();
			src_nbr2_xyz = src.atom( src.atom_index( bb_ca ) ).xyz();

			// nbr1: peptoid CO, src CO ( CO is the atom that prepends the N in the acetylated varient type)
			nbr1_xyz = atom( atom_index( bb_ac ) ).xyz();
			src_nbr1_xyz = src.atom( src.atom_index( bb_ac ) ).xyz();

			// NtermPeptoidFull ( peptoid on peptoid ). All peptoids have at least one sidechain atom called CA1
		} else if ( src.type().is_peptoid() ) {
			// std::cout << "DEBUG LT PEPTOID" << std::endl;
			// cent: peptoid N, src N
			cent_xyz = atom( atom_index( bb_nx ) ).xyz();
			src_cent_xyz = src.atom( src.atom_index( bb_nx ) ).xyz();

			// nbr2; peptoid CA, src CA
			nbr2_xyz = atom( atom_index( bb_ca ) ).xyz();
			src_nbr2_xyz = src.atom( src.atom_index( bb_ca ) ).xyz();

			// nbr1: peptoid CA1, src CA1
			nbr1_xyz = atom( atom_index( sc_ca ) ).xyz();
			src_nbr1_xyz = src.atom( src.atom_index( sc_ca ) ).xyz();

			// NtermPeptoidFull ( peptoid on to peptide ). Nothing to base peptoid sidechain off of so just use ideal internal coords from params
		} else {
			// std::cout << "DEBUG LT ELSE" << std::endl;
			// cent: peptoid CA, src CA
			cent_xyz = atom( atom_index( bb_ca ) ).xyz();
			src_cent_xyz = src.atom( src.atom_index( bb_ca ) ).xyz();

			// nbr1: peptoid N, src N
			nbr1_xyz = atom( atom_index( bb_nx ) ).xyz();
			src_nbr1_xyz = src.atom( src.atom_index( bb_nx ) ).xyz();

			// nbr2; peptoid C, src C
			nbr2_xyz = atom( atom_index( bb_co ) ).xyz();
			src_nbr2_xyz = src.atom( src.atom_index( bb_co ) ).xyz();
		}
	} else { // not first reisdue of chain
		//std::cout << "DEBUG ELSE" << std::endl;
		// cent: peptoid N, src N
		cent_xyz = atom( atom_index( bb_nx ) ).xyz();
		src_cent_xyz = src.atom( src.atom_index( bb_nx ) ).xyz();

		// nbr2; peptoid CA, src CA
		nbr2_xyz = atom( atom_index( bb_ca ) ).xyz();
		src_nbr2_xyz = src.atom( src.atom_index( bb_ca ) ).xyz();

		// nbr1: peptoid lower connect, C of the residue before src
		Residue const prev_src( conformation.residue( src.seqpos() - 1 ) );
		src_nbr1_xyz = prev_src.atom( prev_src.atom_index( bb_co ) ).xyz();

		Stub stub( atom( atom_index( bb_nx ) ).xyz(), atom( atom_index( bb_ca ) ).xyz(), atom( atom_index( bb_co ) ).xyz() );
		AtomICoor icoor ( (*this).type().lower_connect().icoor() );
		nbr1_xyz = stub.spherical(icoor.phi(), icoor.theta(), icoor.d() );
	}

	//  std::cout
	//  << "rot_cent_xyz " << cent_xyz[1] << " " << cent_xyz[2] << " " << cent_xyz[3] << std::endl
	//  << "rot_nbr1_xyz " << nbr1_xyz[1] << " " << nbr1_xyz[2] << " " << nbr1_xyz[3] << std::endl
	//  << "rot_nbr2_xyz " << nbr2_xyz[1] << " " << nbr2_xyz[1] << " " << nbr2_xyz[1] << std::endl
	//  << "src_cent_xyz " << src_cent_xyz[1] << " " << src_cent_xyz[2] << " " << src_cent_xyz[3] << std::endl
	//  << "src_nbr1_xyz " << src_nbr1_xyz[1] << " " << src_nbr1_xyz[2] << " " << src_nbr1_xyz[3] << std::endl
	//  << "src_nbr2_xyz " << src_nbr2_xyz[1] << " " << src_nbr2_xyz[2] << " " << src_nbr2_xyz[3] << std::endl;

	// get midpoint
	Vector rot_midpoint( 0.5 * ( nbr1_xyz + nbr2_xyz ) );
	Vector src_midpoint( 0.5 * ( src_nbr1_xyz + src_nbr2_xyz ) );

	// get stubs
	Stub rot_stub( cent_xyz, rot_midpoint, nbr1_xyz );
	Stub src_stub( src_cent_xyz, src_midpoint, src_nbr1_xyz);

	// orient
	for ( Size i=1; i<= rsd_type_.natoms(); ++i ) {
		Vector const old_xyz( atoms()[i].xyz() );
		Vector const new_xyz( src_stub.local2global( rot_stub.global2local( old_xyz ) ) );
		atoms()[i].xyz( new_xyz );
		// DEBUG DEBUG DEBUG
		//std::cout << "old " << i << " " << old_xyz[1] << " " << old_xyz[2] << " " << old_xyz[3] << std::endl;
		//std::cout << "new " << i << " " << new_xyz[1] << " " << new_xyz[2] << " " << new_xyz[3] << std::endl;
	}

}

/// @details Place/orient "this" Residue onto "src" Residue by backbone superimposition
/// Since rotamer is represented by Residue in mini now, this function is mainly used to place a rotamer
/// onto the backbone of "src" residue. Meanwhile, it can also be used to add sidechains to one pose/conformation
/// from another pose/conformation.\n
/// current logic: find backbone atom with bonded neighbors in sidechain,
/// and which is the base_atom of those neighbors. Take that backbone atom
/// and find two neighboring backbone heavyatoms. The three atoms to be superimposed with
/// are the center/base atom, the backbone neighbor 1 and the mid-point of backbone neighbor 1 and 2. This way,
/// we can avoid large perturbation on backbone neighbor 2 after superimpostion if the two sets of backbone atoms
/// are not perfectly superimposable ( e.g., with slightly different backbone geometry).\n
/// after all atoms in "this" Residue is oriented, copy any corresponding backbone atom coords from "src" and if
/// there are any backbone atom missing from "src" (for example, src is a proline with HN missing), build them using
/// ideal internal coords (that is why "conformation" is needed as an input argument).
/// For residues without any backbone atoms (e.g. some ligands), center on nbr_atom instead
/// and two of its bonded neighbors (preferring heavy atoms to hydrogens if possible).
void
Residue::place( Residue const & src, Conformation const & conformation, bool preserve_c_beta )
{
	using kinematics::Stub;
	using namespace basic::options;

	Size first_scatom( rsd_type_.first_sidechain_atom() );
	if ( first_scatom >= 1 && first_scatom <= rsd_type_.nheavyatoms() ) {
		// not all backbone -- need to do some orienting
		if ( (*this).type().is_peptoid() ) {
			orient_onto_residue_peptoid( src, conformation );
		} else {
			orient_onto_residue( src );
		}
	} // does the residue have any sidechain atoms?

	// now copy the backbone atoms
	//
	utility::vector1< bool > missing( natoms(), false );
	bool any_missing( false );
	for ( Size i=1; i<= natoms(); ++i ) {

		//The O2' is a special case for RNA, because it is a "sidechain" atom that
		// branches off the backbone separately from the base. This could be
		// coded more robustly by modifying orient_onto_residue() properly.

		if ( !rsd_type_.atom_is_backbone(i) &&
				!(rsd_type_.atom_name(i) == " O2'") &&
				!(rsd_type_.atom_name(i) == "HO2'") &&
				!( atom_depends_on_lower(i, rsd_type_.has_polymer_dependent_groups()) ) &&
				!( atom_depends_on_upper(i, rsd_type_.has_polymer_dependent_groups()) )
				) {
			continue;
		}
		if ( src.has( rsd_type_.atom_name(i) ) &&
				( (!rsd_type_.is_polymer() ) ||
				(!rsd_type_.atom_is_hydrogen(i)) ||
				rsd_type_.icoor(i).depends_on_polymer_lower() ||
				rsd_type_.icoor(i).depends_on_polymer_upper() ||
				(rsd_type_.is_lower_terminus() && rsd_type_.icoor(i).stub_atom1().atomno() == 1)
				/*Force rebuild of mainchain hydrogens that don't depend on polymer lower or upper and aren't N-terminal amide protons.  Logic below ensures that hydrogen positions are preserved if they don't deviate by more than a threshhold from ideal. VKM 21 Aug 2015.*/
				)
				) {
			atoms()[i].xyz( src.atom( src.atom_index( rsd_type_.atom_name(i) ) ).xyz() );
		} else {
			missing[i] = true;
			any_missing = true;
		}
	}

	if ( any_missing ) fill_missing_atoms( missing, conformation );

	//Check mainchain hydrogens, which are currently in idealized positions distinct from the input position.
	//If they deviate from less than a threshhold value from the input position, preserve the input position;
	//otherwise, keep the idealized position.
	if ( rsd_type_.is_polymer() ) {
		core::Real const thresh_sq( static_cast<core::Real>( pow(static_cast<core::Real>(option[OptionKeys::packing::mainchain_h_rebuild_threshold]()),2.0) ) );
		for ( core::Size i=1, imax=natoms(); i<=imax; ++i ) {
			if ( !rsd_type_.atom_is_hydrogen(i) ) continue; //Skip non-hydrogen atoms.
			if ( !rsd_type_.atom_is_backbone(i) && !(rsd_type_.atom_name(i) == " O2'") && !(rsd_type_.atom_name(i) == "HO2'") ) continue; //Skip non-backbone atoms, preserving special-case RNA atoms.
			if ( !src.has( rsd_type_.atom_name(i) ) ) continue; //Skip if the source residue doesn't have this atom type.

			// If the idealized position is less than the threshold from the source position, or it's polymer connection-dependent...
			if ( atoms()[i].xyz().distance_squared( src.atom( src.atom_index( rsd_type_.atom_name(i) ) ).xyz() ) < thresh_sq ||
					rsd_type_.icoor(i).depends_on_polymer_lower() ||
					rsd_type_.icoor(i).depends_on_polymer_upper() ||
					(rsd_type_.is_lower_terminus() && rsd_type_.icoor(i).stub_atom1().atomno() == 1) //Damn -- sort of special case: protons on n-terminus
					) {
				atoms()[i].xyz( src.atom( src.atom_index( rsd_type_.atom_name(i) ) ).xyz() ); // ...then move it back to the source position
			}
		}
	}

	if ( preserve_c_beta ) {
		// after superposition, adjust the sidechain atoms by aligning the CA-CB bond
		//  std::cout << "preserving c-beta... " << std::endl;
		std::string root("CA"), mobile_new("CB"), mobile_src("CB");

		if ( is_RNA() ) {
			root = " C1'";
			mobile_new = atom_name( chi_atoms( 1 )[ 3 ] ); //First atom in base...
			mobile_src = src.atom_name( src.chi_atoms( 1 )[ 3 ] ); //First atom in base...
		}

		// only try this when both the src and the new residue types contain both atom types
		if ( type().has( root ) && type().has( mobile_new ) &&
				src.type().has( root ) && src.type().has( mobile_src ) ) {
			debug_assert( xyz( root ) == src.xyz( root ) ); // roots should be aligned by now
			// common 'pseudoatom' vector, perpendicular to the plane defined by the two bonds
			Vector const pseudoatom(
				cross( xyz( mobile_new ) - xyz( root ), src.xyz( mobile_src ) - src.xyz( root ) ) + xyz( root )
			);

			if ( pseudoatom==xyz(root) ) return;

			Stub new_stub(     xyz( root ), pseudoatom,     xyz( mobile_new ) ),
				src_stub( src.xyz( root ), pseudoatom, src.xyz( mobile_src ) );
			// adjust sidechain coordinates by superposition of the bond 'stubs'
			// would need a smarter way to propagate through all child atoms of 'root' to generalize this
			for ( Size atom_index(1); atom_index <= type().natoms(); ++atom_index ) {
				if ( type().atom_is_backbone( atom_index ) ) continue;
				//special case for RNA
				if ( type().atom_name( atom_index ) == " O2'" || type().atom_name( atom_index ) == "HO2'" ) continue;
				Vector const old_xyz( atoms()[ atom_index ].xyz() );
				Vector const new_xyz( src_stub.local2global( new_stub.global2local( old_xyz ) ) );
				atoms()[ atom_index ].xyz( new_xyz );
			}
		}
	}
}

bool
missing_stubs_build(core::Size ii, Residue const & residue, utility::vector1< bool > const & missing, Vector & coordinate ) {
	if ( ! residue.icoor(ii).is_internal() ) { return false; }
	std::string const & ii_name( residue.atom_name(ii) );
	chemical::AtomICoor const & ic( residue.icoor(ii) );
	chemical::ICoorAtomID const & id1( ic.stub_atom(1) );

	// Let's build the icoord from bits and pieces.
	// Note that we only pull from strictly internal coordinates - this ensures that the atomno values corrspond to real atoms.
	// First, possible distances
	utility::vector1< chemical::AtomICoor > distances; // These are valid for distances only
	if ( ! missing[ id1.atomno() ] && id1.atomno() != ii ) { // The current distance could work - put this first to prefer this building method.
		chemical::AtomICoor newicoor( ii_name, 0.0, 0.0, ic.d(),
			residue.atom_name(id1.atomno()), ii_name, ii_name, residue.type() ); // ii's are placeholders
		distances.push_back( newicoor );
	}
	for ( Size jj=1; jj <= residue.natoms(); ++jj ) {
		if ( ! missing[jj] ) {
			chemical::AtomICoor const & icj( residue.icoor(jj) );
			if ( ! icj.is_internal() ) { continue; }
			// Looking for JJ-II distances to build II-JJ
			if ( icj.stub_atom(1).atomno() == ii && jj != ii ) {
				chemical::AtomICoor newicoor( ii_name, 0.0, 0.0, icj.d(),
					residue.atom_name(jj), ii_name, ii_name, residue.type() );
				distances.push_back( newicoor );
			}
		}
	}
	// Now we find possible angles for the distance pairs.
	utility::vector1< chemical::AtomICoor > angles; // These are valid for distance and bondangles only
	for ( core::Size aa=1; aa <= distances.size(); ++aa ) {
		chemical::AtomICoor const & previc( distances[aa] );
		core::Size dd( previc.stub_atom(1).atomno() );
		for ( Size jj=1; jj <= residue.natoms(); ++jj ) {
			chemical::AtomICoor const & icj( residue.icoor(jj) );
			if ( ! icj.is_internal() ) { continue; }
			// We either want II-DD-XX or JJ-DD-II
			if ( icj.stub_atom1().atomno() == dd ) {
				core::Size xx(icj.stub_atom2().atomno());
				if ( jj == ii && ! missing[xx] && xx != ii && xx != dd ) { // II-DD-XX
					chemical::AtomICoor newicoor( ii_name, 0.0, icj.theta(), previc.d(),
						residue.atom_name(dd),
						residue.atom_name(xx),
						ii_name, residue.type() );
					angles.push_back( newicoor );
				} else if ( !missing[jj] && xx == ii && jj != ii && jj != dd ) { // JJ-DD-II
					chemical::AtomICoor newicoor( ii_name, 0.0, icj.theta(), previc.d(),
						residue.atom_name(dd),
						residue.atom_name(jj),
						ii_name, residue.type() );
					angles.push_back( newicoor );
				}
			}
		}
	}
	// Now we find possible dihedrals consistent with valids angles and dihedrals
	utility::vector1< chemical::AtomICoor > dihedrals; // These are fully valid
	for ( core::Size pp=1; pp <= angles.size(); ++pp ) {
		chemical::AtomICoor const & previc( angles[pp] );
		core::Size dd( previc.stub_atom(1).atomno() );
		core::Size aa( previc.stub_atom(2).atomno() );
		for ( Size jj=1; jj <= residue.natoms(); ++jj ) {
			chemical::AtomICoor const & icj( residue.icoor(jj) );
			if ( ! icj.is_internal() ) { continue; }
			// Four cases:
			//    II-DD-AA-XX - this doesn't work, or else we could build it straight
			//    JJ-AA-DD-II - reverse orientation - the same dihedral
			//    JJ-DD-AA-II - We can still steal the dihedral, we just need to negate it.
			//    II-AA-DD-XX - Like the JJ-DD-AA-II case
			if ( jj == ii ) {
				core::Size xx( icj.stub_atom3().atomno() );
				if ( !missing[xx] && icj.stub_atom1().atomno() == aa && icj.stub_atom2().atomno() == dd
						&& xx != ii && xx != dd && xx != aa ) { // II-AA-DD-XX
					chemical::AtomICoor newicoor( ii_name, -1 * icj.phi(), previc.theta(), previc.d(),
						residue.atom_name(dd),
						residue.atom_name(aa),
						residue.atom_name(xx), residue.type() );
					dihedrals.push_back( newicoor );
				}
			} else if ( ! missing[jj] ) {
				if ( icj.stub_atom1().atomno() == aa && icj.stub_atom2().atomno() == dd && icj.stub_atom3().atomno() == ii
						&& jj != ii && jj != dd && jj != aa ) { // JJ-AA-DD-II
					chemical::AtomICoor newicoor( ii_name, icj.phi(), previc.theta(), previc.d(),
						residue.atom_name(dd),
						residue.atom_name(aa),
						residue.atom_name(jj), residue.type() );
					dihedrals.push_back( newicoor );
				} else if ( icj.stub_atom1().atomno() == dd && icj.stub_atom2().atomno() == aa && icj.stub_atom3().atomno() == ii
						&& jj != ii && jj != dd && jj != aa ) { // JJ-DD-AA-II
					chemical::AtomICoor newicoor( ii_name, -1 * icj.phi(), previc.theta(), previc.d(),
						residue.atom_name(dd),
						residue.atom_name(aa),
						residue.atom_name(jj), residue.type() );
					dihedrals.push_back( newicoor );
				}
			}
		}
	}

	if ( dihedrals.size() >= 1 ) {
		// We have multiple ways of building the atom, we can just use the first one.
		TR.Debug << "Building atom " << ii_name << " based on assembled internal coordinates." << std::endl;
		coordinate = dihedrals[1].build(residue);
		return true;
	}

	// Building failed.
	return false;

}

/// @brief Build
void
improper_build(Residue const & residue,
	core::Size missing,
	core::Size parent,
	core::Size sibling1,
	core::Size sibling2 ,
	Vector & coordinate
) {
	core::chemical::ResidueType const & restype( residue.type() );
	core::Vector to_missing( restype.atom(missing).ideal_xyz() - restype.atom(parent).ideal_xyz() );
	core::Vector to_sib1( restype.atom(sibling1).ideal_xyz() - restype.atom(parent).ideal_xyz());
	core::Real d( to_missing.length() );
	core::Real theta( numeric::constants::r::pi - angle_of( to_missing, to_sib1 ) );
	core::Real phi( numeric::dihedral_radians( restype.atom(missing).ideal_xyz(), restype.atom(parent).ideal_xyz(),
		restype.atom(sibling1).ideal_xyz(), restype.atom(sibling2).ideal_xyz() ) );
	chemical::AtomICoor newicoor( residue.atom_name(missing), phi, theta, d,
		residue.atom_name(parent),
		residue.atom_name(sibling1),
		residue.atom_name(sibling2), restype);
	coordinate = newicoor.build(residue);
}

/////////////////////////////////////////////////////////////////////////////
/// @details
/// this uses ideal internal coords to build any missing atom from its three
/// stub atoms. If any of the stub atoms are missing, build them first.
/// Unable to build a missing atom whose stub atoms are from non-existing
/// polymer connection and its input bogus value will not be changed.
bool
Residue::fill_missing_atoms(
	utility::vector1< bool > & missing,
	Conformation const & conformation,
	bool fail /* = true */
)
{
	bool still_missing( true );
	bool progress( false );
	core::Size desperation( 0 );
	// This will never be an infinite loop, because each time through we either
	// 1) turn at least one missing[i] false
	// 2) increase the desperation level
	// and we exit whenever
	// 1) all missing[i] are false, or
	// 2) the desperation level gets too high
	while ( still_missing ) {
		still_missing = false;
		progress = false;
		for ( Size i=1; i<= natoms(); ++i ) {
			if ( missing[i] ) {
				chemical::AtomICoor const & ic( icoor(i) );

				// check to see if any of our stub atoms are missing:
				bool stub_atoms_missing( false );
				bool stubs_buildable( true );
				for ( Size j=1; j<= 3; ++j ) {
					chemical::ICoorAtomID const & id( ic.stub_atom(j) );
					// We assume all connection points to other residues are not missing
					if ( id.type() == chemical::ICoorAtomID::INTERNAL && missing[ id.atomno() ] ) {
						stub_atoms_missing = true;
					}
					if ( ! id.buildable( *this, conformation ) ) {
						stubs_buildable = false;
					}
				}

				if ( !stub_atoms_missing ) {
					// no stub atoms missing: build coordinates for this atom from them
					if ( ! stubs_buildable ) {
						// We are dependant on residue connections which might not have coordinates.
						// (e.g. like non-termini variants at the beginning/end of the pose)
						// NOTE: This used to check for residues next to chainbreaks too, but I (RM) removed that logic,
						// because it made dodgy assumptions about the ordering of residues w/r/t chainbreaks
						TR.Warning << "[ WARNING ] missing an atom: " << seqpos_ << " " << atom_name(i) << " that depends on a nonexistent polymer connection! "
							<< std::endl <<  " --> generating it using idealized coordinates." << std::endl;
						set_xyz( i, ic.build(*this));
					} else {
						TR.Debug << "Building atom " << atom_name(i) << " based on standard internal coordinates." << std::endl;
						set_xyz( i, ic.build(*this, conformation ) );
					}
					missing[i] = false;
					progress = true;
					continue;
				} else if ( !progress && desperation >= 2 ) {
					// Only try building for missing stubs if we have to.
					Vector coordinates;
					if ( missing_stubs_build(i, *this, missing, coordinates ) ) {
						missing[i] = false;
						set_xyz(i, coordinates);
						progress = true;
						desperation = 0; // Try to extend this without further without-stub-building.
						continue;
					}
				}

				still_missing = true;
			} else if ( desperation >= 1 ) {
				// This atom is present. If we have two bonded atoms also present,
				// we can build any missing bonded atoms from the ideal geometry around this atom.
				AtomIndices const & bonded( bonded_neighbor(i) );
				if ( bonded.size() >= 3 ) { // Need at least two present and one not present
					utility::vector1<core::Size> nbr_present;
					utility::vector1<core::Size> nbr_missing;
					for ( core::Size bb(1); bb <= bonded.size(); ++bb ) {
						if ( missing[bonded[bb]] ) {
							nbr_missing.push_back( bonded[bb] );
						} else {
							nbr_present.push_back( bonded[bb] );
						}
					}
					if ( nbr_present.size() >= 2 && nbr_missing.size() >= 1 ) {
						for ( core::Size mm(1); mm <= nbr_missing.size(); ++mm ) {
							Vector coordinates;
							improper_build(*this, nbr_missing[mm], i, nbr_present[1], nbr_present[2], coordinates );
							missing[ nbr_missing[mm] ] = false;
							set_xyz(nbr_missing[mm], coordinates);
						}
						progress = true;
						desperation = 0; // Try to extend this without further improper building
						continue;
					}
				}
			}
		} // end for i in atoms
		if ( ! progress ) {
			++desperation;
			if ( desperation >= 3 ) {
				// Did our best to build the atoms - simply can't.
				if ( fail ) {
					TR.Error << "[ ERROR ] Cannot build coordinates for residue " << name() << " at position " << seqpos() << ": missing too many atoms." << std::endl;
					type().show_all_atom_names(TR.Debug);
					TR.Debug << "Internal coordinate tree:" << std::endl;
					core::chemical::pretty_print_atomicoor( TR.Debug, type().icoor(type().atom_index(type().root_atom())), type());
					TR.Error << "Missing atoms are: ";
					for ( core::Size nn(1); nn <= missing.size(); ++nn ) {
						if ( missing[nn] ) {
							TR.Error << atom_name(nn) << "  ";
						}
					}
					TR.Error << std::endl;
					utility_exit_with_message("Unable to fill in missing atoms.");
				} else {
					return false;
				}
			}
		} // end no progress
	} // end while still missing
	return true;
}

void
Residue::clear_residue_connections()
{
	for ( Size ii = 1; ii <= connect_map_.size(); ++ii ) {
		//connect_map_[ ii ].resid( 0 );
		//connect_map_[ ii ].connid( 0 );
		connect_map_[ ii ].mark_incomplete();
	}
	connections_to_residues_.clear();
	pseudobonds_.clear();
	nonstandard_polymer_ = false;
}

void
Residue::copy_residue_connections_from( Residue const & src )
{
	this->nonstandard_polymer_ = src.nonstandard_polymer_;
	this->connect_map_ = src.connect_map_;
	this->connections_to_residues_ = src.connections_to_residues_;
	this->pseudobonds_ = src.pseudobonds_;
}


bool
Residue::has_incomplete_connection() const
{
	for ( Size ii = 1; ii <= connect_map_.size(); ++ii ) {
		if ( connection_incomplete( ii ) ) return true;
	}
	return false;
}


/// @details
/// determine whether an atom is completely connected to all possible bonded partners
bool
Residue::has_incomplete_connection(
	Size const atomno
) const
{
	Size const num_connections(n_possible_residue_connections());

	for ( Size i = 1; i <= num_connections; ++i ) {
		if ( residue_connect_atom_index(i) == atomno && connection_incomplete(i) ) return true;
	}

	return false;
}


bool
Residue::connection_incomplete( Size resconnid ) const
{
	return connect_map_[ resconnid ].incomplete();
}

id::AtomID
Residue::inter_residue_connection_partner(
	int connid,
	Conformation const & conformation
) const {
	if ( is_lower_terminus() ) --connid;
	if ( is_upper_terminus() ) --connid;

	Size const partner_seqpos( residue_connection_partner( connid ) );
	if ( partner_seqpos < 1 || partner_seqpos > conformation.size() ) {
		TR.Warning << "WARNING Residue::inter_residue_connection_partner: Invalid residue connection, returning BOGUS ID: this_rsd= " << name() <<
			' ' << seqpos() << " connid= " << connid << " partner_seqpos= " << partner_seqpos << std::endl;
		return id::BOGUS_ATOM_ID;
	}

	Size const partner_connid( residue_connection_conn_id( connid ) );
	Size const partner_atomno( conformation.residue_type( partner_seqpos ).residue_connect_atom_index( partner_connid ) );
	return id::AtomID( partner_atomno, partner_seqpos );
}


/// @details
/// set a connection to this residue by adding its partner's residue number
void
Residue::residue_connection_partner(
	Size const resconn_index, // ie, our connid
	Size const otherres,
	Size const other_connid
)
{
	connect_map_[ resconn_index ].resid(otherres);
	connect_map_[ resconn_index ].connid( other_connid );
	update_connections_to_residues();
	//  utility::vector1< Size > newlist;
	//  if (  connections_to_residues_.find( otherres ) != connections_to_residues_.end() ) {
	//   newlist = connections_to_residues_[ otherres ];
	//  }
	//  if ( newlist.size() != 0 ) {
	//   for ( Size ii = 1; ii <= newlist.size(); ++ii ) {
	//    if ( newlist[ ii ] == resconn_index  ) {
	//     //std::cout << "Setting residue connection partner on residue " << seqpos_ << " to residue " << otherres << " twice!" << std::endl;
	//     break;
	//    }
	//    else if ( ii == newlist.size() ) {
	//     newlist.push_back( resconn_index );
	//     connections_to_residues_[ otherres ] = newlist;
	//     break;
	//    }
	//   }
	//  } else {
	//   newlist.push_back( resconn_index );
	//   connections_to_residues_[ otherres ] = newlist;
	//  }
	determine_nonstandard_polymer_status();
}


/// @details  Private function to keep the connections_to_residues_ array up to date
/// @note  This could be made faster -- connections_to_residues_ is a std::map< > so operator[] calls are slow
void
Residue::update_connections_to_residues()
{
	connections_to_residues_.clear();
	for ( Size i=1, i_end = n_possible_residue_connections(); i<= i_end; ++i ) {
		Size const other_resid( connect_map_[ i ].resid() );
		connections_to_residues_[ other_resid ].push_back( i );
	}
}

/// @details update sequence numbers for this residue and
/// the numbers stored about its connections.
/// called by our owning conformation when the
/// sequence numbers are remapped
void
Residue::update_sequence_numbering( utility::vector1< Size > const & old2new )
{
	seqpos_ = old2new[ seqpos_ ];
	//std::map< Size, utility::vector1< Size > > connections_to_residues_copy( connections_to_residues_ );

	//connections_to_residues_.clear();
	for ( Size i=1, ie= connect_map_.size(); i<= ie; ++i ) {
		Size const old_resid = connect_map_[ i ].resid();
		if ( old_resid == 0 ) continue;

		Size const new_resid = old2new[ old_resid ];

		connect_map_[i].resid( new_resid );

		// If the partner disappears, partner atomid should be zero too. Otherwise if you add and
		// then delete a residue to a pose, a neighboring residue does not stay invariant.
		if ( new_resid == 0 ) connect_map_[i].connid( 0 );


		//   if ( new_resid ) {
		//    connections_to_residues_[ new_resid ] = connections_to_residues_copy[ old_resid ];
		//   }
	}
	update_connections_to_residues();
	if ( ! pseudobonds_.empty() ) {
		std::map< Size, PseudoBondCollectionCOP > copy_pseudobonds( pseudobonds_ );
		pseudobonds_.clear();
		for ( std::map< Size, PseudoBondCollectionCOP >::const_iterator
				pb_iter = copy_pseudobonds.begin(),
				pb_iter_end = copy_pseudobonds.end();
				pb_iter != pb_iter_end; ++pb_iter ) {
			Size old_neighbor_resid = pb_iter->first;
			Size new_neighbor_resid = old2new[ old_neighbor_resid ];
			if ( ! new_neighbor_resid ) continue;
			pseudobonds_[ new_neighbor_resid ] = pb_iter->second->clone_with_new_sequence_numbering( old2new );
		}
	}

	determine_nonstandard_polymer_status();
}

void
Residue::update_nus() {
	Size const n_nus( rsd_type_.n_nus() );
	for ( uint i( 1 ); i <= n_nus; ++i ) {
		AtomIndices const & nu_atoms( rsd_type_.nu_atoms( i ) );

		// Calculate the current nu angle from the coordinates.
		Angle const current_nu( numeric::dihedral_degrees(
			atom( nu_atoms[ 1 ] ).xyz(),
			atom( nu_atoms[ 2 ] ).xyz(),
			atom( nu_atoms[ 3 ] ).xyz(),
			atom( nu_atoms[ 4 ] ).xyz() ) );

		nus_[ i ] = current_nu;
	}
}


// Distance between a potential residue connection match and the position of the expected atom
Distance
Residue::connection_distance(
	conformation::Conformation const & conf,
	Size const resconn_index,
	Vector const & matchpoint
) const
{
	Vector ipos = type().residue_connection( resconn_index ).icoor().build( *this, conf );
	//std::cout << "ipos for " << name() << "'s connection atom " << resconn_index;
	//std::cout << ": ( " << ipos.x() << ", " << ipos.y() << ", " << ipos.z() << ")" << std::endl;
	return ipos.distance( matchpoint );
}


/// @details first check if it is polymer upper or lower connected to the other residue.
/// then check if it is boned to the other residue through non-polymer connection.
bool
Residue::is_bonded( Residue const & other ) const
{
	// trying this simpler strategy -- does not require that we keep chain id's in sync with chemical
	// connectivity
	// APL says shouldn't be much slower
	return ( connections_to_residues_.find( other.seqpos() ) != connections_to_residues_.end() );
}

bool
Residue::is_bonded( Size const other_index ) const
{
	return ( connections_to_residues_.find( other_index ) != connections_to_residues_.end() );
}

/// @details  Am I polymer bonded to other.seqpos()?
bool
Residue::is_polymer_bonded( Residue const & other ) const
{
	return is_polymer_bonded( other.seqpos() );
}

/// @details  Am I polymer bonded to other_index?
bool
Residue::is_polymer_bonded( Size const other_index ) const
{
	if ( rsd_type_.is_polymer() ) {
		Size const lower_id( rsd_type_.lower_connect_id() ), upper_id( rsd_type_.upper_connect_id() );
		return ( ( lower_id && residue_connection_partner( lower_id ) == other_index ) ||
			( upper_id && residue_connection_partner( upper_id ) == other_index ) );
	} else return false;
}

/// @brief  Returns the atom-index of my atom which is connected to the other residue
/// @details so long as there is only a single connection to other... if there are multiple
/// connections this will fail.  If there are no connections this will fail.
/// This is a convenience function that can fail; be careful!
/// Fails if I'm not bonded to the other residue.
/// @note not well defined if multiple connections to another residue -- need more general function
Size
Residue::connect_atom( Residue const & other ) const
{
	Size const other_seqpos( other.seqpos() );
	if ( is_polymer() && ! nonstandard_polymer_ ) {
		if ( other_seqpos == Size(seqpos_) + 1 && !is_upper_terminus() ) {
			return upper_connect_atom();
		} else if ( other_seqpos == Size(seqpos_) - 1 && !is_lower_terminus() ) {
			return lower_connect_atom();
		}
	}
	if ( connections_to_residues_.find( Size( other.seqpos()) ) != connections_to_residues_.end() ) {
		return rsd_type_.residue_connection( connections_to_residues_.find( Size( other_seqpos) )->second[ 1 ] ).atomno();
	}

	utility_exit_with_message( "Residue::conect_atom: I'm not bonded to that other residue!!");
	return 0;
}

// Get a list of heavy atoms connected to a given atom.
/// @return The atom indices of all heavy atoms bonded to the given atom (by index)
/// @author Labonte <JWLabonte@jhu.edu>
utility::vector1< uint >
Residue::get_adjacent_heavy_atoms( uint const atom_index ) const
{
	utility::vector1< uint > bonded_heavy_atom_indices;

	// Get list of indices of all atoms connected to given connect atom.
	utility::vector1< uint > const bonded_atom_indices( bonded_neighbor( atom_index ) );

	// Search for heavy atoms.  (A residue connection is not an atom.)
	Size const n_indices( bonded_atom_indices.size() );
	for ( uint i( 1 ); i <= n_indices; ++i ) {
		if ( ! atom_is_hydrogen( bonded_atom_indices[ i ] ) ) {
			bonded_heavy_atom_indices.push_back( bonded_atom_indices[ i ] );
		}
	}
	return bonded_heavy_atom_indices;
}

// Scan through the list of atoms connected to a given atom and return the 1st heavy atom found.
/// @return The atom index of the 1st heavy atom next to the given atom (by index) or 0 if no heavy atom is found
/// @remark This method is crucial for determining atoms defining non-standard torsion angles, such as those found
/// across branch connections or in glycosidic linkages.
/// @author Labonte <JWLabonte@jhu.edu>
uint
Residue::first_adjacent_heavy_atom( uint const atom_index ) const
{
	utility::vector1< uint > const atom_indices( get_adjacent_heavy_atoms( atom_index ) );

	if ( atom_indices.empty() ) {
		TR.Warning << "There are no adjacent heavy atoms to atom index " << atom_index << '!' << std::endl;
		return 0;
	}
	return atom_indices[ 1 ];
}

PseudoBondCollectionCOP
Residue::get_pseudobonds_to_residue( Size resid ) const
{
	auto iter( pseudobonds_.find( resid ) );
	if ( iter != pseudobonds_.end() ) {
		return iter->second;
	}
	return nullptr;
}

void
Residue::set_pseudobonds_to_residue( Size resid, PseudoBondCollectionCOP pbs )
{
	pseudobonds_[ resid ] = pbs;
}


/// @details determine how many atoms n the residue and adjacent residues are bonded to the given atom
/// (by default, intraresidue virtual atoms are excluded)
Size
Residue::n_bonded_neighbor_all_res(
	core::Size const atomno,
	bool virt // = false
) const
{
	Size num_neighbors(0);

	chemical::AtomIndices const & intrares_atomnos(bonded_neighbor(atomno));
	for ( Size i = 1; i <= intrares_atomnos.size(); ++i ) {
		if ( virt || ! is_virtual(intrares_atomnos[i]) ) ++num_neighbors;
	}

	Size const num_connections(n_possible_residue_connections());

	for ( Size i = 1; i <= num_connections; ++i ) {
		// this doesn't check the other residue to see if it is connected to a virtual atom
		if ( residue_connect_atom_index(i) == atomno && ! connection_incomplete(i) ) ++num_neighbors;
	}

	return num_neighbors;
}

/// @brief Does this atom depend on the LOWER_CONNECT?
/// @details If recursive is set to true (the default), the dependencies are traced back and the function returns
/// true if ANY atom that this atom depends on is dependent on LOWER_CONNECT.  If false, only the immediate stub1, stub2, and
/// stub3 atoms are checked.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
Residue::atom_depends_on_lower(
	core::Size const atom_index,
	bool const recursive
) const {
	if ( !recursive ) { return icoor(atom_index).depends_on_polymer_lower(); }

	utility::vector1 < core::Size > visited;
	visited.push_back(atom_index);
	return atom_depends_on_lower( atom_index, visited );
}

/// @brief Does this atom depend on the UPPER_CONNECT?
/// @details If recursive is set to true (the default), the dependencies are traced back and the function returns
/// true if ANY atom that this atom depends on is dependent on UPPER_CONNECT.  If false, only the immediate stub1, stub2, and
/// stub3 atoms are checked.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
Residue::atom_depends_on_upper(
	core::Size const atom_index,
	bool const recursive
) const {
	if ( !recursive ) { return icoor(atom_index).depends_on_polymer_upper(); }

	utility::vector1 < core::Size > visited;
	visited.push_back(atom_index);
	return atom_depends_on_upper( atom_index, visited );
}

// fpd bondlength analog to set_chi
//    like set_chi, assumes changes propagate to atomtree
//    keyed off of chi#, so we only allow distances corresponding to chi angles to refine
//    distance corresponds to the distance between atoms 3 and 4 defining the chi
//    chino==0 ==> CA-CB distance, which allows us to refine ALA CB position for example
void
Residue::set_d( int const chino, Real const setting ) {
	int const effchi = (chino==0)? 1 : 0;
	int const baseatom = (chino==0)? 2 : 3;

	AtomIndices const & chi_atoms( rsd_type_.chi_atoms( effchi ) );

	// get the current d
	Real const current_d( ( atom(chi_atoms[baseatom+1]).xyz() - atom(chi_atoms[baseatom]).xyz() ).length() );

	debug_assert( rsd_type_.atom_base( chi_atoms[baseatom] ) == chi_atoms[baseatom] );
	numeric::xyzMatrix< Real > const R(numeric::xyzMatrix<Real>::identity());

	Vector const axis (( atom(chi_atoms[baseatom+1]).xyz() - atom(chi_atoms[baseatom]).xyz() ).normalized());
	Vector const v( (setting-current_d)*axis );

	// apply the transform to all "downstream" atoms
	apply_transform_downstream( chi_atoms[baseatom+1], R, v );

	ASSERT_ONLY(Real const new_d( ( atom(chi_atoms[baseatom+1]).xyz() - atom(chi_atoms[baseatom]).xyz() ).length() );)
		debug_assert( std::abs( new_d - setting ) < 1e-2 );

	update_actcoord();//ek added 4/28/10
}


// fpd bondangle analog to set_chi (see above for details)
void
Residue::set_theta( int const chino, Real const setting ) {
	int const effchi = (chino==0)? 1 : 0;
	int const baseatom = (chino==0)? 2 : 3;

	AtomIndices const & chi_atoms( rsd_type_.chi_atoms( effchi ) );

	// get the current chi angle
	Real const current_theta
		( numeric::angle_degrees( atom( chi_atoms[baseatom-1] ).xyz(), atom( chi_atoms[baseatom] ).xyz(), atom( chi_atoms[baseatom+1] ).xyz() ) );

	Vector const v12( atom(chi_atoms[baseatom]).xyz() - atom(chi_atoms[baseatom-1]).xyz() );
	Vector const v23( atom(chi_atoms[baseatom+1]).xyz() - atom(chi_atoms[baseatom]).xyz() );
	Vector const axis (v12.cross(v23).normalized());

	// debug ordering of chi atoms
	debug_assert( ( rsd_type_.atom_base( chi_atoms[baseatom] ) == chi_atoms[baseatom-1]  ) &&
		( rsd_type_.atom_base( chi_atoms[baseatom+1] ) == chi_atoms[baseatom]  ) );

	numeric::xyzMatrix< Real > const R
		( numeric::rotation_matrix_degrees( axis, - setting + current_theta ) );

	Vector const chi_atom2_xyz( atom( chi_atoms[baseatom] ).xyz() );
	Vector const v( chi_atom2_xyz - R * chi_atom2_xyz );

	// apply the transform to all "downstream" atoms
	apply_transform_downstream( chi_atoms[baseatom], R, v );

	ASSERT_ONLY(Real const new_th(numeric::angle_degrees(
		atom( chi_atoms[baseatom-1] ).xyz(), atom( chi_atoms[baseatom] ).xyz(), atom( chi_atoms[baseatom+1] ).xyz() )); )
		debug_assert( std::abs( basic::subtract_degree_angles( new_th, setting ) ) < 1e-2 );

	update_actcoord();
}

void
Residue::set_tau( Size const nuno, Real const setting )
{
	Size base_id = nuno > nus_.size() ? 3 : 2;

	AtomIndices const & nu_atoms( rsd_type_.nu_atoms( ( nuno > nus_.size() ? nus_.size() : nuno ) ) );

	// get the current tau angle
	Real const current_tau
		( numeric::angle_degrees( atom( nu_atoms[ base_id-1 ] ).xyz(), atom( nu_atoms[ base_id ] ).xyz(), atom( nu_atoms[ base_id+1 ] ).xyz() ) );

	Vector const v12( atom( nu_atoms[ base_id ] ).xyz() - atom( nu_atoms[ base_id-1 ] ).xyz() );
	Vector const v23( atom( nu_atoms[ base_id+1 ] ).xyz() - atom( nu_atoms[ base_id ] ).xyz() );
	Vector const axis (v12.cross(v23).normalized());

	numeric::xyzMatrix< Real > const R
		( numeric::rotation_matrix_degrees( axis, - setting + current_tau ) );

	Vector const nu_atom2_xyz( atom( nu_atoms[ base_id ] ).xyz() );
	Vector const v( nu_atom2_xyz - R * nu_atom2_xyz );

	// apply the transform to all "downstream" atoms
	apply_transform_downstream( nu_atoms[ base_id+1 ], R, v );

	//ASSERT_ONLY(Real const new_tau(numeric::angle_degrees(
	//  atom( nu_atoms[ base_id-1 ] ).xyz(), atom( nu_atoms[ base_id ] ).xyz(), atom( nu_atoms[ base_id+1 ] ).xyz() )); )
	//debug_assert( std::abs( basic::subtract_degree_angles( new_tau, setting ) ) < 1e-2 );

	update_actcoord();
}

void
Residue::set_all_nu( utility::vector1< Real > const & nus, utility::vector1< Real > const & taus )
{
	debug_assert( nus.size() == nus_.size() );
	debug_assert( taus.size() == nus_.size()+1 );

	set_all_ring_nu( 1, nus_.size(), nus, taus );
}

void
Residue::set_all_ring_nu( Size first, Size last, utility::vector1< Real > const & nus, utility::vector1< Real > const & taus )
{

	debug_assert( nus.size() == last-first+1 );
	debug_assert( taus.size() == last-first+2 );

	for ( Size nuno = first; nuno <= last; ++nuno ) {

		nus_[ nuno ] = nus[ nuno ];

		AtomIndices const & nu_atoms( rsd_type_.nu_atoms( nuno ) );

		// get the current nu angle
		Real const current_nu
			( numeric::dihedral_degrees( atom( nu_atoms[1] ).xyz(),
			atom( nu_atoms[2] ).xyz(),
			atom( nu_atoms[3] ).xyz(),
			atom( nu_atoms[4] ).xyz() ) );

		Vector const axis
			(( atom(nu_atoms[3]).xyz() - atom(nu_atoms[2]).xyz() ).normalized());

		numeric::xyzMatrix< Real > const R
			( numeric::rotation_matrix_degrees( axis, nus[ nuno ] - current_nu ) );

		Vector const nu_atom3_xyz( atom( nu_atoms[3] ).xyz() );
		Vector const v( nu_atom3_xyz - R * nu_atom3_xyz );

		// apply the transform to all "downstream" atoms
		apply_transform_downstream( nu_atoms[3], R, v );

		ASSERT_ONLY(Real const new_nu
			( numeric::dihedral_degrees( atom( nu_atoms[1] ).xyz(),
			atom( nu_atoms[2] ).xyz(),
			atom( nu_atoms[3] ).xyz(),
			atom( nu_atoms[4] ).xyz() ) );)
			debug_assert( std::abs( basic::subtract_degree_angles( new_nu, nus[ nuno ] ) ) <
			1e-2 );

		update_actcoord();//ek added 4/28/10

	}

	for ( Size nuno = 1; nuno <= nus.size(); ++nuno ) {
		set_tau( nuno, taus[ nuno ] );
	}
	set_tau( nus.size()+1, taus[ nus.size()+1 ] );
}

/////////////////////////////////////////////////////////////////////////////
/// @details this assumes that change propagates according to the information from
/// atom_base array, not from atom tree. So be sure not to get into an
/// endless loop.
void
Residue::set_chi( int const chino, Real const setting )
{
	chi_[ chino ] = setting;

	AtomIndices const & chi_atoms( rsd_type_.chi_atoms( chino ) );

	// get the current chi angle
	Real const current_chi
		( numeric::dihedral_degrees( atom( chi_atoms[1] ).xyz(),
		atom( chi_atoms[2] ).xyz(),
		atom( chi_atoms[3] ).xyz(),
		atom( chi_atoms[4] ).xyz() ) );

	Vector const axis
		(( atom(chi_atoms[3]).xyz() - atom(chi_atoms[2]).xyz() ).normalized());
	// debug ordering of chi atoms
	debug_assert( ( rsd_type_.atom_base( chi_atoms[3] ) == chi_atoms[2]  ) &&
		( rsd_type_.atom_base( chi_atoms[4] ) == chi_atoms[3]  ) );

	numeric::xyzMatrix< Real > const R
		( numeric::rotation_matrix_degrees( axis, setting - current_chi ) );

	Vector const chi_atom3_xyz( atom( chi_atoms[3] ).xyz() );
	Vector const v( chi_atom3_xyz - R * chi_atom3_xyz );

	// apply the transform to all "downstream" atoms
	apply_transform_downstream( chi_atoms[3], R, v );


	ASSERT_ONLY(Real const new_chi
		( numeric::dihedral_degrees( atom( chi_atoms[1] ).xyz(),
		atom( chi_atoms[2] ).xyz(),
		atom( chi_atoms[3] ).xyz(),
		atom( chi_atoms[4] ).xyz() ) );)
		debug_assert( std::abs( basic::subtract_degree_angles( new_chi, setting ) ) <
		1e-2 );

	update_actcoord();//ek added 4/28/10
}


void
Residue::set_all_chi( utility::vector1< Real > const & chis )
{
	// This works for now, but there's probably a faster implementation which only runs the coordinate update once.
	for ( Size i=1; i<= nchi(); ++i ) {
		set_chi( i, chis[i] );
	}
}


/////////////////////////////////////////////////////////////////////////////
/// @details xyz --> R * xyz + v \n
/// this uses information from atom_base array to transform all the downstream atoms
/// along the side chain recursively. it assumes that the atom_base array will not get
/// us into any infinite loops!
///
/// @note this is not for general atom tree folding. only used in set_chi in which
/// changes for a chi angle is fast propagated within one residue and not to invoke
/// folding the whole atom tree.
void
Residue::apply_transform_downstream(
	int const atomno,
	numeric::xyzMatrix< Real > const & R,
	Vector const & v
)
{
	// transform my coordinates:: xyz -> R * xyz + v
	//
	atom( atomno ).xyz( R * atom( atomno ).xyz() + v );

	// now apply recursively to my downstream nbrs:
	AtomIndices const & nbrs( rsd_type_.bonded_neighbor( atomno ) );
	int const my_atom_base( rsd_type_.atom_base( atomno ) );
	for ( Size i=1; i<= nbrs.size(); ++i ) {

		int const nbr( nbrs[i] );
		int const nbr_base( rsd_type_.atom_base( nbr ) );

		if ( nbr_base == atomno ) {
			if ( my_atom_base != nbr ) {
				apply_transform_downstream( nbr, R, v );
			} else {
				if ( atomno!=1 && atomno!=2 ) TR.Warning << "DANGER: almost got stuck in infinite loop!  Atom " << atomno << " is apparently a parent AND child of atom " << nbr << "." << std::endl;
				//Note: atom 1 and atom 2 have a weird relationship, and that leads to bogus warning messages here.
			}
		}
	}
}

void
Residue::apply_transform_Rx_plus_v(
	numeric::xyzMatrix< Real > R,
	Vector v
) {
	for ( Size atom_idx = 1; atom_idx <= type().natoms(); ++atom_idx ) {
		set_xyz( atom_idx, R * xyz(atom_idx) + v );
	}
}


void
Residue::determine_nonstandard_polymer_status()
{
	if ( is_polymer() ) {
		if ( ! is_upper_terminus() &&
				( type().upper_connect_id() == 0 ||
				connect_map_[ type().upper_connect_id() ].incomplete() ||
				connect_map_[ type().upper_connect_id() ].resid() != seqpos() + Size( 1 )) ) {
			nonstandard_polymer_ = true;
			return;
		}
		if ( ! is_lower_terminus() &&
				( type().lower_connect_id() == 0 ||
				connect_map_[ type().lower_connect_id() ].incomplete() ||
				connect_map_[ type().lower_connect_id() ].resid() != seqpos() - Size( 1 )) ) {
			nonstandard_polymer_ = true;
			return;
		}
	}
	nonstandard_polymer_ = false;
}

/// AMW TODO: somehow SWA spends literally 3.7% of its time calling this function.
/// @note A misnomer; this should really be called "is_virtual_atom()". ~Labonte
bool
Residue::is_virtual( Size const & atomno ) const
{
	return rsd_type_.atom_type( atomno ).is_virtual();
}

/// @note A misnomer; this should really be called "is_repulsive_atom()". ~Labonte
bool
Residue::is_repulsive( Size const & atomno ) const
{
	return rsd_type_.atom_type( atomno ).is_repulsive();
}


/// @details Return a COP to the data cache
/// @note Might be a null pointer if the cache has not been initialized
basic::datacache::BasicDataCacheCOP
Residue::data_ptr() const
{
	return data_cache_;
}

/// @details Return an OP to the datacache
/// @note Will create new one if not already initialized
basic::datacache::BasicDataCacheOP
Residue::nonconst_data_ptr()
{
	if ( data_cache_ == nullptr ) {
		data_cache_ = basic::datacache::BasicDataCacheOP( new basic::datacache::BasicDataCache( residue_datacache::n_cacheable_types ) );
	}
	return data_cache_;
}

#ifdef    SERIALIZATION

template < class Archive >
void
Residue::save( Archive & arc ) const
{
	// with the ResidueTypeCOP, find out if the ResidueType is held in one of the (global)
	// ResidueTypeSets and if so, then Residue will serialize the name of the RTS and the
	// name of the ResidueType -- otherwise, it will serialize the RT itself.
	core::chemical::serialize_residue_type( arc, rsd_type_ptr_ );
	// EXEMPT rsd_type_

	arc( atoms_, orbitals_ );
	arc( seqpos_, chain_ );
	arc( mirrored_relative_to_type_ );
	arc( chi_, nus_, mainchain_torsions_, actcoord_ );
	arc( data_cache_ );
	arc( nonstandard_polymer_, connect_map_ );
	arc( connections_to_residues_ );
	arc( pseudobonds_ );
}

template < class Archive >
void
Residue::load_and_construct(
	Archive & arc,
	cereal::construct< Residue > & construct
)
{
	ResidueTypeCOP restype;
	core::chemical::deserialize_residue_type( arc, restype );
	construct( restype, true );
	// EXEMPT rsd_type_ptr_ rsd_type_

	arc( construct->atoms_, construct->orbitals_ );
	arc( construct->seqpos_, construct->chain_ );
	arc( construct->mirrored_relative_to_type_ );
	arc( construct->chi_, construct->nus_, construct->mainchain_torsions_, construct->actcoord_ );
	arc( construct->data_cache_ );
	arc( construct->nonstandard_polymer_, construct->connect_map_ );
	arc( construct->connections_to_residues_ );

	// deserialization of constant owning pointers requires first deserializing
	// into locally-declared non-constant owning pointers and then assigning
	// from the non-const versions into const versions.  Sadly, this means
	// iterating across the non-const map and inserting elements into the const-map
	// one at a time.
	// too bad the following code doesn't work! construct->pseudobonds_ = nonconst_pbs;
	std::map< Size, PseudoBondCollectionOP > nonconst_pbs;
	arc( nonconst_pbs );
	for ( std::map< Size, PseudoBondCollectionOP >::const_iterator
			ncpb = nonconst_pbs.begin(), ncpb_end = nonconst_pbs.end();
			ncpb != ncpb_end; ++ncpb ) {
		construct->pseudobonds_[ ncpb->first ] = ncpb->second;
	}
}

SAVE_AND_LOAD_AND_CONSTRUCT_SERIALIZABLE( Residue );

#endif // SERIALIZATION

void Residue::assign_orbitals() {
	// Assign orbitals.
	for ( core::Size const atom_with_orbitals : rsd_type_.atoms_with_orb_index() ) {
		utility::vector1<core::Size> const & orbital_indices(rsd_type_.bonded_orbitals(atom_with_orbitals));
		for ( core::Size const orbital_index : orbital_indices ) {
			Vector orb_xyz(this->build_orbital_xyz(orbital_index));
			core::Size type = rsd_type_.orbital(orbital_index).orbital_type_index();
			orbitals_.push_back(orbitals::OrbitalXYZCoords(orb_xyz, type));
		}
	}
}

/// @brief Does this atom depend on the LOWER_CONNECT?  Checks recursively.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
Residue::atom_depends_on_lower(
	core::Size const atom_index,
	utility::vector1< core::Size > visited_atoms
) const {
	if ( icoor(atom_index).depends_on_polymer_lower() ) return true;

	visited_atoms.push_back( atom_index );

	core::Size const stub1_index( icoor( atom_index ).stub_atom1().atomno() );
	bool const stub1( stub1_index == 0 || in_list( stub1_index, visited_atoms ) ? false : atom_depends_on_lower( stub1_index, visited_atoms )  );
	if ( stub1 ) return true;

	core::Size const stub2_index( icoor( atom_index ).stub_atom2().atomno() );
	bool const stub2( stub2_index == 0 || in_list( stub2_index, visited_atoms )  ? false : atom_depends_on_lower( stub2_index, visited_atoms )  );
	if ( stub2 ) return true;

	core::Size const stub3_index( icoor( atom_index ).stub_atom3().atomno() );
	bool const stub3( stub3_index == 0 || in_list( stub3_index, visited_atoms )  ? false : atom_depends_on_lower( stub3_index, visited_atoms )  );
	return stub3;
}

/// @brief Does this atom depend on the UPPER_CONNECT?  Checks recursively.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
Residue::atom_depends_on_upper(
	core::Size const atom_index,
	utility::vector1< core::Size > visited_atoms
) const {
	if ( icoor(atom_index).depends_on_polymer_upper() ) return true;

	visited_atoms.push_back( atom_index );

	core::Size const stub1_index( icoor( atom_index ).stub_atom1().atomno() );
	bool const stub1( stub1_index == 0 || in_list( stub1_index, visited_atoms ) ? false : atom_depends_on_upper( stub1_index, visited_atoms )  );
	if ( stub1 ) return true;

	core::Size const stub2_index( icoor( atom_index ).stub_atom2().atomno() );
	bool const stub2( stub2_index == 0 || in_list( stub2_index, visited_atoms )  ? false : atom_depends_on_upper( stub2_index, visited_atoms )  );
	if ( stub2 ) return true;

	core::Size const stub3_index( icoor( atom_index ).stub_atom3().atomno() );
	bool const stub3( stub3_index == 0 || in_list( stub3_index, visited_atoms )  ? false : atom_depends_on_upper( stub3_index, visited_atoms )  );
	return stub3;
}

/// @brief Is a value in a list?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
Residue::in_list(
	core::Size const val,
	utility::vector1< core::Size > const &vect
) const {
	for ( core::Size i=1, imax=vect.size(); i<=imax; ++i ) { if ( vect[i] == val ) return true; }
	return false;
}

////////////////////////////////////////////////////////////////////////////////
//ja
std::ostream & operator << ( std::ostream & os, Residue const & res )
{
	res.show(os);
	return os;
}

} // conformation
} // core
