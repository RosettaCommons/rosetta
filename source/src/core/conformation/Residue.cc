// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   Residue.cc
/// @brief
/// @author Phil Bradley

// Unit header
#include <core/conformation/Residue.hh>

// Package headers
#include <core/conformation/PseudoBond.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>

// Project headers
#include <core/kinematics/Stub.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#ifdef USEBOOSTSERIALIZE
#include <core/pack/dunbrack/SingleLigandRotamerLibrary.hh>
#endif

// Basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>

// Utility Headers
#include <utility/assert.hh>
#include <utility/PyAssert.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <iostream>

// Boost headers
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH


namespace core {
namespace conformation {

static basic::Tracer TR("core.conformation.Residue");

/// @details Constructor from residue type; sets coords to ideal values
/// create a residue of type residue_type_in.
/// @note Dummmy arg to prevent secret type conversions from ResidueType to Residue
Residue::Residue( ResidueType const & rsd_type_in, bool const /*dummy_arg*/ ):
	utility::pointer::ReferenceCount(),
	rsd_type_( rsd_type_in ),
	seqpos_( 0 ),
	chain_( 0 ),
	chi_( rsd_type_.nchi(), 0.0 ), // uninit
	mainchain_torsions_( rsd_type_.mainchain_atoms().size(), 0.0 ),
	actcoord_( 0.0 ),
	nonstandard_polymer_( false ),
	connect_map_( rsd_type_in.n_residue_connections() )
{

	for ( Size i=1; i<= rsd_type_.natoms(); ++i ) {
		atoms_.push_back( Atom( rsd_type_.atom(i).ideal_xyz(), rsd_type_.atom(i).atom_type_index(),  rsd_type_.atom(i).mm_atom_type_index() ) );
		//std::cout << this->atom_name(i) << std::endl;
	}


	foreach(core::Size atom_with_orbitals, rsd_type_.atoms_with_orb_index()){
		utility::vector1<core::Size> const & orbital_indices(rsd_type_.bonded_orbitals(atom_with_orbitals));
		foreach(core::Size orbital_index, orbital_indices){
			Vector orb_xyz(this->build_orbital_xyz(orbital_index));
			core::Size type = rsd_type_.orbital(orbital_index).orbital_type_index();
			orbitals_.push_back(orbitals::OrbitalXYZCoords(orb_xyz, type));
		}
	}



}

/// @details Create a residue/rotamer of type rsd_type_in placed at the position occupied by current_rsd
/// Used primarily in rotamer building. The newly created Residue has the same sequence postion, chain id
/// and mainchain torsion angles as current_rsd. It has a ResidueType as defined by rsd_type_in. Its sidechain
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
	rsd_type_( rsd_type_in ),
	seqpos_( current_rsd.seqpos() ),
	chain_( current_rsd.chain() ),
	chi_( rsd_type_.nchi(), 0.0 ), // uninit
	mainchain_torsions_( current_rsd.mainchain_torsions() ),
	actcoord_( 0.0 ),
	nonstandard_polymer_( current_rsd.nonstandard_polymer_ ),
	connect_map_( current_rsd.connect_map_ ),
	connections_to_residues_( current_rsd.connections_to_residues_ ),
	pseudobonds_( current_rsd.pseudobonds_ )
{



	for ( Size i=1; i<= rsd_type_.natoms(); ++i ) {
		atoms_.push_back( Atom( rsd_type_.atom(i).ideal_xyz(), rsd_type_.atom(i).atom_type_index(), rsd_type_.atom(i).mm_atom_type_index() ));
	}

	assert( current_rsd.mainchain_torsions().size() == rsd_type_.mainchain_atoms().size() );

	// now orient
	place( current_rsd, conformation, preserve_c_beta );

	// Assumption: if two residue types have the same number of residue connections,
	// then their residue connections are "the same" residue connections.
	// This assumption works perfectly for all amino acids, except CYD.
	//
	// THIS REALLY NEEDS TO BE FIXED AT SOME POINT
	//

	if ( rsd_type_in.n_residue_connections() != current_rsd.type().n_residue_connections() ) {
		if ( ! current_rsd.pseudobonds_.empty() ) {
			std::cerr << "Unable to handle change in the number of residue connections in the presence of pseudobonds!" <<
				std::endl;
			utility_exit();
		}

		copy_residue_connections( current_rsd );

	}

	// This seems a little silly, but the update of chi's doesn't seem to occur automatically in
	// any of the functions above.
	for ( Size chino = 1; chino <= rsd_type_.nchi(); chino++ ) {
		AtomIndices const & chi_atoms( rsd_type_.chi_atoms( chino ) );

		// get the current chi angle
		Real const current_chi
		(
				numeric::dihedral_degrees(
						atom( chi_atoms[1] ).xyz(),
						atom( chi_atoms[2] ).xyz(),
						atom( chi_atoms[3] ).xyz(),
						atom( chi_atoms[4] ).xyz()
				)
		);
		chi_[ chino ] = current_chi;
	}

	foreach(core::Size atom_with_orbitals, rsd_type_.atoms_with_orb_index()){
		utility::vector1<core::Size> const & orbital_indices(rsd_type_.bonded_orbitals(atom_with_orbitals));
		foreach(core::Size orbital_index, orbital_indices){
			Vector orb_xyz(this->build_orbital_xyz(orbital_index));
			core::Size type = rsd_type_.orbital(orbital_index).orbital_type_index();
			orbitals_.push_back(orbitals::OrbitalXYZCoords(orb_xyz, type));
		}
	}

}


Residue::Residue( Residue const & src )
:
	utility::pointer::ReferenceCount(),
	rsd_type_(src.rsd_type_),
	atoms_(src.atoms_),
	orbitals_(src.orbitals_),
	seqpos_(src.seqpos_),
	chain_(src.chain_),
	chi_(src.chi_),
	mainchain_torsions_(src.mainchain_torsions_),
	actcoord_(src.actcoord_),
	nonstandard_polymer_(src.nonstandard_polymer_),
	connect_map_(src.connect_map_),
	connections_to_residues_(src.connections_to_residues_),
	pseudobonds_(src.pseudobonds_)
{

}

Residue::~Residue() {}


///@details make a copy of this residue( allocate actual memory for it )
ResidueOP
Residue::clone() const
{
	return new Residue( *this );
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
	assert(rsd_type_.is_carbohydrate());
	PyAssert(rsd_type_.is_carbohydrate(), "Residue::carbohydrate_info(): This residue is not a carbohydrate!");

	return rsd_type_.carbohydrate_info();
}


bool Residue::connections_match( Residue const & other ) const
{
	if ( connect_map_.size() != other.connect_map_.size() ) return false;
	//if ( connections_to_residues_.size() != other.connections_to_residues_.size()) return false; // duplicate data
	if ( pseudobonds_.size() != other.pseudobonds_.size() ) return false;

	for ( Size ii = 1; ii <= connect_map_.size(); ++ii ) {
		if ( connect_map_[ ii ] != other.connect_map_[ ii ] ) return false;
	}
	for ( std::map< Size, PseudoBondCollectionCOP >::const_iterator
			iter = pseudobonds_.begin(), iter_end = pseudobonds_.end(),
			other_iter_end = other.pseudobonds_.end();
			iter != iter_end; ++iter ) {
		std::map< Size, PseudoBondCollectionCOP >::const_iterator other_iter = other.pseudobonds_.find( iter->first );
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
	if (chi_.size() != other_chi.size() || rsd_type_.aa() != other.aa() || rsd_type_.name3() != other.name3() ){
		return false;
	}
	else {
		for (Size i = 1; i<= chi_.size(); i++){
				if ( std::abs( this_chi[i] - other_chi[i]) >= 5){
					match = false;
				}
		}
	}
	return match;
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
	//

	if ( type().n_residue_connections() == src_rsd.type().n_residue_connections() ) {

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

		connect_map_.resize( type().n_residue_connections() );

		// Find correspondence between src_rsd's connection atoms and atoms on *this.
		for ( Size ii = 1; ii <= src_rsd.type().n_residue_connections(); ++ii ) {
			Size const ii_connatom = src_rsd.type().residue_connection( ii ).atomno();
			if ( has( src_rsd.atom_name( ii_connatom ) )) {

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
							TR.Debug << "WARNING: ResConnID info stored on residue " << src_rsd.connect_map( ii ).resid();
							TR.Debug << " is now out of date!" << std::endl;
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
						///	<< " to " << name() << " at position " << seqpos() << std::endl;
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
	//	std::cout << " CENTER " << atom_name( center ) << "   NBR1 " << atom_name( nbr1 ) << "    NBR2 " << atom_name( nbr2 ) << std::endl;

	assert( center && nbr1 && nbr2 );
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
	if( atom_pairs.size() != 3 ){
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


/// @details place/orient "this" Residue onto "src" Residue by backbone superimposition
///	Since rotamer is represented by Residue in mini now, this function is mainly used to place a rotamer
///	onto the backbone of "src" residue. Meanwhile, it can also be used to add sidechains to one pose/conformation
///	from another pose/conformation.\n
///	current logic: find backbone atom with bonded neighbors in sidechain,
///	and which is the base_atom of those neighbors. Take that backbone atom
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
	Size first_scatom( rsd_type_.first_sidechain_atom() );
	if ( first_scatom >= 1 && first_scatom <= rsd_type_.nheavyatoms() ) {
		// not all backbone -- need to do some orienting
		orient_onto_residue( src );
	} // does the residue have any sidechain atoms?

	// now copy the backbone atoms
	//
	utility::vector1< bool > missing( natoms(), false );
	bool any_missing( false );
	for ( Size i=1; i<= natoms(); ++i ) {

		//The O2' is a special case for RNA, because it is a "sidechain" atom that
		// branches off the backbone separately from the base. This could be
		// coded more robustly by modifying orient_onto_residue() properly.

		if ( !rsd_type_.atom_is_backbone(i) && !(rsd_type_.atom_name(i) == " O2'") && !(rsd_type_.atom_name(i) == "HO2'") ) continue;
		if ( src.has( rsd_type_.atom_name(i) ) ) {
			atoms()[i].xyz( src.atom( src.atom_index( rsd_type_.atom_name(i) ) ).xyz() );
		} else {
			missing[i] = true;
			any_missing = true;
		}
	}

	if ( any_missing ) fill_missing_atoms( missing, conformation );

	if ( preserve_c_beta ) {
		// after superposition, adjust the sidechain atoms by aligning the CA-CB bond
		//		std::cout << "preserving c-beta... " << std::endl;
		std::string root("CA"), mobile_new("CB"), mobile_src("CB");

		if ( is_RNA() ){
			root = " C1'";
			mobile_new = atom_name( chi_atoms( 1 )[ 3 ] ); //First atom in base...
			mobile_src = src.atom_name( src.chi_atoms( 1 )[ 3 ] ); //First atom in base...
		}

		// only try this when both the src and the new residue types contain both atom types
		if ( type().has( root ) && type().has( mobile_new ) &&
				 src.type().has( root ) && src.type().has( mobile_src ) ) {
			assert( xyz( root ) == src.xyz( root ) ); // roots should be aligned by now
			// common 'pseudoatom' vector, perpendicular to the plane defined by the two bonds
			Vector const pseudoatom(
				cross( xyz( mobile_new ) - xyz( root ), src.xyz( mobile_src ) - src.xyz( root ) ) + xyz( root )
			);

			if (pseudoatom==xyz(root)) return;

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


/////////////////////////////////////////////////////////////////////////////
/// @details
/// this uses ideal internal coords to build any missing atom from its three
/// stub atoms. If any of the stub atoms are missing, build them first.
/// Unable to build a missing atom whose stub atoms are from non-existing
/// polymer connection and its input bogus value will not be changed.
void
Residue::fill_missing_atoms(
	utility::vector1< bool > missing, // make local copy
	Conformation const & conformation
)
{
	bool still_missing( true );
	while ( still_missing ) {
		still_missing = false;
		for ( Size i=1; i<= natoms(); ++i ) {
			if ( missing[i] ) {
				chemical::AtomICoor const & ic( icoor(i) );
				if ( (seqpos_ == 1                   && ic.depends_on_polymer_lower()) ||
					(Size(seqpos_) == conformation.size() && ic.depends_on_polymer_upper()) ) {
					missing[i] = false;
					TR.Warning << "[ WARNING ] missing an atom: " << seqpos_ << " " << atom_name(i) << " that depends on a nonexistent polymer connection! "
						<< std::endl <<  " --> generating it using idealized coordinates." << std::endl;
					set_xyz( i, ic.build(*this));
					continue;
				}
				still_missing = true;
				// check to see if any of our stub atoms are missing:
				bool stub_atoms_missing( false );
				for ( Size j=1; j<= 3; ++j ) {
					chemical::ICoorAtomID const & id( ic.stub_atom(j) );
					if ( id.type() == chemical::ICoorAtomID::INTERNAL && missing[ id.atomno() ] ) {
						stub_atoms_missing = true;
						if ( id.atomno() == i ) {
							TR.Error << "[ ERROR ] missing atom " << i << " (" << atom_name(i) << ") in " << type().name() << " is its own stub" << std::endl;
							utility_exit_with_message("Endless loop in fill_missing_atoms()");
						}
						break;
					}
				}

				if ( !stub_atoms_missing ) {
					// no stub atoms missing: build our ideal coordinates
					missing[i] = false;
					//std::cout << "Residue::fill_missing_atoms: rebuild backbone atom: " << name() << ' ' <<
					//	atom_name(i) << std::endl;
					set_xyz( i, build_atom_ideal( i, conformation ) );
				}
			}
		}
	}
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
	Size const num_connections(n_residue_connections());

	for (Size i = 1; i <= num_connections; ++i) {
		if (residue_connect_atom_index(i) == atomno && connection_incomplete(i)) return true;
	}

	return false;
}


bool
Residue::connection_incomplete( Size resconnid ) const
{
	return connect_map_[ resconnid ].incomplete();
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
// 	utility::vector1< Size > newlist;
// 	if (  connections_to_residues_.find( otherres ) != connections_to_residues_.end() ) {
// 		newlist = connections_to_residues_[ otherres ];
// 	}
// 	if ( newlist.size() != 0 ) {
// 		for ( Size ii = 1; ii <= newlist.size(); ++ii ) {
// 			if ( newlist[ ii ] == resconn_index  ) {
// 				//std::cout << "Setting residue connection partner on residue " << seqpos_ << " to residue " << otherres << " twice!" << std::endl;
// 				break;
// 			}
// 			else if ( ii == newlist.size() ) {
// 				newlist.push_back( resconn_index );
// 				connections_to_residues_[ otherres ] = newlist;
// 				break;
// 			}
// 		}
// 	} else {
// 		newlist.push_back( resconn_index );
// 		connections_to_residues_[ otherres ] = newlist;
// 	}
	determine_nonstandard_polymer_status();
}


/// @details  Private function to keep the connections_to_residues_ array up to date
/// @note  This could be made faster -- connections_to_residues_ is a std::map< > so operator[] calls are slow
void
Residue::update_connections_to_residues()
{
	connections_to_residues_.clear();
	for ( Size i=1, i_end = n_residue_connections(); i<= i_end; ++i ) {
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
		if( new_resid == 0 ) connect_map_[i].connid( 0 );


// 		if ( new_resid ) {
// 			connections_to_residues_[ new_resid ] = connections_to_residues_copy[ old_resid ];
// 		}
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

Distance
Residue::connection_distance(
	conformation::Conformation const & conf,
	Size const resconn_index,
	Vector const matchpoint
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
// 	if ( is_polymer() && ! nonstandard_polymer_ ) {
// 		if ( polymeric_sequence_distance( other ) == 1 ) {
// 			// confirm that termini status is consistent with sequence_distance, which depends on chain
// 			assert( ( other.seqpos() == seqpos() + 1 && !is_upper_terminus() ) ||
// 				( other.seqpos() == seqpos() - 1 && !is_lower_terminus() ) );
// 			return true;
// 		} else if ( rsd_type_.n_non_polymeric_residue_connections() == 0 ) {
// 			return false; // generic case
// 		}
// 	}
// 	return ( connections_to_residues_.find( Size(other.seqpos()) ) != connections_to_residues_.end() );
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

PseudoBondCollectionCOP
Residue::get_pseudobonds_to_residue( Size resid ) const
{
	std::map< Size, PseudoBondCollectionCOP >::const_iterator iter( pseudobonds_.find( resid ) );
	if ( iter != pseudobonds_.end() ) {
		return iter->second;
	}
	return 0;
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
	for (Size i = 1; i <= intrares_atomnos.size(); ++i) {
		if (virt || ! is_virtual(intrares_atomnos[i]) ) ++num_neighbors;
	}

	Size const num_connections(n_residue_connections());

	for (Size i = 1; i <= num_connections; ++i) {
		// this doesn't check the other residue to see if it is connected to a virtual atom
		if (residue_connect_atom_index(i) == atomno && ! connection_incomplete(i)) ++num_neighbors;
	}

	return num_neighbors;
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

	assert( rsd_type_.atom_base( chi_atoms[baseatom] ) == chi_atoms[baseatom] );
	numeric::xyzMatrix< Real > const R(numeric::xyzMatrix<Real>::identity());

	Vector const axis (( atom(chi_atoms[baseatom+1]).xyz() - atom(chi_atoms[baseatom]).xyz() ).normalized());
	Vector const v( (setting-current_d)*axis );

	// apply the transform to all "downstream" atoms
	apply_transform_downstream( chi_atoms[baseatom+1], R, v );

	ASSERT_ONLY(Real const new_d( ( atom(chi_atoms[baseatom+1]).xyz() - atom(chi_atoms[baseatom]).xyz() ).length() );)
	assert( std::abs( new_d - setting ) < 1e-2 );

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
	Vector const axis	(v12.cross(v23).normalized());

	// debug ordering of chi atoms
	assert( ( rsd_type_.atom_base( chi_atoms[baseatom] ) == chi_atoms[baseatom-1]  ) &&
		( rsd_type_.atom_base( chi_atoms[baseatom+1] ) == chi_atoms[baseatom]  ) );

	numeric::xyzMatrix< Real > const R
		( numeric::rotation_matrix_degrees( axis, - setting + current_theta ) );

	Vector const chi_atom2_xyz( atom( chi_atoms[baseatom] ).xyz() );
	Vector const v( chi_atom2_xyz - R * chi_atom2_xyz );

	// apply the transform to all "downstream" atoms
	apply_transform_downstream( chi_atoms[baseatom], R, v );

	ASSERT_ONLY(Real const new_th(numeric::angle_degrees(
	              atom( chi_atoms[baseatom-1] ).xyz(), atom( chi_atoms[baseatom] ).xyz(), atom( chi_atoms[baseatom+1] ).xyz() )); )
	assert( std::abs( basic::subtract_degree_angles( new_th, setting ) ) < 1e-2 );

	update_actcoord();
}


/////////////////////////////////////////////////////////////////////////////
/// @details this assumes that change propagates according to the information from
/// atom_base array, not from atom tree. So be sure not to get into an
/// endless loop.
void
Residue::set_chi( int const chino, Real const setting )
{

//#ifdef NDEBUG
//	bool const debug( false );
//#else
//	bool const debug( true );
//#endif

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
	assert( ( rsd_type_.atom_base( chi_atoms[3] ) == chi_atoms[2]  ) &&
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
	assert( std::abs( basic::subtract_degree_angles( new_chi, setting ) ) <
			1e-2 );

	update_actcoord();//ek added 4/28/10
}


void
Residue::set_all_chi( utility::vector1< Real > const & chis )
{
	// This works for now, but there's probably a faster implementation which only runs the coordinate update once.
	for(Size i=1; i<= nchi(); i++) {
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
				TR.Warning << "DANGER: almost got stuck in infinite loop!" << std::endl;
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


/// @note A misnomer; this should really be called "is_virtual_atom()". ~Labonte
bool
Residue::is_virtual( Size const & atomno ) const
{
	return rsd_type_.atom_type( atomno ).is_virtual();
}


////////////////////////////////////////////////////////////////////////////////
//ja
std::ostream & operator << ( std::ostream & os, Residue const & res )
{
	os << res.name() << ' ' << res.seqpos() << ": \n";
	for ( Size j=1; j<=res.natoms(); ++j ) {
		Atom const & atom ( res.atom(j) );
		os << res.atom_name(j) << ": ";
		os << atom.xyz().x() << ' ' << atom.xyz().y() << ' '<< atom.xyz().z();
		if (res.is_virtual(j)) {
			os << " (virtual)";
		}
		os << std::endl;
	}
	if (res.is_carbohydrate()) {
		os << std::endl;
		res.carbohydrate_info()->show(os);
	}
	return os;
}

std::ostream & operator << ( std::ostream & os, Atom const & atom )
{
	os << "Atom type:" << atom.type() << " xyz:" << atom.xyz().x() << ' ' << atom.xyz().y() << ' ' << atom.xyz().z();
	return os;
}


#ifdef USEBOOSTSERIALIZE
// this function takes the old res, clones it rotamer library by applying it to the new res
// and adds it to the rotamerlibrary singleton
void add_cloned_ligand_rotamer_library( core::chemical::ResidueType & new_res, core::chemical::ResidueType const & base_res ) {
	using namespace core::pack::dunbrack;

	SingleLigandRotamerLibraryOP new_lrots = new SingleLigandRotamerLibrary;
	SingleLigandRotamerLibraryCAP old_lrots(
		static_cast< SingleLigandRotamerLibrary const * >
		( RotamerLibrary::get_instance().get_rsd_library( base_res )() ));
	if( old_lrots != 0 ) {
		utility::vector1< ResidueOP > new_rotamers;
		utility::vector1< ResidueOP > const old_rotamers = old_lrots->get_rotamers();
		for( utility::vector1< ResidueOP>::const_iterator oldrot_it = old_rotamers.begin(); oldrot_it != old_rotamers.end(); ++oldrot_it){
			ResidueOP new_rot_res = new Residue( new_res, true);
			for( core::Size at_ct = 1; at_ct <= new_rot_res->natoms(); at_ct++){
				if( !(*oldrot_it)->has( new_rot_res->atom_name( at_ct ) ) ){
					std::cerr << "Unexpected ERROR: when regenerating ligand rotamer library (for covalent constraints), one atom wasn't found in a template rotamer." << std::endl;
					utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
				}
				else{
					new_rot_res->set_xyz( at_ct, (*oldrot_it)->xyz( new_rot_res->atom_name( at_ct ) ) );
				}
			}
			new_rot_res->chi( (*oldrot_it)->chi() );
			new_rotamers.push_back( new_rot_res );
		}
		new_lrots->set_reference_energy( old_lrots->get_reference_energy() );
		new_lrots->set_rotamers( new_rotamers );
	} // no fallback if there isnt a reference rotamer library
	RotamerLibrary::get_instance().add_residue_library( new_res, new_lrots );
}
#endif

} // conformation
} // core

