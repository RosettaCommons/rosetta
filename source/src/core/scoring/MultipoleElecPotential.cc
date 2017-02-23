// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/MultipoleElecPotential.cc
/// @brief
/// @author Jim Havranek

// Project headers
#include <core/scoring/MultipoleElecPotential.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/prof.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>

#include <numeric/constants.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.io.hh>
#include <numeric/xyz.functions.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <boost/lexical_cast.hpp>

#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <map>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.MultipoleElecPotential" );

typedef numeric::xyzVector< core::Real > Vector;
typedef numeric::xyzMatrix< core::Real > Matrix;

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

MultipoleElecResidueInfo::~MultipoleElecResidueInfo() = default;

bool
ShouldItCount(
	conformation::Residue const & rsd,
	Size const & atm
) {
	if ( rsd.seqpos() == 2 ) return false;

	if ( rsd.atom_name( atm ) == " CB " ) {
		return true;
	} else {
		return false;
	}
}

///
void
MultipoleElecResidueInfo::initialize( conformation::Residue const & rsd )
{
	Size const natoms( rsd.natoms() );
	type_.resize( natoms );
	coord_frame_ref_.resize( natoms );
	monopole_.resize( natoms );
	rKirkwood_.resize( natoms );
	dipole_.resize( natoms );
	induced_dipole_.resize( natoms );
	stored_induced_dipole_.resize( natoms );
	induced_rf_dipole_.resize( natoms );
	stored_induced_rf_dipole_.resize( natoms );
	Efield_fixed_.resize( natoms );
	Efield_induced_.resize( natoms );
	Efield_rf_fixed_.resize( natoms );
	Efield_rf_induced_.resize( natoms );
	quadrupole_.resize( natoms );
	local_coord_matrix_.resize( natoms );
	my_group_.resize( natoms );
	my_local_coord_frame_.resize( natoms );
	mp_param_.resize( natoms );
	rosetta_res_type_ = rsd.name();

	// Assign each atom their Amoeba type, which is the atom type
	// assigned by Amoeba.  It does not have a one-to-one mapping with
	// Rosetta atom types.

	// Look up values by Amoeba type
	for ( Size i=1; i<= natoms; ++i ) {
		monopole_[ i ] = 0.0;
		rKirkwood_[ i ] = 0.0;
		dipole_[ i ] = 0.0;
		induced_dipole_[ i ] = 0.0;
		stored_induced_dipole_[ i ] = 0.0;
		induced_rf_dipole_[ i ] = 0.0;
		stored_induced_rf_dipole_[ i ] = 0.0;
		Efield_fixed_[ i ] = 0.0;
		Efield_induced_[ i ] = 0.0;
		Efield_rf_fixed_[ i ] = 0.0;
		Efield_rf_induced_[ i ] = 0.0;
		quadrupole_[ i ] = 0.0;
		my_group_[ i ].clear();
		if ( rsd.xyz(i).x() != rsd.xyz(i).x() ) {
			TR << "initialize stuff Problem in coord! " << rsd.seqpos() << " res " << rsd.name() << " atom " << rsd.atom_name( i ) << std::endl;
			std::exit(1);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MultipoleElecPoseInfo::MultipoleElecPoseInfo( MultipoleElecPoseInfo const & src ):
	CacheableData()
{
	Size const src_size( src.size() );

	residue_info_.resize( src_size );
	placeholder_residue_.resize( src_size );
	placeholder_info_.resize( src_size );

	for ( Size i=1; i<= src_size; ++i ) {
		residue_info_[i] = src.residue_info_[i]->copy_clone();
		if ( src.placeholder_residue_[i] ) {
			placeholder_residue_[i] = src.placeholder_residue_[i]->clone();
			placeholder_info_[i] = src.placeholder_info_[i]->copy_clone();
		} else {
			placeholder_residue_[i] = nullptr;
			placeholder_info_[i] = nullptr;
		}
	}
	being_packed_ = src.being_packed_;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
void
MultipoleElecPoseInfo::initialize( pose::Pose const & pose )
{
	Size const nres( pose.size() );

	residue_info_.resize( nres, nullptr );
	placeholder_residue_.resize( nres, nullptr );
	placeholder_info_.resize( nres, nullptr );

	for ( Size i=1; i<= nres; ++i ) {
		if ( !residue_info_[i] ) residue_info_[i] = MultipoleElecResidueInfoOP( new MultipoleElecResidueInfo( pose.residue(i) ) );
		else  residue_info_[i]->initialize( pose.residue(i) );
	}

	being_packed_.clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
MultipoleElecPoseInfo::set_placeholder( Size const i, ResidueOP rsd, MultipoleElecResidueInfoOP info )
{
	placeholder_residue_[ i ] = rsd;
	placeholder_info_[ i ] = info;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
MultipoleElecPoseInfo::set_repack_list( utility::vector1< bool > const & repacking_residues )
{
	being_packed_.resize( size(), false );
	for ( Size i=1; i<= size(); ++i ) {
		being_packed_[i] = repacking_residues[ i ];
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
MultipoleElecRotamerSetInfo::initialize( RotamerSetBase const & rotamer_set )
{
	Size const nrot( rotamer_set.num_rotamers() );
	residue_info_.resize( nrot );
	for ( Size i=1; i<= nrot; ++i ) {
		residue_info_[i] = MultipoleElecResidueInfoOP( new MultipoleElecResidueInfo( *rotamer_set.rotamer(i) ) );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultipoleElecPotential::read_in_amoeba_parameters() {

	/////////////////////////////////////
	// First read in multipole parameters
	/////////////////////////////////////
	std::string type_file_name = basic::database::full_name( "chemical/amoeba/amoeba09_types.txt" );
	std::ifstream types_file( type_file_name.c_str() );
	std::string input_line;
	while ( getline( types_file, input_line ) ) {
		//TR << "Processing line: " << input_line << std::endl;
		utility::vector1< std::string > tokens( utility::split_whitespace( input_line ) );
		Size const num_tokens( tokens.size() );
		Size const num_variant_qualifiers( num_tokens - 5 );
		//TR << "Found " << num_tokens << " tokens " << std::endl;
		std::string atomname( tokens[2] );
		utility::trim( atomname, " " );
		//TR << "Atom name " << atomname << " Residue name " << tokens[3];
		std::string new_key = atomname + tokens[3];
		for ( Size ivar = 1 ; ivar <= num_variant_qualifiers ; ++ivar ) {
			new_key = new_key + tokens[ ivar + 3 ];
			//TR << "   " << tokens[ ivar + 3 ];
		}
		//TR << " type index is:  " << tokens[ num_tokens ] << std::endl;
		type_lookup_[ new_key ] = static_cast<core::Size>( boost::lexical_cast< core::Size >( tokens[ num_tokens - 1 ] ) );
		//TR << "Stashing key X" << new_key << "X" << std::endl;
	}
	types_file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultipoleElecPotential::read_in_multipole_parameters() {

	std::string multipole_file_name = basic::database::full_name( "chemical/amoeba/amoeba09_multipole_params.txt" );
	std::ifstream multipole_file( multipole_file_name.c_str() );
	std::string input_line;

	MultipoleAxisType axis;
	Real chirality_check;
	Real mono_mom;
	Vector di_mom;
	Matrix quad_mom;

	utility::vector1< int > type_keys( 4, 0 );

	while ( getline( multipole_file, input_line ) ) {

		//TR << "Processing line: " << input_line << std::endl;

		utility::vector1< std::string > tokens( utility::split_whitespace( input_line ) );
		Size const num_tokens( tokens.size() );
		// Get the monopole moment
		mono_mom = static_cast<core::Real>( boost::lexical_cast< core::Real >( tokens[ num_tokens ] ) );

		// default value
		chirality_check = 1.0;

		//TR << "Found " << num_tokens << " tokens in line: " << "\n" << input_line << "\n";
		for ( Size k = 1 ; k <= type_keys.size(); ++k ) {
			type_keys[k] = 0.0;
		}

		// Read in the type keys
		for ( Size k = 1 ; k <= (num_tokens - 2) ; ++k ) {
			type_keys[ k ] = static_cast<core::Size>( boost::lexical_cast< core::Size >( tokens[ k + 1 ] ) );
		}

		//TR << "Keys " << type_keys[1] << "   " << type_keys[2] << "   " << type_keys[3] << "   " << type_keys[4] << std::endl;

		if ( type_keys[4] == 0 && type_keys[3] == 0 && type_keys[2] == 0 ) {
			//TR << "No multipole moments beyond monopole" << std::endl;
			axis = none;
			di_mom = 0.0;
			quad_mom = 0.0;
		} else if ( type_keys[4] == 0 && type_keys[3] == 0 ) {
			axis = z_axis_only;
			//TR << "Z-axis definition only" << std::endl;
		} else if ( type_keys[2] < 0 && type_keys[3] < 0 && type_keys[4] == 0 ) {
			axis = bisector;
			//TR << "Bisector axis definition" << std::endl;
			type_keys[2] = -type_keys[2];
			type_keys[3] = -type_keys[3];
		} else if ( type_keys[2] > 0 && type_keys[3] < 0 && type_keys[4] < 0 ) {
			axis = z_then_bisector;
			//TR << "Z-then-Bisector axis definition" << std::endl;
			type_keys[3] = -type_keys[3];
			type_keys[4] = -type_keys[4];
		} else if ( type_keys[2] < 0 && type_keys[3] < 0 && type_keys[4] < 0  ) {
			axis = three_fold;
			type_keys[2] = -type_keys[2];
			type_keys[3] = -type_keys[3];
			type_keys[4] = -type_keys[4];
			//TR << "Three-fold axis definition" << std::endl;
		} else {
			axis = z_then_x;
			if ( type_keys[4] < 0 ) {
				chirality_check = -1.0;
				type_keys[4] = -type_keys[4];
			}
			//TR << "Z-then-X axis definition" << std::endl;
		}

		multipole_file >> di_mom.x();
		multipole_file >> di_mom.y();
		multipole_file >> di_mom.z();

		multipole_file >> quad_mom.xx();
		multipole_file >> quad_mom.yx();
		multipole_file >> quad_mom.yy();
		multipole_file >> quad_mom.zx();
		multipole_file >> quad_mom.zy();
		multipole_file >> quad_mom.zz();

		// Set the symmetrically related elements
		quad_mom.xy() = quad_mom.yx();
		quad_mom.xz() = quad_mom.zx();
		quad_mom.yz() = quad_mom.zy();

		// Convert to Angstroms
		di_mom *= bohr;
		quad_mom *= (bohr*bohr/3.0);

		//TR << "Dipole components " << di_mom << std::endl;
		//TR << "Quadrupole components " << quad_mom << std::endl;

		// This is necessary to chomp to the end of the line
		getline( multipole_file, input_line );

		utility::vector1< Size > amoeba_type_keys( 4, 0 );
		for ( Size k = 1 ; k <= type_keys.size(); ++k ) {
			amoeba_type_keys[k] = static_cast< Size >( type_keys[k] );
		}

		// now stash in the multimap
		MultipoleParameter::MultipoleParameterOP new_mp_param( new MultipoleParameter( axis, amoeba_type_keys, chirality_check, mono_mom, di_mom, quad_mom ) );
		multipole_info_.insert( std::pair< Size, MultipoleParameter::MultipoleParameterOP >( amoeba_type_keys[1], new_mp_param ) );
	}

	multipole_file.close();

	///////////////////////////////////////
	// Next read in polarization parameters
	///////////////////////////////////////
	std::string polarization_file_name = basic::database::full_name( "chemical/amoeba/amoeba09_polarize_params.txt" );
	std::ifstream polarization_file( polarization_file_name.c_str() );

	while ( getline( polarization_file, input_line ) ) {
		utility::vector1< std::string > tokens( utility::split_whitespace( input_line ) );
		Size const num_tokens( tokens.size() );

		//TR << "Processing line: " << input_line << std::endl;
		//TR << "Found " << num_tokens << " tokens in line: " << "\n" << input_line << std::endl;

		Size const this_type( static_cast<core::Size>( boost::lexical_cast< core::Size >( tokens[2] ) ) );
		Real const polarity( static_cast<core::Real>( boost::lexical_cast< core::Real >( tokens[3] ) ) );
		Real const thole( static_cast<core::Real>( boost::lexical_cast< core::Real >( tokens[4] ) ) );
		utility::vector1< Size > group_members;

		if ( num_tokens > 4 ) {
			for ( Size k = 5 ; k <= num_tokens ; ++k ) {
				group_members.push_back(  static_cast<core::Size>( boost::lexical_cast< core::Size >( tokens[ k ] ) ) );
			}
			//TR << "This type has " << group_members.size() << " defined group types" << std::endl;
			//for( Size k = 1 ; k <= group_members.size() ; ++k ) {
			// TR << " Type " << group_members[k] << "   ";
			//}
			//TR << std::endl;
		}

		// Attach these values to all relevant type entries in the multimap
		// Get iterators for the range of possible multipole parameters
		std::pair<
			std::multimap< Size, MultipoleParameter::MultipoleParameterOP >::const_iterator,
			std::multimap< Size, MultipoleParameter::MultipoleParameterOP >::const_iterator >
			key_range = multipole_info_.equal_range( this_type );

		for ( auto it = key_range.first ;
				it != key_range.second; ++it ) {
			MultipoleParameter::MultipoleParameterOP this_mp( it->second );
			this_mp->polarity() = polarity;
			this_mp->thole() = thole;
			this_mp->pdamp() = std::pow( polarity, (1.0/6.0) );
			this_mp->my_group_members() = group_members;
		}
	}
	polarization_file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultipoleElecPotential::find_params_and_neighbors(
	core::pose::Pose const & pose,
	MultipoleParameterOP & mp_param,
	MultipoleElecResidueInfo & mp,
	core::conformation::Residue const & rsd,
	Size const j,
	Size const this_type
) const {
	using namespace id;
	//if ( rsd.is_virtual( j ) ) return;

	// I feel dirty doing this, but I really need to know if this
	// residue is part of the pose or just a rotamer. This checks whether
	// the raw address of this residue is the same as that of the residue
	// at the corresponding position in the pose.  If not, the residue is
	// taken to be a rotamer.  This is not an pretty way to hack this.
	bool const is_a_rotamer( &rsd != &(pose.residue( rsd.seqpos() )) );

	// if( is_a_rotamer ) {
	//  TR << "I must be a rotamer!" << std::endl;
	// } else {
	//  TR << "I am not a rotamer!" << std::endl;
	// }

	// MultipoleElecPoseInfoCOP multipole_info;
	MultipoleElecPoseInfo const & multipole_info( static_cast< MultipoleElecPoseInfo const & >( pose.data().get( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO)));

	// Get iterators for the range of possible multipole parameters
	std::pair<
		std::multimap< Size, MultipoleParameter::MultipoleParameterOP >::const_iterator,
		std::multimap< Size, MultipoleParameter::MultipoleParameterOP >::const_iterator >
		key_range = multipole_info_.equal_range( this_type );

	mp.my_local_coord_frame( j ).clear();

	// TR << "This type is " << this_type << std::endl;

	for ( auto it = key_range.first ;
			it != key_range.second; ++it ) {
		MultipoleParameter::MultipoleParameterOP this_mp( it->second );
		//  TR << "Type is " << this_mp->coord_type() << "  Keys are " << it->second->atom_type()[1] << "   " << it->second->atom_type()[2] << "   " << it->second->atom_type()[3] << "   " << it->second->atom_type()[4] << std::endl;
		// This part depends on the type of axis definition

		// Make a vector of AtomIDs to hold the bonded neighbors for this atom
		AtomID this_id( j, rsd.seqpos() );
		//  utility::vector1< AtomID > const & bonded( pose.conformation().bonded_neighbor_all_res( this_id ) );
		utility::vector1< AtomID > bonded;
		if ( is_a_rotamer ) {
			for ( Size idummy = 1 ; idummy <= rsd.nbrs( j ).size() ; ++idummy ) {
				bonded.push_back( AtomID( rsd.nbrs( j )[ idummy ], rsd.seqpos() ) );
			}
		} else {
			bonded = pose.conformation().bonded_neighbor_all_res( this_id );
		}

		// Try to be clever here for corner cases.  For rare rotamers that require an atom
		// in the neighboring residue in the pose (but that are not 'connected'), adding
		// an AtomID to the bonded vector may work.  The cost is that 'bonded' can no longer
		// be merely a reference to a vector, since it will be modified.  This involves copying.

		if ( is_a_rotamer &&
				rsd.has_upper_connect() &&
				( rsd.upper_connect_atom() == j ) &&
				pose.residue( rsd.seqpos() + 1 ).has_lower_connect() ) {
			//   TR << "Adding a connection neighbor" << std::endl;
			Size const nxt_atm( pose.residue( rsd.seqpos() + 1 ).lower_connect_atom() );
			Residue const & nxt_rsd( pose.residue( rsd.seqpos() + 1 ) );
			bonded.push_back( AtomID( nxt_atm, nxt_rsd.seqpos() ) );
		}

		//  TR << "Atom has " << bonded.size() << " bonded neighbors" << std::endl;
		//  for( Size idummy = 1 ; idummy <= bonded.size() ; ++idummy ) {
		//   MultipoleElecResidueInfo const & chk_mp( multipole_info.residue_info( bonded[idummy].rsd() ) );
		//   if( is_a_rotamer && ( rsd.seqpos() == bonded[ idummy ].rsd() ) ) {
		//    TR << "Doing the rotamer thing" << std::endl;
		//    TR << "Atom name " << rsd.atom_name( bonded[idummy].atomno() )  << std::endl;
		//    TR << "Atom index " << bonded[idummy].atomno() << std::endl;
		//    TR << "Rot style Rsd, atom " << rsd.seqpos() << ", " << j << " bonded neighbor " << idummy << " is rsd " << bonded[idummy].rsd() << " atom " << rsd.atom_name( bonded[idummy].atomno() ) << " type " << mp.type( bonded[idummy].atomno() ) << std::endl;
		//   } else {
		//    TR << "Doing the residue thing" << std::endl;
		//    TR << "Rsd, atom " << rsd.seqpos() << ", " << j << " bonded neighbor " << idummy << " is rsd " << bonded[idummy].rsd() << " atom " << pose.residue(bonded[idummy].rsd()).atom_name( bonded[idummy].atomno() ) << " type " << chk_mp.type( bonded[idummy].atomno() ) << std::endl;
		//   }
		//  }

		// First handle ions etc., that have no higher moments
		if ( this_mp->coord_type() == none ) {
			mp_param = this_mp;
			return;
			// Now handle atoms that have rotational symmetry
		} else if ( this_mp->coord_type() == z_axis_only ) {
			//TR << "Trying Z-axis only system" << std::endl;
			// Find the first atom with the required amoeba type in the residue
			for ( Size ichk = 1 ; ichk <= bonded.size() ; ichk++ ) {
				Size const chk_atom( bonded[ ichk ].atomno() );
				Size const chk_res( bonded[ ichk ].rsd() );
				bool use_local( is_a_rotamer && ( chk_res == rsd.seqpos() ) );
				MultipoleElecResidueInfo const & chk_mp( use_local ? mp : multipole_info.residue_info( chk_res ) );
				if ( this_mp->atom_type()[2] == chk_mp.type( chk_atom ) ) {
					mp_param = this_mp;
					mp.my_local_coord_frame( j ).push_back( bonded[ ichk ] );
					//TR << "Found match!" << std::endl;
					return;
				}
			}

			continue;
			// Now handle bisector systems
		} else if ( this_mp->coord_type() == bisector ) {
			//TR << "Trying Bisector system" << std::endl;
			// Find the first atom with the required amoeba type in the residue
			for ( Size ichk = 1 ; ichk <= bonded.size() ; ichk++ ) {
				Size const chk_atom( bonded[ ichk ].atomno() );
				Size const chk_res( bonded[ ichk ].rsd() );
				bool use_local( is_a_rotamer && ( chk_res == rsd.seqpos() ) );
				MultipoleElecResidueInfo const & chk_mp( use_local ? mp : multipole_info.residue_info( chk_res ) );
				if ( this_mp->atom_type()[2] != chk_mp.type( chk_atom ) ) continue;

				for ( Size ichk2 = ichk + 1 ; ichk2 <= bonded.size() ; ichk2++ ) {
					Size const chk_atom2( bonded[ ichk2 ].atomno() );
					Size const chk_res2( bonded[ ichk2 ].rsd() );
					bool use_local_2( is_a_rotamer && ( chk_res2 == rsd.seqpos() ) );
					MultipoleElecResidueInfo const & chk_mp2( use_local_2 ? mp : multipole_info.residue_info( chk_res2 ) );
					if ( this_mp->atom_type()[3] != chk_mp2.type( chk_atom2 ) ) continue;

					mp_param = this_mp;
					mp.my_local_coord_frame( j ).push_back( bonded[ ichk ] );
					mp.my_local_coord_frame( j ).push_back( bonded[ ichk2 ] );
					//TR << "Found match!" << std::endl;
					return;
				}
			}
			continue;
		} else if ( this_mp->coord_type() == z_then_bisector ) {
			//TR << "Trying Z-Then-Bisector system" << std::endl;
			// In the current set, all involved atoms are one bond away
			for ( Size ichk = 1 ; ichk <= bonded.size() ; ichk++ ) {
				Size const chk_atom( bonded[ ichk ].atomno() );
				Size const chk_res( bonded[ ichk ].rsd() );
				bool use_local( is_a_rotamer && ( chk_res == rsd.seqpos() ) );
				MultipoleElecResidueInfo const & chk_mp( use_local ? mp : multipole_info.residue_info( chk_res ) );
				if ( this_mp->atom_type()[2] != chk_mp.type( chk_atom ) ) continue;

				for ( Size ichk2 = 1 ; ichk2 <= bonded.size() ; ichk2++ ) {
					Size const chk_atom2( bonded[ ichk2 ].atomno() );
					Size const chk_res2( bonded[ ichk2 ].rsd() );
					bool use_local_2( is_a_rotamer && ( chk_res2 == rsd.seqpos() ) );
					MultipoleElecResidueInfo const & chk_mp2( use_local_2 ? mp : multipole_info.residue_info( chk_res2 ) );
					if ( this_mp->atom_type()[3] != chk_mp2.type( chk_atom2 ) ) continue;

					for ( Size ichk3 = ichk2 + 1 ; ichk3 <= bonded.size() ; ichk3++ ) {
						Size const chk_atom3( bonded[ ichk3 ].atomno() );
						Size const chk_res3( bonded[ ichk3 ].rsd() );
						bool use_local_3( is_a_rotamer && ( chk_res3 == rsd.seqpos() ) );
						MultipoleElecResidueInfo const & chk_mp3( use_local_3 ? mp : multipole_info.residue_info( chk_res3 ) );
						if ( this_mp->atom_type()[4] != chk_mp3.type( chk_atom3 ) ) continue;

						mp_param = this_mp;
						mp.my_local_coord_frame( j ).push_back( bonded[ ichk ] );
						mp.my_local_coord_frame( j ).push_back( bonded[ ichk2 ] );
						mp.my_local_coord_frame( j ).push_back( bonded[ ichk3 ] );
						//TR << "Found match!" << std::endl;
						return;
					}
				}
			}
		} else if ( this_mp->coord_type() == three_fold ) {
			//TR << "Trying Three-Fold system" << std::endl;
			// In the current set, all involved atoms are one bond away
			for ( Size ichk = 1 ; ichk <= bonded.size() ; ichk++ ) {
				Size const chk_atom( bonded[ ichk ].atomno() );
				Size const chk_res( bonded[ ichk ].rsd() );
				bool use_local( is_a_rotamer && ( chk_res == rsd.seqpos() ) );
				MultipoleElecResidueInfo const & chk_mp( use_local ? mp : multipole_info.residue_info( chk_res ) );
				if ( this_mp->atom_type()[2] != chk_mp.type( chk_atom ) ) continue;

				for ( Size ichk2 = 1 ; ichk2 <= bonded.size() ; ichk2++ ) {
					Size const chk_atom2( bonded[ ichk2 ].atomno() );
					Size const chk_res2( bonded[ ichk2 ].rsd() );
					bool use_local_2( is_a_rotamer && ( chk_res2 == rsd.seqpos() ) );
					MultipoleElecResidueInfo const & chk_mp2( use_local_2 ? mp : multipole_info.residue_info( chk_res2 ) );
					if ( this_mp->atom_type()[3] != chk_mp2.type( chk_atom2 ) ) continue;

					for ( Size ichk3 = 1 ; ichk3 <= bonded.size() ; ichk3++ ) {
						Size const chk_atom3( bonded[ ichk3 ].atomno() );
						Size const chk_res3( bonded[ ichk3 ].rsd() );
						bool use_local_3( is_a_rotamer && ( chk_res3 == rsd.seqpos() ) );
						MultipoleElecResidueInfo const & chk_mp3( use_local_3 ? mp : multipole_info.residue_info( chk_res3 ) );
						if ( this_mp->atom_type()[4] != chk_mp3.type( chk_atom3 ) ) continue;

						mp_param = this_mp;
						mp.my_local_coord_frame( j ).push_back( bonded[ ichk ] );
						mp.my_local_coord_frame( j ).push_back( bonded[ ichk2 ] );
						mp.my_local_coord_frame( j ).push_back( bonded[ ichk3 ] );
						//TR << "Found match!" << std::endl;
						return;
					}
				}
			}
		} else if ( this_mp->coord_type() == z_then_x ) {
			// This comes in two flavors.  The simple subtype only
			// specifies the z and x directions.  The other also
			// specifies a y atom, and is used at chiral sites.

			// First try checking for only directly bonded.  This can use
			// the bonded vector.  Below handles cases where some of the
			// axis-determining atoms are two bonds away, but in that case
			// all the neighboring atoms have to be in the same residue.
			// If you have a case where the atoms are more than one bond
			// away and in different residues, you are hosed.
			// As it stands, for proteins, all required atoms are in the
			// same residue.  For nucleic acids, you need an atom in
			// a different residue for  O3', but all required atoms are
			// one bond away.

			//TR << "Trying Z-Then-X system" << std::endl;

			for ( Size ichk = 1 ; ichk <= bonded.size() ; ichk++ ) {
				Size const chk_atom( bonded[ ichk ].atomno() );
				Size const chk_res( bonded[ ichk ].rsd() );
				bool use_local( is_a_rotamer && ( chk_res == rsd.seqpos() ) );
				MultipoleElecResidueInfo const & chk_mp( use_local ? mp : multipole_info.residue_info( chk_res ) );
				if ( this_mp->atom_type()[2] != chk_mp.type( chk_atom ) ) continue;

				for ( Size ichk2 = 1 ; ichk2 <= bonded.size() ; ichk2++ ) {
					Size const chk_atom2( bonded[ ichk2 ].atomno() );
					Size const chk_res2( bonded[ ichk2 ].rsd() );
					bool use_local_2( is_a_rotamer && ( chk_res2 == rsd.seqpos() ) );
					MultipoleElecResidueInfo const & chk_mp2( use_local_2 ? mp : multipole_info.residue_info( chk_res2 ) );
					if ( ( this_mp->atom_type()[3] != chk_mp2.type( chk_atom2 ) ) ||
							( ichk2 == ichk ) ) continue;

					if ( this_mp->atom_type()[4] == 0 ) {
						// This is the simple, achiral case
						mp_param = this_mp;
						mp.my_local_coord_frame( j ).push_back( bonded[ ichk ] );
						mp.my_local_coord_frame( j ).push_back( bonded[ ichk2 ] );
						//TR << "Found match!" << std::endl;
						return;
					} else {
						// Chiral case - look for the third match
						for ( Size ichk3 = 1 ; ichk3 <= bonded.size() ; ichk3++ ) {
							Size const chk_atom3( bonded[ ichk3 ].atomno() );
							Size const chk_res3( bonded[ ichk3 ].rsd() );
							bool use_local_3( is_a_rotamer && ( chk_res3 == rsd.seqpos() ) );
							MultipoleElecResidueInfo const & chk_mp3( use_local_3 ? mp : multipole_info.residue_info( chk_res3 ) );
							if ( this_mp->atom_type()[4] != chk_mp3.type( chk_atom3 ) ) continue;

							mp_param = this_mp;
							mp.my_local_coord_frame( j ).push_back( bonded[ ichk ] );
							mp.my_local_coord_frame( j ).push_back( bonded[ ichk2 ] );
							mp.my_local_coord_frame( j ).push_back( bonded[ ichk3 ] );
							//TR << "Found match!" << std::endl;
							return;
						}
					}
				}
			}

			for ( Size iatom = 1 ; iatom <= rsd.natoms() ; ++iatom ) {
				if ( ( rsd.path_distance( j, iatom ) != 1 ) ||
						( mp.type( iatom ) != this_mp->atom_type()[2] ) ) continue;

				for ( Size iatom2 = 1 ; iatom2 <= rsd.natoms() ; ++iatom2 ) {
					if ( ( ( rsd.path_distance( j, iatom2 )  != 1 ) &&
							( rsd.path_distance( j, iatom2 ) != 2 ) ) ||
							( mp.type( iatom2 ) != this_mp->atom_type()[3] ) ||
							( iatom == iatom2 ) ) continue;

					if ( this_mp->atom_type()[4] == 0 ) {
						// This is the simple, achiral case
						mp_param = this_mp;
						mp.my_local_coord_frame( j ).push_back( AtomID( iatom, rsd.seqpos()  ) );
						mp.my_local_coord_frame( j ).push_back( AtomID( iatom2, rsd.seqpos() ) );
						//TR << "Found match!" << std::endl;
						return;
					} else {
						// Chiral case - look for the third match
						for ( Size iatom3 = 1 ; iatom3 <= rsd.natoms() ; ++iatom3 ) {
							if ( ( ( rsd.path_distance( j, iatom3 )  != 1 ) &&
									( rsd.path_distance( j, iatom3 ) != 2 ) ) ||
									( mp.type( iatom3 ) != this_mp->atom_type()[4] ) ) continue;

							mp_param = this_mp;
							mp.my_local_coord_frame( j ).push_back( AtomID( iatom, rsd.seqpos() ) );
							mp.my_local_coord_frame( j ).push_back( AtomID( iatom2, rsd.seqpos() ) );
							mp.my_local_coord_frame( j ).push_back( AtomID( iatom3, rsd.seqpos() ) );
							//TR << "Found match!" << std::endl;
							return;
						}
					}
				}
			}

		} else {
			TR << "Error - multipole type not found!" << std::endl;
		}

	}

	TR << "Error - All matching multipole parameters exhausted without finding a match!" << std::endl;
	TR << "Residue " << rsd.name() << " atom " << rsd.atom_name( j ) << " position " << rsd.seqpos() << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

core::Size
MultipoleElecPotential::amoeba_type_lookup(
	std::string const & atomname,
	std::string const & resname,
	std::string const & variantname,
	std::string const & variantname2
) const {

	std::string const not_variant( "NONE" );
	core::Size type( 0 );

	// Lookup default value w/o variant info
	std::string tmp_key = atomname + resname;
	tmp_key.erase( std::remove_if( tmp_key.begin(), tmp_key.end(), ::isspace ), tmp_key.end() );
	std::map< std::string, core::Size >::const_iterator map_it;

	//TR << "Querying key X" << tmp_key << "X" << std::endl;
	map_it = type_lookup_.find( tmp_key );
	if ( map_it != type_lookup_.end() ) {
		type = map_it->second;
	}

	// Check for variants
	if ( !utility::trimmed_compare( variantname, not_variant ) ) {
		std::string tmp_key2 = atomname + resname + variantname;
		tmp_key2.erase( std::remove_if( tmp_key2.begin(), tmp_key2.end(), ::isspace ), tmp_key2.end() );
		// Only overwrite default if an entry is found

		//TR << "Querying key X" << tmp_key2 << "X" << std::endl;
		map_it = type_lookup_.find( tmp_key2 );
		if ( map_it != type_lookup_.end() ) {
			type = map_it->second;
		}

		// Try for double variants - test both orderings of variant strings
		if ( !utility::trimmed_compare( variantname2, not_variant ) ) {
			std::string single_key2 = atomname + resname + variantname2;
			single_key2.erase( std::remove_if( single_key2.begin(), single_key2.end(), ::isspace ), single_key2.end() );
			map_it = type_lookup_.find( single_key2 );
			if ( map_it != type_lookup_.end() ) {
				type = map_it->second;
			}

			std::string double_key1 = atomname + resname + variantname + variantname2;
			std::string double_key2 = atomname + resname + variantname2 + variantname;
			double_key1.erase( std::remove_if( double_key1.begin(), double_key1.end(), ::isspace ), double_key1.end() );
			double_key2.erase( std::remove_if( double_key2.begin(), double_key2.end(), ::isspace ), double_key2.end() );
			//TR << "Querying key X" << double_key1 << "X" << std::endl;
			map_it = type_lookup_.find( double_key1 );
			if ( map_it != type_lookup_.end() ) {
				type = map_it->second;
			}
			//TR << "Querying key X" << double_key2 << "X" << std::endl;
			map_it = type_lookup_.find( double_key2 );
			if ( map_it != type_lookup_.end() ) {
				type = map_it->second;
			}
		}
	}

	// Give a holler if nothing has been found
	if ( type == 0 ) {
		TR << "PROBLEM - NO BIOTYPE FOUND FOR " << atomname << " " << resname << " " << variantname << std::endl;
	}

	return type;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultipoleElecPotential::align_residue_multipole_axes(
	core::pose::Pose const & pose,
	Residue const & rsd,
	MultipoleElecResidueInfo & mp
) const {
	// First find a reference atom for each in the residue.  This atom
	// must be a direct bonded neighbor, but must not be in the torsion
	// path.  The easiest guess is an attached hydrogen.  After that,
	// find an attached heavy atom that isn't in the torsion path.  After
	// this, find out what's left.  I imagine the sulfur in the methionine
	// side chain may be problematic, but this will only be a problem when
	// a jump is defined involving the terminal carbon on the Met sidechain,
	// and when folding proceeds inward from the jump.

	// Find all disallowed atoms in backbone
	utility::vector1< Size > torsion_atoms;
	chemical::AtomIndices mainchain( rsd.type().mainchain_atoms() );
	for ( Size iat = 1 ; iat <= mainchain.size() ; ++iat ) {
		//TR << "Disallowing MC " << rsd.atom_name( mainchain[iat] ) << std::endl;
		torsion_atoms.push_back( mainchain[iat] );
	}

	// Find all disallowed atoms in sidechain
	for ( Size ichi = 1 ; ichi <= rsd.nchi() ; ++ichi ) {
		// The atoms that could be a problem are the central two.
		//TR << "Disallowing SC " << rsd.atom_name( rsd.chi_atoms( ichi )[2] ) << std::endl;
		torsion_atoms.push_back( rsd.chi_atoms( ichi )[2] );
		//TR << "Disallowing SC " << rsd.atom_name( rsd.chi_atoms( ichi )[3] ) << std::endl;
		torsion_atoms.push_back( rsd.chi_atoms( ichi )[3] );
	}

	for ( Size iat = 1 ; iat <= rsd.natoms() ; ++iat ) {
		//   if( rsd.is_virtual( iat ) ) continue;
		if ( rsd.xyz(iat).x() != rsd.xyz(iat).x() ) {
			TR << "align_mult_axes top stuff Problem in coord! " << rsd.seqpos() << " res " << rsd.name() << " atom " << rsd.atom_name( iat ) << std::endl;
			std::exit(1);
		}
		// Test for hydrogens and atoms that are not in the torsion path
		if ( rsd.atom_is_hydrogen( iat ) ||
				( std::find( torsion_atoms.begin(), torsion_atoms.end(), iat ) == torsion_atoms.end() )
				) {
			//mp.set_coord_frame_ref( iat, rsd.atom_base( iat ) );
			// A zero implies that no coordinate frame correction is required
			// in the derivative calculation.
			id::AtomID dummy_atom( 0, rsd.seqpos() );
			mp.set_coord_frame_ref( iat, dummy_atom );
			//TR << "Coord frame ref for " << rsd.atom_name( iat ) << " is base atom " << rsd.atom_name( rsd.atom_base( iat ) ) << " skip frame derivs!" << std::endl;
		} else {
			// Atoms in the torsion path.
			// Loop backwards through bonded neighbors, the better to find hydrogens first
			chemical::AtomIndices const my_nbrs( rsd.bonded_neighbor( iat ) );
			for ( Size inbr = my_nbrs.size() ; inbr >= 1 ; --inbr ) {
				Size const nbr_index( my_nbrs[ inbr ] );
				if ( rsd.is_virtual( nbr_index ) ) continue;
				if ( rsd.atom_is_hydrogen( nbr_index ) ) {
					id::AtomID nbr_atom( nbr_index, rsd.seqpos() );
					mp.set_coord_frame_ref( iat, nbr_atom );
					//TR << "Coord frame ref for " << rsd.atom_name( iat ) << " is " << rsd.atom_name( nbr_index ) << std::endl;
					break;
				} else if ( std::find( torsion_atoms.begin(), torsion_atoms.end(), nbr_index ) == torsion_atoms.end() ) {
					id::AtomID nbr_atom( nbr_index, rsd.seqpos() );
					mp.set_coord_frame_ref( iat, nbr_atom );
					//TR << "Coord frame ref for " << rsd.atom_name( iat ) << " is " << rsd.atom_name( nbr_index ) << std::endl;
					break;
				}
			}
			// Bad news.  We haven't found a good reference atom.  We will need
			// to use an atom in the torsion path.  We will have to peek at the
			// fold tree to see which direction along the path we should use -
			// should be the downstream atom with respect to folding direction.
			// This is especially a problem when the downstream atom is in another
			// residue.  As is always the case for non-terminal DNA residues.

			// Ok, can I find my options - that is, the upstream and downstream atoms
			// in the torsion path?

			if ( mp.coord_frame_ref( iat ).atomno() == 0 ) {

				//TR << std::endl;
				//TR << "Looking for ref atom for residue " << rsd.seqpos() << " index " << iat << std::endl;
				//TR << "Details " << rsd.name3() << " position " << rsd.seqpos() << " and " << rsd.atom_name( iat ) << std::endl;

				Size this_mc_index( 0 );
				for ( Size ichk = 1 ; ichk <= rsd.mainchain_atoms().size() ; ++ichk ) {
					if ( iat == rsd.mainchain_atom( ichk ) ) {
						this_mc_index = ichk;
						//TR << "This is mainchain atom " << ichk << " out of " << rsd.mainchain_atoms().size() << " total " << std::endl;
					}
				}

				if ( this_mc_index > 0 ) {

					// If we are the start, we need to go one residue earlier
					Size lower_index;
					Size lower_res;
					if ( this_mc_index == 1 ) {
						Residue const & lower_rsd( pose.residue( rsd.seqpos() - 1 ) );
						lower_index = lower_rsd.mainchain_atom( lower_rsd.mainchain_atoms().size() );
						lower_res = rsd.seqpos() - 1;
						//TR << "lower atom is " << lower_rsd.atom_name( lower_index ) << " in previous residue" << std::endl;
					} else {
						lower_index = rsd.mainchain_atom( this_mc_index - 1 );
						lower_res = rsd.seqpos();
						//TR << "lower atom is " << rsd.atom_name( lower_index ) << std::endl;
					}

					// If we are the end, we need to go one residue farther
					Size upper_index;
					Size upper_res;
					if ( this_mc_index == rsd.mainchain_atoms().size() ) {
						Residue const & upper_rsd( pose.residue( rsd.seqpos() + 1 ) );
						upper_index = upper_rsd.mainchain_atom( 1 );
						upper_res = rsd.seqpos() + 1;
						//TR << "upper atom is " << upper_rsd.atom_name( upper_index ) << " in next residue" << std::endl;
					} else {
						upper_index = rsd.mainchain_atom( this_mc_index + 1 );
						upper_res = rsd.seqpos();
						//TR << "upper atom is " << rsd.atom_name( upper_index ) << std::endl;
					}

					// Now that we know the two choices, need to figure out which we need,
					// which is equivalent to knowing the direction of the fold tree in
					// this segment.

					id::AtomID this_atom_id( iat, rsd.seqpos() );
					kinematics::tree::Atom const & this_atom( pose.atom_tree().atom( this_atom_id ) );
					//this_atom.show( 1 );
					kinematics::tree::AtomCOP parent_atom( this_atom.parent() );
					id::AtomID const & parent_atom_id( parent_atom->id() );
					//TR << "Parent atom id info " << parent_atom_id << std::endl;
					//Residue const & parent_rsd( pose.residue( parent_atom_id.rsd() ) );
					Size const parent_res( parent_atom_id.rsd() );
					Size const parent_atm( parent_atom_id.atomno() );
					//     TR << "This atom's parent is " << parent_rsd.name3() << " position " << parent_rsd.seqpos() << " name " << parent_rsd.atom_name( parent_atom_id.atomno() ) << std::endl;

					if ( parent_res == upper_res && parent_atm == upper_index ) {
						id::AtomID nbr_atom( lower_index, lower_res );
						mp.set_coord_frame_ref( iat, nbr_atom );
						//      TR << "Coord frame ref for " << rsd.atom_name( iat ) << " is " << nbr_atom << std::endl;
					} else if (  parent_res == lower_res && parent_atm == lower_index ) {
						id::AtomID nbr_atom( upper_index, upper_res );
						mp.set_coord_frame_ref( iat, nbr_atom );
						//      TR << "Coord frame ref for " << rsd.atom_name( iat ) << " is " << nbr_atom << std::endl;
					} else {
						TR << "Parent doesn't match either option!" << std::endl;
						std::exit( 1 );
					}
				}
			}
		}
	}

	for ( Size iat = 1 ; iat <= rsd.natoms() ; ++iat ) {
		// assert( mp.coord_frame_ref( iat ) != 0 );
		if ( mp.coord_frame_ref( iat ).atomno() == 0 ) {
			//   TR << "Whoops - no suitable neighbors for " << rsd.name() << " atom " << rsd.atom_name( iat ) << std::endl;
		}
	}




	// Rotate all dipole and quadrupole moments into the global frame
	for ( Size j = 1 ; j <= rsd.natoms() ; j++ ) {
		if ( rsd.xyz(j).x() != rsd.xyz(j).x() ) {
			TR << "align_mult_axes pre stuff Problem in coord! " << rsd.seqpos() << " res " << rsd.name() << " atom " << rsd.atom_name( j ) << std::endl;
			std::exit(1);
		}
		Size this_type( mp.type( j ) );
		//TR << "Residue " << rsd.name3() << " Atom " << rsd.atom_name( j ) << " amoeba type is: " << this_type << std::endl;
		//  TR << "Residue " << rsd.name() << " Atom " << rsd.atom_name( j ) << " amoeba type is: " << this_type << std::endl;
		// This function finds multipole parameter set and atom numbers for
		// the z, x, and y determining neighbor atoms - note that all neighbors
		// are assumed to be in the same residue.

		// This information could/should be cached, but must know when to be updated
		MultipoleParameter::MultipoleParameterOP mp_param;
		find_params_and_neighbors( pose, mp_param, mp, rsd, j, this_type );

		// Build a matrix for the coordinate system - for the three-fold systems,
		// this also involves a chirality check.  Then use this to rotate the dipole
		// and quadrupole moments into the global frame and store.
		build_frame_and_rotate( pose, mp_param, j, mp, rsd );

		if ( rsd.xyz(j).x() != rsd.xyz(j).x() ) {
			TR << "align_mult_axes post stuff Problem in coord! " << rsd.seqpos() << " res " << rsd.name() << " atom " << rsd.atom_name( j ) << std::endl;
			std::exit(1);
		}

		// Store the mp_param we found
		mp.nonconst_mp_param(j) = mp_param;
	}
}

void
MultipoleElecPotential::align_multipole_axes(
	pose::Pose & pose
) const {
	PROF_START( basic::MULTIPOLE_SETUP );

	Size const nres( pose.size() );

	//  * Identify the neighboring atoms that define the coord set
	//  * Build the orthonormal axis
	//  * Rotate the dipole and quadrupole info into the global frame
	//
	MultipoleElecPoseInfoOP multipole_info;
	multipole_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );

	// Rotate all dipole and quadrupole moments into the global frame
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		align_residue_multipole_axes( pose, rsd, mp1 );
	}

	// Check for bad dipole
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd.natoms() ; ++atm1 ) {

			Vector const & dipole( mp1.dipole( atm1 ) );
			if ( dipole != dipole ) {
				TR << "align_multi_axis Problem in dipole! " << rsd.seqpos() << " res " << rsd.name() << " atom " << rsd.atom_name( atm1 ) << std::endl;
				std::exit(1);
			}
		}
	}

	pose.data().set( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO, multipole_info );

	PROF_STOP( basic::MULTIPOLE_SETUP );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultipoleElecPotential::build_frame_and_rotate(
	core::pose::Pose const & pose,
	MultipoleParameterOP & mp_param,
	Size orig_atom,
	MultipoleElecResidueInfo & mp,
	core::conformation::Residue const & rsd
) const {
	//if ( rsd.is_virtual( orig_atom ) ) return;

	bool const is_a_rotamer( &rsd != &(pose.residue( rsd.seqpos() )) );

	utility::vector1< id::AtomID > & atms( mp.my_local_coord_frame( orig_atom ) );

	Size const atom1( atms.size() > 0 ? atms[1].atomno() : 0 );
	Size const atom2( atms.size() > 1 ? atms[2].atomno() : 0 );
	Size const atom3( atms.size() > 2 ? atms[3].atomno() : 0 );
	Size const res1( atms.size() > 0 ? atms[1].rsd() : 0 );
	Size const res2( atms.size() > 1 ? atms[2].rsd() : 0 );
	Size const res3( atms.size() > 2 ? atms[3].rsd() : 0 );

	// TR << "atms size is " << atms.size() << std::endl;
	// TR << "atom numbers are  " << atom1 << "  " << atom2 << "  " << atom3 << std::endl;
	// TR << "res  numbers are  " << res1 << "  " << res2 << "  " << res3 << std::endl;

	MultipoleAxisType const axis_type( mp_param->coord_type() );
	if ( axis_type == none ) {
		mp.nonconst_monopole( orig_atom ) = mp_param->monopole();
		mp.nonconst_dipole( orig_atom ) = mp_param->dipole();
		mp.nonconst_quadrupole( orig_atom ) = mp_param->quadrupole();
		if ( mp.dipole( orig_atom ) != mp.dipole( orig_atom ) ) {
			TR << "Bad dipole for the none case" << std::endl;
		}
		return;
	}

	Residue const & r1( is_a_rotamer && ( res1 == rsd.seqpos() ) ? rsd : pose.residue(res1) );
	Residue const & r2( is_a_rotamer && ( res2 == rsd.seqpos() ) ? rsd : pose.residue(res2) );
	Residue const & r3( ( is_a_rotamer && ( res3 == rsd.seqpos() ) ) || ( res3 == 0 ) ? rsd : pose.residue(res3) );

	// Matrix construction has some system-dependent behavior
	Vector const v1( ( r1.xyz( atom1 ) - rsd.xyz( orig_atom ) ).normalize() );
	Vector const v2( ( r2.xyz( atom2 ) - rsd.xyz( orig_atom ) ).normalize() );
	Vector const v3( ( atom3 == 0 ? Vector( 0.0 ) : ( r3.xyz( atom3 ) - rsd.xyz( orig_atom ) ).normalize() ) );

	//TR << "Axis type is " << axis_type << std::endl;

	Vector vz, vx, vy;

	// With the exception of the strictly monopole case, this
	// is where the coordinate axes are frobbed from the input
	// atom vectors.
	switch( axis_type ) {
	case z_axis_only :
		vz = v1;
		vx = Vector( numeric::random::uniform(), numeric::random::uniform(), numeric::random::uniform() );
		vx.project_normal( vz ).normalized();
		break;
	case bisector :
		vz = ( v1 + v2 ).normalized();
		vx = v2.projected_normal( vz ).normalized();
		break;
	case z_then_bisector :
		vz = v1;
		vx = ( v2 + v3 ).normalized();
		vx.project_normal( vz ).normalized();
		break;
	case z_then_x :
		vz = v1;
		vx = v2.projected_normal( vz ).normalized();
		break;
	case three_fold :
		vz = ( v1 + v2 + v3 ).normalized();
		vx = v2.projected_normal( vz ).normalized();
		break;
	default :
		TR << "Error - unhandled axis type" << std::endl;
		return;
	}

	// Now get the y-axis
	vy = vz.cross( vx );

	// Build the rotation matrix from the axes

	mp.local_coord_matrix( orig_atom ) = ( Matrix::cols_constructor( vx, vy, vz ) );

	Matrix const & rotate( mp.local_coord_matrix( orig_atom ) );

	// Rotate and store
	mp.nonconst_monopole( orig_atom ) = mp_param->monopole();
	mp.nonconst_dipole( orig_atom ) = rotate * mp_param->dipole();
	mp.nonconst_induced_dipole( orig_atom ) = rotate * mp.stored_induced_dipole( orig_atom );
	mp.nonconst_induced_rf_dipole( orig_atom ) = rotate * mp.stored_induced_rf_dipole( orig_atom );
	mp.nonconst_quadrupole( orig_atom ) = rotate * mp_param->quadrupole() * rotate.transposed();

	if ( mp.dipole( orig_atom ) != mp.dipole( orig_atom ) ) {
		TR << "Residue: " << rsd.name3() << " Atom: " << rsd.atom_name( orig_atom ) << std::endl;
		TR << "Bad dipole from type dipole " << mp_param->dipole() << " and rotate matrix " << rotate << std::endl;
		TR << "Other atoms res,atm: res " << res1 << " atm " << atom1 << " res2 " << res2 << " atm2 " << atom2 << " res3 " << res3 << " atm3 " << atom3 << std::endl;
		std::exit(1);
	}

	// Check chirality

	if ( axis_type == z_then_x && mp_param->atom_type()[4] != 0 ) {
		//TR << "Checking chirality!" << std::endl;
		// Check chirality
		Vector const chk1( rsd.xyz( orig_atom ) - pose.residue(res3).xyz( atom3 ) );
		Vector const chk2( pose.residue(res1).xyz( atom1 ) - pose.residue(res3).xyz( atom3 ) );
		Vector const chk3( pose.residue(res2).xyz( atom2 ) - pose.residue(res3).xyz( atom3 ) );

		Vector const chk4( ( chk2.y()*chk3.z() - chk2.z()*chk3.y() ),
			( chk3.y()*chk1.z() - chk3.z()*chk1.y() ),
			( chk1.y()*chk2.z() - chk1.z()*chk2.y() ) );

		core::Real volume_check( chk1.x()*chk4.x() + chk2.x()*chk4.y() + chk3.x()*chk4.z() );

		//TR << "Chirality check volume " << volume_check << " sign " << mp_param->chirality_sign() << std::endl;
		//TR << "Check1 " << chk1 << std::endl;
		//TR << "Check2 " << chk2 << std::endl;
		//TR << "Check3 " << chk3 << std::endl;
		//TR << "Cross " << chk4 << std::endl;

		if ( ( volume_check < 0.0 && mp_param->chirality_sign() > 0.0 ) ||
				( volume_check > 0.0 && mp_param->chirality_sign() < 0.0 ) ) {
			// Invert all multipole components involving the y-axis
			//TR << "Inverting chirality!" << std::endl;
			mp.nonconst_dipole( orig_atom ).y() *= -1.0;
			mp.nonconst_induced_dipole( orig_atom ).y() *= -1.0;
			mp.nonconst_induced_rf_dipole( orig_atom ).y() *= -1.0;
			mp.nonconst_quadrupole( orig_atom ).xy() *= -1.0;
			mp.nonconst_quadrupole( orig_atom ).yy() *= -1.0;
			mp.nonconst_quadrupole( orig_atom ).zy() *= -1.0;
			mp.nonconst_quadrupole( orig_atom ).yx() *= -1.0;
			mp.nonconst_quadrupole( orig_atom ).zx() *= -1.0;
		}
	}

	// Diagnostic output for comparison to tinker

	//TR << "Residue: " << rsd.name3() << " Atom: " << rsd.atom_name( orig_atom ) << " monopole " << mp.monopole( orig_atom ) << std::endl;
	// TR << "Dipole: " << mp.dipole( orig_atom ) << std::endl;
	//Matrix const & m( mp.quadrupole( orig_atom ));
	//TR << "Quadrupole: " << m.xx() << ' ' <<  m.xy() << ' ' <<  m.xz() << ' '
	// <<  m.yx() << ' ' <<  m.yy() << ' ' <<  m.yz() << ' '
	// <<  m.zx() << ' ' <<  m.zy() << ' ' <<  m.zz() << std::endl;

	return;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultipoleElecPotential::assign_residue_amoeba_type(
	Residue const & rsd,
	MultipoleElecResidueInfo & mp
) const {
	utility::vector1< std::string > parsed_resname( utility::string_split_simple( rsd.name(), ':' ) );
	std::string resname( parsed_resname[1] );

	// Handle variants
	std::string variantname( "NONE" );
	std::string variantname2( "NONE" );
	if ( parsed_resname.size() == 2 ) {
		//TR << "Using variant specialization " << parsed_resname[2] << std::endl;
		variantname = parsed_resname[2];
	} else if ( parsed_resname.size() > 2 ) {
		//   TR << "DOUBLE VARIANT SITUATION!!!" << std::endl;
		//   TR << "Full resname is " << rsd.name() <<  std::endl;
		//   TR << "Using parsed resname " << parsed_resname[2] <<  std::endl;
		//   TR << "Using first variant name only!!!" << std::endl;
		variantname = parsed_resname[2];
		variantname2 = parsed_resname[3];
		//for( Size ivar = 1 ; ivar <= parsed_resname.size() ; ++ivar ) {
		//TR << "Residue info " << parsed_resname[ ivar ] << std::endl;
		//}
	} else if ( parsed_resname.size() > 3 ) {
		TR << "TRIPLE  VARIANT SITUATION -> NOT YET CODED!!!" << std::endl;
		TR << "Full resname is " << rsd.name() <<  std::endl;
		TR << "Using parsed resname " << parsed_resname[2] <<  std::endl;
		TR << "Using first variant name only!!!" << std::endl;
		variantname = parsed_resname[2];
		for ( Size ivar = 1 ; ivar <= parsed_resname.size() ; ++ivar ) {
			TR << "Residue info " << parsed_resname[ ivar ] << std::endl;
		}
	}

	for ( Size j = 1 ; j <= rsd.natoms() ; j++ ) {
		if ( rsd.xyz(j).x() != rsd.xyz(j).x() ) {
			//TR << "assign_type stuff Problem in coord! " << rsd.seqpos() << " res " << rsd.name() << " atom " << rsd.atom_name( j ) << std::endl;
			std::exit(1);
		}
		//   TR << "Residue " << resname << " Atom " << rsd.atom_name( j )   << std::endl;
		// Lookup amoeba type
		core::Size this_type( amoeba_type_lookup( rsd.atom_name(j), resname, variantname, variantname2 ) );
		//   TR << "Found Amoeba type " << this_type << std::endl;
		mp.set_type( j, this_type );
	}
}

void
MultipoleElecPotential::assign_all_amoeba_types(
	pose::Pose & pose
) const {
	PROF_START( basic::MULTIPOLE_SETUP );

	Size const nres( pose.size() );

	// Look up the Amoeba type information
	MultipoleElecPoseInfoOP multipole_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) ) {
		multipole_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );
	} else {
		multipole_info = MultipoleElecPoseInfoOP( new MultipoleElecPoseInfo() );
	}

	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		assign_residue_amoeba_type( rsd, mp1 );
	}

	pose.data().set( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO, multipole_info );
	PROF_STOP( basic::MULTIPOLE_SETUP );
}



void
MultipoleElecPotential::determine_polarization_groups(
	pose::Pose & pose
) const {
	PROF_START( basic::MULTIPOLE_SETUP );

	Size const nres( pose.size() );

	// Look up the Amoeba type information
	MultipoleElecPoseInfoOP multipole_info;
	multipole_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );

	// Temporary storage of directly bonded group neighbors
	std::map< id::AtomID, utility::vector1< id::AtomID > > direct_group_neighbors;

	// Determine all of the directly bonded atoms that are in the same group,
	// then spread outward from there.
	for ( Size res1 = 1 ; res1 <= nres; ++res1 ) {
		Residue const & rsd1( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd1.natoms() ; ++atm1 ) {
			if ( rsd1.xyz(atm1).x() != rsd1.xyz(atm1).x() ) {
				TR << "det_pol_groups Problem in coord! " << rsd1.seqpos() << " res " << rsd1.name() << " atom " << rsd1.atom_name( atm1 ) << std::endl;
				std::exit(1);
			}
			Vector const old_value( mp1.induced_dipole( atm1 ) );
			utility::vector1< id::AtomID > group;

			if ( rsd1.is_virtual( atm1 ) ) continue;

			// Seed with directly bonded atoms
			id::AtomID this_atom( atm1, res1 );
			utility::vector1<core::id::AtomID> const & bonded_neighbors(
				pose.conformation().bonded_neighbor_all_res( this_atom ) );
			for ( Size inbr = 1 ; inbr <= bonded_neighbors.size() ; inbr++ ) {
				id::AtomID nbr_atom( bonded_neighbors[inbr] );
				Size atm2( nbr_atom.atomno() );
				Size res2( nbr_atom.rsd() );
				//Residue const & nbr_res( pose.residue( res2 ) );
				//if( nbr_res.is_virtual( atm2 ) ) continue;
				MultipoleElecResidueInfo & nbr_mp( multipole_info->residue_info( res2 ) );
				core::Size const nbr_type( nbr_mp.type( atm2 ) );
				if ( std::find( mp1.mp_param(atm1)->my_group_members().begin(), mp1.mp_param(atm1)->my_group_members().end(), nbr_type ) !=
						mp1.mp_param(atm1)->my_group_members().end() ) {
					group.push_back( nbr_atom );
				}
			}
			direct_group_neighbors[ this_atom ] = group;
		}
	}

	// Now expand
	for ( Size res1 = 1 ; res1 <= nres; ++res1 ) {
		Residue const & rsd1( pose.residue( res1 ) );
		//MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd1.natoms() ; ++atm1 ) {
			if ( rsd1.is_virtual( atm1 ) ) continue;

			id::AtomID this_atom( atm1, res1 );
			utility::vector1< id::AtomID > new_group( direct_group_neighbors[ this_atom ] );

			// Walk atoms until exhausted
			while ( new_group.size() != 0 ) {
				utility::vector1< id::AtomID > work_group( new_group );
				new_group.clear();
				for ( Size cur_nbr = 1 ; cur_nbr <= work_group.size() ; ++cur_nbr ) {
					id::AtomID cur_nbr_atm( work_group[ cur_nbr ] );
					//TR << "Looking at bonded neighbors of " << pose.residue(cur_nbr_atm.rsd()).name() << " atom " << pose.residue(cur_nbr_atm.rsd()).atom_name( cur_nbr_atm.atomno() ) << std::endl;
					// Now check out all of this atom's directly bonded neighbors
					utility::vector1<core::id::AtomID> const & nbr_neighbors( direct_group_neighbors[ cur_nbr_atm ]  );
					for ( Size inbr = 1 ; inbr <= nbr_neighbors.size() ; inbr++ ) {
						id::AtomID nbr_atom( nbr_neighbors[inbr] );
						// To be added now, the type must be correct, and the AtomID can't
						// already be in the current vector.
						//TR << "Trying " << nbr_res.name() << " atom " << nbr_res.atom_name( atm2 ) << std::endl;
						if ( std::find( direct_group_neighbors[ this_atom ].begin(),
								direct_group_neighbors[ this_atom ].end(), nbr_atom ) ==
								direct_group_neighbors[ this_atom ].end() ) {
							//TR << "Found an indirect group member!" << std::endl;
							direct_group_neighbors[ this_atom ].push_back( nbr_atom );
							new_group.push_back( nbr_atom );
						}
					}
				}
			}
		}
	}


	for ( Size res1 = 1 ; res1 <= nres; ++res1 ) {
		Residue const & rsd1( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd1.natoms() ; ++atm1 ) {
			id::AtomID this_atom( atm1, res1 );
			//utility::vector1< id::AtomID > & group( direct_group_neighbors[ this_atom ] );

			mp1.my_group( atm1 ) = direct_group_neighbors[ this_atom ];

			//TR << "Polarization group members for " << rsd1.name() << " atom " << rsd1.atom_name( atm1 ) << " are: " << std::endl;
			// Silly, temporary output
			//  for( Size inbr = 1 ; inbr <= group.size() ; ++inbr ) {
			//   Size res2( group[ inbr ].rsd() );
			//   Size atm2( group[ inbr ].atomno() );
			//   Residue const & rsd2( pose.residue( res2 ) );
			//   TR << " " << rsd2.name() << " atom " << rsd2.atom_name( atm2 ) << std::endl;
			//  }
		}
	}
	PROF_STOP( basic::MULTIPOLE_SETUP );
}

void
MultipoleElecPotential::calculate_fixed_fields_for_polarization(
	pose::Pose & pose
) const {
	PROF_START( basic::MULTIPOLE_SETUP );
	Size const nres( pose.size() );

	MultipoleElecPoseInfoOP multipole_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) ) {
		multipole_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );
	} else {
		multipole_info = MultipoleElecPoseInfoOP( new MultipoleElecPoseInfo() );
	}

	// Zero out the fixed field members
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd.natoms() ; ++atm1 ) {
			mp1.Efield_fixed( atm1 ) = 0.0;
			mp1.Efield_rf_fixed( atm1 ) = 0.0;
		}
	}

	// Double loop to accumulate E field contributions
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd1( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size res2 = res1; res2 <= nres; ++res2 ) {
			Residue const & rsd2( pose.residue( res2 ) );
			MultipoleElecResidueInfo & mp2( multipole_info->residue_info( res2 ) );
			calculate_res_res_fixed_fields_for_polarization( rsd1, mp1, rsd2, mp2 );
		}
	}

	// Print out energy field values
	// for ( Size res1 = 1; res1 <= nres; ++res1 ) {
	//  Residue const & rsd1( pose.residue( res1 ) );
	//  MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
	//  for ( Size atm1 = 1 ; atm1 <= rsd1.natoms() ; ++atm1 ) {
	//   TR << "Fixed Electric field for res " << res1 << " atom " << rsd1.atom_name( atm1 ) << " at position " << rsd1.xyz( atm1 ) << " is " << mp1.Efield_fixed( atm1 ) << std::endl;
	//   TR << "Fixed RF Electric field for res " << res1 << " atom " << rsd1.atom_name( atm1 ) << " at position " << rsd1.xyz( atm1 ) << " is " << mp1.Efield_rf_fixed( atm1 ) << " rkirk " << mp1.rKirkwood( atm1 ) << std::endl;
	//  }
	// }

	PROF_STOP( basic::MULTIPOLE_SETUP );
}

void
MultipoleElecPotential::calculate_induced_fields_for_polarization(
	pose::Pose & pose
) const {
	PROF_START( basic::MULTIPOLE_SETUP );
	Size const nres( pose.size() );

	MultipoleElecPoseInfoOP multipole_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) ) {
		multipole_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );
	} else {
		multipole_info = MultipoleElecPoseInfoOP( new MultipoleElecPoseInfo() );
	}

	// Zero out the induced field members
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd.natoms() ; ++atm1 ) {
			mp1.Efield_induced( atm1 ) = 0.0;
		}
	}

	// Double loop to accumulate E field contributions
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd1( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size res2 = res1; res2 <= nres; ++res2 ) {
			Residue const & rsd2( pose.residue( res2 ) );
			MultipoleElecResidueInfo & mp2( multipole_info->residue_info( res2 ) );
			calculate_res_res_induced_fields_for_polarization( rsd1, mp1, rsd2, mp2 );
		}
	}

	// Print out energy field values
	// for ( Size res1 = 1; res1 <= nres; ++res1 ) {
	//  Residue const & rsd1( pose.residue( res1 ) );
	//  MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
	//  for ( Size atm1 = 1 ; atm1 <= rsd1.natoms() ; ++atm1 ) {
	//   TR << "Induced Electric field for res " << res1 << " atom " << rsd1.atom_name( atm1 ) << " at position " << rsd1.xyz( atm1 ) << " is " << mp1.Efield_induced( atm1 ) << std::endl;
	//  }
	// }

	PROF_STOP( basic::MULTIPOLE_SETUP );
}

void
MultipoleElecPotential::clear_induced_fields(
	pose::Pose & pose
) const {
	PROF_START( basic::MULTIPOLE_SETUP );
	Size const nres( pose.size() );

	MultipoleElecPoseInfoOP multipole_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) ) {
		multipole_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );
	} else {
		multipole_info = MultipoleElecPoseInfoOP( new MultipoleElecPoseInfo() );
	}

	// Zero out the induced field members
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd.natoms() ; ++atm1 ) {
			mp1.Efield_induced( atm1 ) = 0.0;
		}
	}

	PROF_STOP( basic::MULTIPOLE_SETUP );
}

void
MultipoleElecPotential::store_induced_dipoles(
	pose::Pose & pose
) const {
	PROF_START( basic::MULTIPOLE_SETUP );
	Size const nres( pose.size() );

	MultipoleElecPoseInfoOP multipole_info;

	multipole_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );

	// Zero out the induced field members
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd.natoms() ; ++atm1 ) {
			mp1.nonconst_stored_induced_dipole( atm1 ) = mp1.induced_dipole( atm1 );
		}
	}

	PROF_STOP( basic::MULTIPOLE_SETUP );
}

void
MultipoleElecPotential::get_polarization_from_fields(
	pose::Pose & pose
) const {
	PROF_START( basic::MULTIPOLE_SETUP );
	Size const nres( pose.size() );

	MultipoleElecPoseInfoOP multipole_info;

	multipole_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );

	// Zero out the induced field members
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd.natoms() ; ++atm1 ) {
			if ( rsd.is_virtual( atm1 ) ) {
				mp1.nonconst_induced_dipole( atm1 ) = 0.0;
				continue;
			}
			Real const pol( mp1.mp_param( atm1 )->polarity() );
			if ( pol == 0.0 ) {
				mp1.nonconst_induced_dipole( atm1 ) = 0.0;
			} else {
#ifdef NOTDEF
				if( mp1.Efield_fixed( atm1 ) != mp1.Efield_fixed( atm1 ) ) {
					TR << "Bad field fixed for " << res1 << " atom " << rsd.atom_name( atm1 ) << std::endl;
				}
				if( mp1.Efield_rf_fixed( atm1 ) != mp1.Efield_rf_fixed( atm1 ) ) {
					TR << "Bad field rf fixed for " << res1 << " atom " << rsd.atom_name( atm1 ) << std::endl;
				}
				if( mp1.Efield_induced( atm1 ) != mp1.Efield_induced( atm1 ) ) {
					TR << "Bad field induced for " << res1 << " atom " << rsd.atom_name( atm1 ) << std::endl;
				}
#endif
				mp1.nonconst_induced_dipole( atm1 ) =  mp1.mp_param( atm1 )->polarity() *
					( mp1.Efield_fixed( atm1 ) + mp1.Efield_rf_fixed( atm1 ) + mp1.Efield_induced( atm1 ) );
			}
		}
	}

#ifdef NOTDEF
	// Print out polarization values
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd1( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd1.natoms() ; ++atm1 ) {
			if ( rsd1.is_virtual( atm1 ) ) continue;
			TR << "Induced polarization for res " << res1 << " atom " << rsd1.atom_name( atm1 ) << " at position " << rsd1.xyz( atm1 ) << " is " << mp1.induced_dipole( atm1 ) << " with polarization " << mp1.mp_param( atm1 )->polarity() << std::endl;
		}
	}
#endif

	PROF_STOP( basic::MULTIPOLE_SETUP );
}

void
MultipoleElecPotential::get_effective_radii(
	pose::Pose & pose
) const {
	PROF_START( basic::MULTIPOLE_SETUP );
	Size const nres( pose.size() );

	MultipoleElecPoseInfoOP multipole_info;

	multipole_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );

	Size const GK_RADIUS_INDEX( pose.residue(1).atom_type_set().extra_parameter_index( "GK_RADIUS" ) );

	Real const third( 1.0/3.0 );
	Real const pi43( 4.0*third*numeric::constants::d::pi );
	Real const hct_scale( 0.69 );
	Real const bondi_scale( 1.03 );

	// Accumulate overlap between atoms

	for ( Size res1 = 1 ; res1 <= nres; ++res1 ) {
		Residue const & rsd1( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd1.natoms() ; ++atm1 ) {
			Real const rsolv1( bondi_scale*rsd1.atom_type( atm1 ).extra_parameter( GK_RADIUS_INDEX ) );
			Real accum( pi43/std::pow( rsolv1, 3.0 ) );
			if ( rsolv1 == 0.0 || rsd1.is_virtual( atm1 ) ) continue;
			for ( Size res2 = 1 ; res2 <= nres; ++res2 ) {
				Residue const & rsd2( pose.residue( res2 ) );
				//MultipoleElecResidueInfo & mp2( multipole_info->residue_info( res2 ) );
				for ( Size atm2 = 1 ; atm2 <= rsd2.natoms() ; ++atm2 ) {
					Real const rsolv2( bondi_scale*rsd2.atom_type( atm2 ).extra_parameter( GK_RADIUS_INDEX ) );
					if ( rsolv2 == 0.0 || rsd2.is_virtual( atm2 ) || ( res1 == res2 && atm1 == atm2 ) ) continue;
					Vector const Rij( rsd2.xyz( atm2 ) - rsd1.xyz( atm1 ) );
					Real const dist( Rij.magnitude() );
					Real const dist_squared( Rij.magnitude_squared() );
					Real const sr2( hct_scale * rsolv2 );
					Real const sr2_sqr( sr2 * sr2 );
					Real lower( 0.0 );
					if ( (rsolv1 + dist) < sr2 ) {
						lower = rsolv1;
						Real const upper( sr2 - dist );
						accum += pi43*(1.0/std::pow(upper,3.0) - 1.0/std::pow(lower, 3.0));
					}
					Real const upper( dist + sr2 );
					if ( (rsolv1 + dist) < sr2 ) {
						lower = sr2 - dist;
					} else if ( dist < (rsolv1 + sr2) ) {
						lower = rsolv1;
					} else {
						lower = dist - sr2;
					}
					Real const lower2 ( lower*lower );
					Real const lower4 ( lower2*lower2 );
					Real const lowerXdist ( lower * dist );
					Real const lower4Xdist ( lower4 * dist );
					Real const upper2 ( upper*upper );
					Real const upper4 ( upper2*upper2 );
					Real const upperXdist ( upper * dist );
					Real const upper4Xdist ( upper4 * dist );
					accum -= numeric::constants::d::pi*(1.0/12.0)*
						( ( 3.0*(dist_squared-sr2_sqr) + 6.0*upper2 - 8.0*upperXdist)/upper4Xdist -
						( 3.0*(dist_squared-sr2_sqr) + 6.0*lower2 - 8.0*lowerXdist)/lower4Xdist );

				}
			}
			accum = std::pow( accum/pi43, third );
			accum = ( accum < 0.0 ? 0.0001 : accum );
			mp1.nonconst_rKirkwood( atm1 ) = 1.0 / accum;
		}
	}

	// Print out effective radii
	// for ( Size res1 = 1; res1 <= nres; ++res1 ) {
	//  Residue const & rsd1( pose.residue( res1 ) );
	//  MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
	//  for ( Size atm1 = 1 ; atm1 <= rsd1.natoms() ; ++atm1 ) {
	//   Real const rsolv1( bondi_scale*rsd1.atom_type( atm1 ).extra_parameter( GK_RADIUS_INDEX ) );
	//   TR << "Effective radii for res " << res1 << " atom " << rsd1.atom_name( atm1 ) << " at position " << rsd1.xyz( atm1 ) << " is " << mp1.rKirkwood( atm1 ) << " from rsolv " << rsolv1 << std::endl;
	//  }
	// }

	PROF_STOP( basic::MULTIPOLE_SETUP );
}


//////////////////////////////////////////////////////////////////////////////

core::Real
MultipoleElecPotential::relax_induced_dipoles(
	pose::Pose & pose,
	core::Real const relax
) const {
	PROF_START( basic::MULTIPOLE_SETUP );
	Size const nres( pose.size() );

	MultipoleElecPoseInfoOP multipole_info;

	multipole_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );

	core::Real max_diff( 0.0 );

	// Copy over the induced dipole field members
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd.natoms() ; ++atm1 ) {
			Vector const old_value( mp1.induced_dipole( atm1 ) );
			mp1.nonconst_induced_dipole( atm1 ) = (1.0-relax)*mp1.stored_induced_dipole( atm1 ) +  relax*mp1.induced_dipole( atm1 );

			core::Real const this_diff( old_value.distance( mp1.induced_dipole( atm1 ) ) );
			if ( this_diff > max_diff ) {
				max_diff = this_diff;
			}
		}
	}

	PROF_STOP( basic::MULTIPOLE_SETUP );
	return max_diff;
}

//////////////////////////////////////////////////////////////////////////////

void
MultipoleElecPotential::induce_polarizable_dipoles(
	pose::Pose & pose
) const {
	PROF_START( basic::MULTIPOLE_SETUP );
	Size const nres( pose.size() );

	//TR << "In induce_polarize_dipoles()" << std::endl;

	// Early exit if polarization not wanted
	if ( !use_polarization ) {
		clear_induced_fields( pose );
		return;
	}

	MultipoleElecPoseInfoOP multipole_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) ) {
		multipole_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );
	} else {
		multipole_info = MultipoleElecPoseInfoOP( new MultipoleElecPoseInfo() );
	}

	// Will need functions that do the following:
	// 1.  Calculate fields due to fixed multipole and store
	// 2.  Calculate fields due to induced dipoles and store
	// 3.  Update all induced dipoles using the sum of the previous two
	// 4.  A crude optimization algorithm to converge iteratively

	// Initialization
	calculate_fixed_fields_for_polarization( pose );
	clear_induced_fields( pose );
	get_polarization_from_fields( pose );
	calculate_induced_fields_for_polarization( pose );

	// Do a simple relaxation for now.  Tinker uses conjugate gradient.

	Real max_diff( 1.0 );
	Real const test_diff( 1.0e-5 );
	Size const max_iters( 50 );
	// SOR factor
	Real const relax( 0.40 );
	Size iter( 0 );
	while ( ( max_diff > test_diff ) && ( iter < max_iters ) ) {
		iter++;
		// Store old values
		store_induced_dipoles( pose );
		// Calculate new values
		calculate_induced_fields_for_polarization( pose );
		get_polarization_from_fields( pose );
		// Relax towards new values from old, returning largest difference
		max_diff = relax_induced_dipoles( pose, relax );
		//  if( ( iter % 10 ) == 1 ) {
		//   TR << "SOR iteration " << iter << " max_diff is " << max_diff << std::endl;
		//  }
	}

	if ( iter == max_iters ) {
		TR << "Warning in induce_polarizable_dipoles:  Exited SOR after " << max_iters <<
			" iterations, without convergence!" << std::endl;
	}

	// Store the induced dipole in the local coord frame
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd( pose.residue( res1 ) );
		MultipoleElecResidueInfo & mp1( multipole_info->residue_info( res1 ) );
		for ( Size atm1 = 1 ; atm1 <= rsd.natoms() ; ++atm1 ) {
			if ( rsd.is_virtual( atm1 ) ) continue;

			Vector const old_value( mp1.induced_dipole( atm1 ) );
			//   Real check_pdamp( mp1.mp_param( atm1 )->pdamp() );
			//   TR << "Induced dipole " << rsd.xyz( atm1 ) << " is " << mp1.induced_dipole( atm1 ) << " pdamp " << check_pdamp << std::endl;
			mp1.nonconst_stored_induced_dipole( atm1 ) = (mp1.local_coord_matrix( atm1 ).inverse())*mp1.induced_dipole( atm1 );
		}
	}

	PROF_STOP( basic::MULTIPOLE_SETUP );

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultipoleElecPotential::setup_for_scoring(
	pose::Pose & pose
) const {
	PROF_START( basic::MULTIPOLE_SETUP );

	// TR << "in setup_for_scoring use_nblist is " << pose.energies().use_nblist() << std::endl;

	MultipoleElecPoseInfoOP multipole_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) ) {
		multipole_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );
	} else {
		//TR << "Allocating fresh multipole info objects" << std::endl;
		multipole_info = MultipoleElecPoseInfoOP( new MultipoleElecPoseInfo() );
		multipole_info->initialize( pose );
		pose.data().set( pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO, multipole_info );
		assign_all_amoeba_types( pose );
		align_multipole_axes( pose );
		determine_polarization_groups( pose );
		get_effective_radii( pose );
		induce_polarizable_dipoles( pose );
	}

	// TR << "In setup_for_scoring, past initialization" << std::endl;

	// These shouldn't be necessary, but I don't know how to
	// transfer from RotamerSetInfo to ResidueInfo
	// The lines with the extra indent are only necessary because
	// RotamerSetInfo isn't automatically transferred to ResidueInfo
	// after a repack.

	//  multipole_info->initialize( pose );
	//  assign_all_amoeba_types( pose );
	align_multipole_axes( pose );
	//  determine_polarization_groups( pose );
	//  get_effective_radii( pose );
	//  induce_polarizable_dipoles( pose );

	if ( !pose.energies().use_nblist() ) {
		//
		//  //TR << "Getting induced dipoles" << std::endl;
		//  multipole_info->initialize( pose );
		//  assign_all_amoeba_types( pose );
		align_multipole_axes( pose );
		determine_polarization_groups( pose );
		get_effective_radii( pose );
		induce_polarizable_dipoles( pose );
	}

	// TR << "Exiting setup_for_scoring" << std::endl;

	PROF_STOP( basic::MULTIPOLE_SETUP );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief
/// Note: when called at the beginning of rotamer_trials, task.being_packed(i) will be false for all i
/// this ensures that we use all the information we have to compute the current set of radii

void
MultipoleElecPotential::setup_for_packing(
	pose::Pose & , // ellided
	utility::vector1< bool > const & // repacking_residues
) const {
	// ////using core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO;

	// jjh Commenting out for now.  First need to get regular scoring done.  Worry
	// about packing later.
#ifdef NOTDEF
	PROF_START( basic::GB_SETUP_FOR_PACKING );

	MultipoleElecPoseInfoOP mp_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) ) {
		mp_info = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );
	} else {
		mp_info = new MultipoleElecPoseInfo();
	}

	//jjh zero out arrays
	mp_info->initialize( pose );

	/// store info about which positions are moving
	mp_info->set_repack_list( repacking_residues );
	build_placeholders( pose, *mp_info );
	get_template_born_radii( pose, *mp_info );

	pose.data().set( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO, mp_info );

	PROF_STOP( basic::GB_SETUP_FOR_PACKING );
#endif
}

void
MultipoleElecPotential::get_rotamers_multipole_info(
	core::pose::Pose const & pose,
	conformation::RotamerSetBase & rotamer_set
) const {
	MultipoleElecRotamerSetInfoOP mp_info_rotamers( new MultipoleElecRotamerSetInfo( rotamer_set ) );

	for ( Size n=1; n<= rotamer_set.num_rotamers(); ++n ) {
		assign_residue_amoeba_type( *rotamer_set.rotamer(n), mp_info_rotamers->residue_info( n ) );
		align_residue_multipole_axes( pose, *rotamer_set.rotamer(n), mp_info_rotamers->residue_info( n ) );
	}

	rotamer_set.data().set( core::conformation::RotamerSetCacheableDataType::MULTIPOLE_ELEC_ROTAMER_SET_INFO, mp_info_rotamers );
}

void
MultipoleElecPotential::get_rotamers_effective_radii(
	pose::Pose const & pose,
	conformation::RotamerSetBase & rotamer_set
) const {
	using core::conformation::RotamerSetCacheableDataType::MULTIPOLE_ELEC_ROTAMER_SET_INFO;

	MultipoleElecPoseInfoCOP mp_info_pose;

	mp_info_pose = utility::pointer::static_pointer_cast< MultipoleElecPoseInfo const >( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) );

	// this will get cached in the rotamer set
	// this call should initialize the residue_info objects with the appropriate Residue info
	MultipoleElecRotamerSetInfo & mp_info_rotamers
		( rotamer_set.data().get< MultipoleElecRotamerSetInfo >( MULTIPOLE_ELEC_ROTAMER_SET_INFO ) );

	for ( Size n=1; n<= rotamer_set.num_rotamers(); ++n ) {
		get_single_rotamer_effective_radii( *rotamer_set.rotamer(n), pose, mp_info_pose, mp_info_rotamers.residue_info( n ) );
	}
}



void
MultipoleElecPotential::get_single_rotamer_effective_radii(
	Residue const & rsd1,
	pose::Pose const & pose,
	MultipoleElecPoseInfoCOP , //mp_info,
	MultipoleElecResidueInfo & mp1
) const {
	PROF_START( basic::MULTIPOLE_SETUP );
	Size const nres( pose.size() );

	Size const GK_RADIUS_INDEX( pose.residue(1).atom_type_set().extra_parameter_index( "GK_RADIUS" ) );

	Real const third( 1.0/3.0 );
	Real const pi43( 4.0*third*numeric::constants::d::pi );
	Real const hct_scale( 0.69 );
	Real const bondi_scale( 1.03 );

	// Accumulate overlap between atoms
	for ( Size atm1 = 1 ; atm1 <= rsd1.natoms() ; ++atm1 ) {
		Real const rsolv1( bondi_scale*rsd1.atom_type( atm1 ).extra_parameter( GK_RADIUS_INDEX ) );
		Real accum( pi43/std::pow( rsolv1, 3.0 ) );
		if ( rsolv1 == 0.0 || rsd1.is_virtual( atm1 ) ) {
			mp1.nonconst_rKirkwood( atm1 ) = 0.0;
			continue;
		}
		for ( Size res2 = 1 ; res2 <= nres; ++res2 ) {
			// Skip the residue at the same position the rotamer is in.
			// Should replace this with the rotamer itself.
			if ( res2 == rsd1.seqpos() ) continue;
			Residue const & rsd2( pose.residue( res2 ) );
			//MultipoleElecResidueInfo const & mp2( mp_info->residue_info( res2 ) );
			for ( Size atm2 = 1 ; atm2 <= rsd2.natoms() ; ++atm2 ) {
				Real const rsolv2( bondi_scale*rsd2.atom_type( atm2 ).extra_parameter( GK_RADIUS_INDEX ) );
				if ( rsolv2 == 0.0 || rsd2.is_virtual( atm2 ) ) continue;
				Vector const Rij( rsd2.xyz( atm2 ) - rsd1.xyz( atm1 ) );
				Real const dist( Rij.magnitude() );
				Real const dist_squared( Rij.magnitude_squared() );
				Real const sr2( hct_scale * rsolv2 );
				Real const sr2_sqr( sr2 * sr2 );
				Real lower( 0.0 );
				if ( (rsolv1 + dist) < sr2 ) {
					lower = rsolv1;
					Real const upper( sr2 - dist );
					accum += pi43*(1.0/std::pow(upper,3.0) - 1.0/std::pow(lower, 3.0));
				}
				Real const upper( dist + sr2 );
				if ( (rsolv1 + dist) < sr2 ) {
					lower = sr2 - dist;
				} else if ( dist < (rsolv1 + sr2) ) {
					lower = rsolv1;
				} else {
					lower = dist - sr2;
				}

				Real const lower2 ( lower*lower );
				Real const lower4 ( lower2*lower2 );
				Real const lowerXdist ( lower * dist );
				Real const lower4Xdist ( lower4 * dist );
				Real const upper2 ( upper*upper );
				Real const upper4 ( upper2*upper2 );
				Real const upperXdist ( upper * dist );
				Real const upper4Xdist ( upper4 * dist );
#ifdef NOTDEF
				if( dist != dist ) {
					TR << "Problem with dist" << std::endl;
				} else if ( lower != lower ) {
					TR << "Problem with lower" << std::endl;
				} else if ( upper != upper ) {
					TR << "Problem with upper" << std::endl;
				} else if ( sr2_sqr != sr2_sqr ) {
					TR << "Problem with sr2_sqr" << std::endl;
				} else if ( dist_squared != dist_squared ) {
					TR << "Problem with dist_squared" << std::endl;
				}
#endif

				accum -= numeric::constants::d::pi*(1.0/12.0)*
					( ( 3.0*(dist_squared-sr2_sqr) + 6.0*upper2 - 8.0*upperXdist)/upper4Xdist -
					( 3.0*(dist_squared-sr2_sqr) + 6.0*lower2 - 8.0*lowerXdist)/lower4Xdist );

			}
		}
		Real const prior_accum( accum );
		if ( accum < 0.0 ) {
			accum = 0.0001;
		} else {
			accum = std::pow( accum/pi43, third );
		}
		if ( accum != accum ) {
			TR << "Accum fail 2" << std::endl;
			TR << "Prior accum was " << prior_accum << " pi43 is " << pi43 << " third is " << third << std::endl;
		}
		accum = ( accum < 0.0 ? 0.0001 : accum );
		mp1.nonconst_rKirkwood( atm1 ) = 1.0 / accum;
		if ( mp1.nonconst_rKirkwood( atm1 ) != mp1.nonconst_rKirkwood( atm1 ) ) {
			TR << "Problem in rotamer rkirk, accum was " << accum << std::endl;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// called eg after a rotamer substitution is accepted during rotamer trials
void
MultipoleElecPotential::update_residue_for_packing(
	pose::Pose & pose,
	Size const seqpos
) const {
	// ////using core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO;

	MultipoleElecPoseInfo & multipole_info( static_cast< MultipoleElecPoseInfo & >( pose.data().get( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO ) ) );
	MultipoleElecResidueInfo & mp_residue_info( multipole_info.residue_info( seqpos ) );

	Residue const & rsd( pose.residue( seqpos ) );

	mp_residue_info.initialize( rsd );
	assign_residue_amoeba_type( rsd, mp_residue_info );
	align_residue_multipole_axes( pose, rsd, mp_residue_info );
	determine_polarization_groups( pose );
	get_effective_radii( pose );

}

inline void
get_damped_scale_factors(
	MultipoleParameter::MultipoleParameterOP const & mp_param1,
	MultipoleParameter::MultipoleParameterOP const & mp_param2,
	core::Real dist,
	core::Real & scale3,
	core::Real & scale5,
	core::Real & scale7
) {
	// Damping for scale factors
	Real damp( mp_param1->pdamp() *  mp_param2->pdamp() );
	if ( damp == 0.0 ) return;

	Real const pgamma( std::min( mp_param1->thole(), mp_param2->thole() ) );
	damp = -pgamma * std::pow((dist/damp), 3);
	if ( damp > -50.0 ) {
		Real const expdamp( std::exp( damp ) );
		scale3 = 1.0 - expdamp;
		scale5 = 1.0 - expdamp*(1.0-damp);
		scale7 = 1.0 - expdamp*(1.0-damp+0.6*damp*damp);
	}
}

inline void
get_damped_scale_factors_with_derivs(
	MultipoleParameter::MultipoleParameterOP const & mp_param1,
	MultipoleParameter::MultipoleParameterOP const & mp_param2,
	core::Real dist,
	core::Real & scale3,
	core::Real & scale5,
	core::Real & scale7,
	core::Real & dscale3_dr,
	core::Real & dscale5_dr,
	core::Real & dscale7_dr
) {
	// Damping for scale factors
	Real damp( mp_param1->pdamp() *  mp_param2->pdamp() );
	if ( damp == 0.0 ) return;

	Real const pgamma( std::min( mp_param1->thole(), mp_param2->thole() ) );
	damp = -pgamma * std::pow((dist/damp), 3);
	if ( damp > -50.0 ) {
		Real const expdamp( std::exp( damp ) );
		scale3 = 1.0 - expdamp;
		dscale3_dr = -3.0*damp*expdamp/dist;
		scale5 = 1.0 - expdamp*(1.0-damp);
		dscale5_dr = 3.0*damp*damp*expdamp/dist;
		scale7 = 1.0 - expdamp*(1.0-damp+0.6*damp*damp);
		dscale7_dr = -3.0*expdamp*damp*damp*( 0.2 + 0.6*damp )/dist;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline bool
same_polarization_group(
	core::conformation::Residue const & rsd1,
	MultipoleElecResidueInfo const & mp1,
	core::Size atm1,
	core::conformation::Residue const & rsd2,
	core::Size atm2
) {
	if ( abs( static_cast<int>(rsd1.seqpos()) - static_cast<int>(rsd2.seqpos()) ) > 1 ) {
		return false;
	} else {
		utility::vector1< id::AtomID > const & group( mp1.const_my_group( atm1 ) );
		id::AtomID atm2_id( atm2, rsd2.seqpos() );
		if ( std::find( group.begin(), group.end(), atm2_id ) != group.end()  ) {
			return true;
		} else {
			return false;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
MultipoleElecPotential::get_res_res_elecE(
	Residue const & rsd1,
	MultipoleElecResidueInfo const & mp1,
	Residue const & rsd2,
	MultipoleElecResidueInfo const & mp2
) const {
	using namespace etable::count_pair;

	Size natoms1 = rsd1.natoms();
	Size natoms2 = rsd2.natoms();

	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );

	etable::count_pair::CountPairFunctionOP cpfxn( nullptr );

	if ( same_res ) {
		cpfxn = etable::count_pair::CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_34 );
	} else {
		cpfxn = etable::count_pair::CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_34 );
	}

	Real RFelecE = 0.0;
	Real PolelecE = 0.0;
	Real elecE = 0.0;

	// For a description of the polytensor formulation, see
	// the book "The Theory of Intermolecular Forces" by Stone.
	//
	// Following the treatment of Ren and Ponder, "Polarizable
	// Atomic Multipole Water Model for Molecular Mechanics
	// Simulation", J. Phys. Chem. B (2003) v.107 pp 5933-5947

	Real weight( 1.0 );
	Size path_dist( 0 );
	Real const gkc( 2.455 );
	Real const dprotein( Ep );
	Real const dwater( Ew );
	Real const coul( 332.063714 );
	Real const fm( coul * (dprotein-dwater) / dwater );
	Real const fd( coul * 2.0*(dprotein-dwater) / (1.0 + 2.0*dwater) );
	Real const fq( coul * 3.0*(dprotein-dwater) / (2.0 + 3.0*dwater) );

	for ( Size atm1 = 1 ; atm1 <= natoms1 ; ++atm1 ) {
		if ( rsd1.is_virtual( atm1 ) ) continue;

		//  if( !ShouldItCount( rsd1, atm1 ) ) continue;

		Real const q1( mp1.monopole( atm1 ) );
		Vector const & p1( mp1.dipole( atm1 ) );
		Vector const & ip1( mp1.induced_dipole( atm1 ) );
		Matrix const & quad1( mp1.quadrupole( atm1 ) );
		Real const rkirk1( mp1.rKirkwood( atm1 ) );

		for ( Size atm2 = (same_res ? atm1 : 1 ) ; atm2 <= natoms2 ; ++atm2 ) {
			if ( rsd2.is_virtual( atm2 ) ) continue;

			Real atom_atomE(0.0);

			//   if( !ShouldItCount( rsd2, atm2 ) ) continue;

			// Assuming all the checks tell us to calculate the interaction
			Vector rvec( rsd2.xyz( atm2 ) - rsd1.xyz( atm1 ) );
			Real dist2( rvec.magnitude_squared() );

			Real dist = sqrt( dist2 );
			Real const q2( mp2.monopole( atm2 ) );
			Matrix const & quad2( mp2.quadrupole( atm2 ) );
			Real const rkirk2( mp2.rKirkwood( atm2 ) );
			Vector const & p2( mp2.dipole( atm2 ) );
			Vector const & ip2( mp2.induced_dipole( atm2 ) );

			// Intermediate vectors
			Vector int1( mp1.quadrupole( atm1 ) * rvec );
			Vector int2( mp2.quadrupole( atm2 ) * rvec );

			Real const sc2( p1.dot( p2 ) );
			Real const sc3( p1.dot( rvec ) );
			Real const sc4( p2.dot( rvec ) );

			Real const psc2( ip1.dot( p2 ) + p1.dot( ip2 ) );
			Real const psc3( ip1.dot( rvec ) );
			Real const psc4( ip2.dot( rvec ) );
			Real const sc5( int1.dot( rvec ) );
			Real const sc6( int2.dot( rvec ) );
			Real const sc7( int1.dot( p2 ) );
			Real const sc8( int2.dot( p1 ) );
			Real const sc9( int1.dot( int2 ) );
			Real const sc10( 2.0 * ( quad1.xy()*quad2.xy() + quad1.xz()*quad2.xz() + quad1.yz()*quad2.yz() ) +
				quad1.xx()*quad2.xx() + quad1.yy()*quad2.yy() + quad1.zz()*quad2.zz() );

			// Scalar products involving polarization
			Real const psc7( int1.dot( ip2 ) );
			Real const psc8( int2.dot( ip1 ) );


			if ( cpfxn->count( atm1, atm2, weight, path_dist ) &&
					(!same_res || (atm1 != atm2) ) ) {

				Real polar_weight( weight > 0.0 ? 1.0 : 0.0 );
				if ( same_polarization_group( rsd1, mp1, atm1, rsd2, atm2 ) && weight == 0.4 ) {
					//     TR << "ENERGY USING SPECIAL POLAR WEIGHT" << std::endl;
					polar_weight = 0.5;
				}

				// gl functions
				Real const gl0( q1*q2 );
				Real const gl1( q2*sc3 - q1*sc4 + sc2 );
				Real const gl2( q1*sc6 + q2*sc5 - sc3*sc4 + 2.0*( sc7 - sc8 + sc10 ) );
				Real const gl3( sc3*sc6 - sc4*sc5 - 4.0*sc9 );
				Real const gl4( sc5*sc6 );

				// inverse distances
				Real const rr1 = 1.0 / dist;
				Real const rr3 = rr1 / dist2;
				Real const rr5 = 3.0 * rr3 / dist2;
				Real const rr7 = 5.0 * rr5 / dist2;
				Real const rr9 = 7.0 * rr7 / dist2;

				// Damping for scale factors
				Real scale3( 1.0 );
				Real scale5( 1.0 );
				Real scale7( 1.0 );
				get_damped_scale_factors( mp1.mp_param( atm1 ), mp2.mp_param( atm2 ),
					dist, scale3, scale5, scale7 );

				// gl functions involving polarization
				Real const pgl1( q2*psc3 - q1*psc4 + psc2 );
				Real const pgl2( -psc3*sc4 - sc3*psc4 + 2.0*( psc7 - psc8 ) );
				Real const pgl3( psc3*sc6 - psc4*sc5 );

				Real const this_energy( gl0*rr1 + gl1*rr3 + gl2*rr5 + gl3*rr7 + gl4*rr9 );
				Real const this_polarization_energy( scale3*pgl1*rr3 + scale5*pgl2*rr5 + scale7*pgl3*rr7 );

				//    TR << rsd1.xyz( atm1 ) << " and " << rsd2.xyz( atm2 ) << "\t" << 0.5*this_polarization_energy*polar_weight*332.063714 << "  " << polar_weight << std::endl;

				//    TR << "Interaction between " << rsd1.xyz( atm1 ) << " and " << rsd2.xyz( atm2 ) << std::endl;
				//    TR << "Energy " << this_energy << " will be weighted by " << weight << std::endl;
				//TR << "Polarization " << this_polarization_energy << " will be weighted by " << 0.5 * polar_weight << std::endl;

				PolelecE += 0.5 * polar_weight * this_polarization_energy;

				atom_atomE = 332.063714 * ( weight * this_energy + 0.5 * polar_weight * this_polarization_energy );
				//    atom_atomE = 332.063714 * ( 0.5 * polar_weight * this_polarization_energy );
				//    TR << rsd1.xyz( atm1 ) << " and " << rsd2.xyz( atm2 ) << "\t" << 0.5*this_polarization_energy*polar_weight*332.063714 << std::endl;
			}

			if ( use_gen_kirkwood ) {
				if ( rkirk1 == 0.0 || rkirk2 == 0.0 ) continue;
				Real const rs1_rs2( rkirk1 * rkirk2 );
				Real const expterm( std::exp( -dist2/( gkc* rs1_rs2 ) ) );
				Real const expc( expterm / gkc );
				Real const dexpc( -2.0/(gkc*rs1_rs2) );
				Real const gf2( 1.0 / (dist2 + rs1_rs2*expterm ) );
				Real const gf( std::sqrt( gf2 ) );
				Real const gf3( gf2 * gf );
				Real const gf5( gf3 * gf2 );
				Real const gf7( gf5 * gf2 );
				Real const gf9( gf7 * gf2 );

				Real const expc1( 1.0 - expc );

				Real const expcdexpc( -expc * dexpc );

				//Real const xr( rvec.x() );
				//Real const yr( rvec.y() );
				//Real const zr( rvec.z() );
				//    Real const xr2( xr*xr );
				//    Real const yr2( yr*yr );
				//    Real const zr2( zr*zr );

				Real const rf_term_m_m( q1 * q2 * fm * gf );

				Real const rf_term_m_d( -0.5* (
					(fm*expc1+fd)*gf3*q1*sc4 - (fm*expc1+fd)*gf3*q2*sc3 +
					0.5*( (fm*expc1+fd)*gf3*q1*psc4 - (fm*expc1+fd)*gf3*q2*psc3 )));

				Real const rf_term_d_d( fd*( +gf3*sc2 + 0.5*gf3*psc2
					- 3.0*gf5*expc1*sc3*sc4
					- 1.5*gf5*expc1*sc3*psc4
					- 1.5*gf5*expc1*psc3*sc4 ) );

				Real const rf_term_m_q( 0.5*(
					q1*sc6*(3.0*fq*gf5 + 3.0*expc1*expc1*fm*gf5 - fm*expcdexpc*gf3) +
					q2*sc5*(3.0*fq*gf5 + 3.0*expc1*expc1*fm*gf5 - fm*expcdexpc*gf3)
					) );

				Real const rf_term_d_q( 0.5*(
					( sc7 - sc8 )*( 6.0*fd*gf5*expc1 + 6.0*fq*gf5 )  +
					( sc4*sc5 - sc3*sc6 )*( +3.0*fd*gf5*expcdexpc - 15.0*fd*gf7*expc1*expc1 -15.0*fq*gf7*expc1 )
					) );

				Real const rf_term_id_q( 0.25*(
					( psc7 - psc8 )*( 6.0*fd*gf5*expc1 + 6.0*fq*gf5 )  +
					( psc4*sc5 - psc3*sc6 )*( +3.0*fd*gf5*expcdexpc - 15.0*fd*gf7*expc1*expc1 -15.0*fq*gf7*expc1 )
					) );

				Real const rf_term_q_q(
					fq*sc5*sc6*( 105.0*gf9*expc1*expc1 - 15.0*gf7*expcdexpc )   // quads fully dotted with rvec
					-60.0*fq*sc9*gf7*expc1      // product with quads touching
					+ 6.0*fq*sc10*gf5   // element-wise product
				);

				Real rf_energy( rf_term_m_m + rf_term_m_d + rf_term_d_d + rf_term_m_q + rf_term_d_q + rf_term_id_q + rf_term_q_q );

				if ( same_res && (atm1 == atm2 ) ) {
					rf_energy *= 0.5;
				}

				//    TR << "RFE " << rsd1.xyz(atm1) << " " << rsd2.xyz(atm2) << " " << rf_energy << std::endl;

				atom_atomE += rf_energy;
				RFelecE += rf_energy;
			}

			elecE += atom_atomE;
		}
	}

	// TR << "res-res energy between " << rsd1.seqpos() << " and " << rsd2.seqpos() << " is " << elecE << std::endl;
	// TR << "RF res-res energy between " << rsd1.seqpos() << " and " << rsd2.seqpos() << " is " << RFelecE << std::endl;
	// TR << "Pol res-res energy between " << rsd1.seqpos() << " and " << rsd2.seqpos() << " is " << PolelecE << std::endl;

	// TR << "res-res interaction between residue " << rsd1.seqpos() << " and " << rsd2.seqpos() << std::endl;
	// TR << "Returning elecE of " << elecE << std::endl;

	return elecE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
MultipoleElecPotential::calculate_res_res_fixed_fields_for_polarization(
	Residue const & rsd1,
	MultipoleElecResidueInfo & mp1,
	Residue const & rsd2,
	MultipoleElecResidueInfo & mp2
) const {
	using namespace etable::count_pair;

	Size natoms1 = rsd1.natoms();
	Size natoms2 = rsd2.natoms();

	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );

	// Members of a group
	// do not polarize each other via their permanent multipole, although
	// they do mutually polarized each other via their induced dipoles.
	// The information determining which other atoms are in the same group
	// as a given atom is in the polarization parameter file.

	Real const gkc( 2.455 );
	Real const dprotein( Ep );
	Real const dwater( Ew );
	Real const fm( (dprotein-dwater) / dwater );
	Real const fd( 2.0*(dprotein-dwater) / (1.0 + 2.0*dwater) );
	Real const fq( 3.0*(dprotein-dwater) / (2.0 + 3.0*dwater) );

	for ( Size atm1 = 1 ; atm1 <= natoms1 ; ++atm1 ) {
		Real const q1( mp1.monopole( atm1 ) );
		if ( rsd1.is_virtual( atm1 ) ) continue;
		Vector const & p1( mp1.dipole( atm1 ) );
		Matrix const & quad1( mp1.quadrupole( atm1 ) );
		Real const rkirk1( mp1.rKirkwood( atm1 ) );
		for ( Size atm2 = (same_res ? atm1 : 1 ) ; atm2 <= natoms2 ; ++atm2 ) {

			// Assuming the checks tell us to calculate the interaction
			Vector rvec( rsd2.xyz( atm2 ) - rsd1.xyz( atm1 ) );
			Real dist2( rvec.magnitude_squared() );

			Real dist = sqrt( dist2 );
			Real const q2( mp2.monopole( atm2 ) );
			if ( rsd2.is_virtual( atm2 ) ) continue;
			Vector const & p2( mp2.dipole( atm2 ) );
			Matrix const & quad2( mp2.quadrupole( atm2 ) );
			Real const rkirk2( mp2.rKirkwood( atm2 ) );

			// This is the direct field from permanent multipoles
			if ( !same_polarization_group( rsd1, mp1, atm1, rsd2, atm2 ) ) {

				if ( !same_res || ( atm1 != atm2 ) ) {

					// Intermediate vectors
					Vector const int1( mp1.quadrupole( atm1 ) * rvec );
					Vector const int2( mp2.quadrupole( atm2 ) * rvec );

					// Scalar products
					Real const sc3( p1.dot( rvec ) );
					Real const sc4( p2.dot( rvec ) );
					Real const sc5( int1.dot( rvec ) );
					Real const sc6( int2.dot( rvec ) );

					// Damping for scale factors
					Real scale3( 1.0 );
					Real scale5( 1.0 );
					Real scale7( 1.0 );
					get_damped_scale_factors( mp1.mp_param( atm1 ), mp2.mp_param( atm2 ),
						dist, scale3, scale5, scale7 );

					// inverse distances
					Real const rr1( 1.0 / dist );
					Real rr3 = rr1 / dist2;
					Real rr5 = 3.0 * rr3 / dist2;
					Real rr7 = 5.0 * rr5 / dist2;
					rr3 *= scale3;
					rr5 *= scale5;
					rr7 *= scale7;

					Vector const Efield_at_1_due_to_2( -rr3*p2 + 2.0*rr5*int2 - ( q2*rr3 - sc4*rr5 + sc6*rr7 )*rvec );
					Vector const Efield_at_2_due_to_1( -rr3*p1 - 2.0*rr5*int1 + ( q1*rr3 + sc3*rr5 + sc5*rr7 )*rvec );

					// Accumulate
					mp1.Efield_fixed( atm1 ) += Efield_at_1_due_to_2;
					mp2.Efield_fixed( atm2 ) += Efield_at_2_due_to_1;
				}
			}

			if ( use_gen_kirkwood ) {

				if ( rkirk1 == 0.0 || rkirk2 == 0.0 ) continue;

				// This is the reaction field due to permanent multipoles
				Real const rs1_rs2( rkirk1 * rkirk2 );
				Real const expterm( std::exp( -dist2/( gkc* rs1_rs2 ) ) );
				Real const expc( expterm / gkc );
				Real const dexpc( -2.0/(gkc*rs1_rs2) );
				//TR << "expc " << expc << " dexpc " << dexpc << " rk1 " << rkirk1 << " rk2 " << rkirk2 << std::endl;
				Real const gf2( 1.0 / (dist2 + rs1_rs2*expterm ) );
				Real const gf( std::sqrt( gf2 ) );
				Real const gf3( gf2 * gf );
				Real const gf5( gf3 * gf2 );
				Real const gf7( gf5 * gf2 );
				Real A10( -gf3 );
				Real A20( 3.0*gf5 );
				Real A30( -15.0*gf7 );

				Real const expc1( 1.0 - expc );
				Real A01( expc1 * A10 );
				Real A11( expc1 * A20 );
				Real A21( expc1 * A30 );

				Real const expcdexpc( -expc * dexpc );
				Real A12( expc1*A21 + expcdexpc*A20 );

				A01 *= fm;
				A10 *= fd; A11 *= fd; A12 *= fd;
				A20 *= fq; A21 *= fq;

				Real const xr( rvec.x() );
				Real const yr( rvec.y() );
				Real const zr( rvec.z() );
				Real const xr2( xr*xr );
				Real const yr2( yr*yr );
				Real const zr2( zr*zr );

				Real const gux1 = xr * A10;
				Real const guy1 = yr * A10;
				Real const guz1 = zr * A10;

				Real const gc2 = xr * A01;
				Real const gc3 = yr * A01;
				Real const gc4 = zr * A01;
				Real const gux2 = A10 + xr2*A11;
				Real const gux3 = xr * yr * A11;
				Real const gux4 = xr * zr * A11;
				Real const guy2 = gux3;
				Real const guy3 = A10 + yr2*A11;
				Real const guy4 = yr * zr * A11;
				Real const guz2 = gux4;
				Real const guz3 = guy4;
				Real const guz4 = A10 + zr2*A11;
				Real const gqxx2 = xr * (2.0*A20+xr2*A21);
				Real const gqxx3 = yr * xr2*A21;
				Real const gqxx4 = zr * xr2*A21;
				Real const gqyy2 = xr * yr2*A21;
				Real const gqyy3 = yr * (2.0*A20+yr2*A21);
				Real const gqyy4 = zr * yr2 * A21;
				Real const gqzz2 = xr * zr2 * A21;
				Real const gqzz3 = yr * zr2 * A21;
				Real const gqzz4 = zr * (2.0*A20+zr2*A21);
				Real const gqxy2 = yr * (A20+xr2*A21);
				Real const gqxy3 = xr * (A20+yr2*A21);
				Real const gqxy4 = zr * xr * yr * A21;
				Real const gqxz2 = zr * (A20+xr2*A21);
				Real const gqxz3 = gqxy4;
				Real const gqxz4 = xr * (A20+zr2*A21);
				Real const gqyz2 = gqxy4;
				Real const gqyz3 = zr * (A20+yr2*A21);
				Real const gqyz4 = yr * (A20+zr2*A21);

				Real const gux5 = xr * (3.0*A11+xr2*A12);
				Real const gux6 = yr * (A11+xr2*A12);
				Real const gux7 = zr * (A11+xr2*A12);
				Real const gux8 = xr * (A11+yr2*A12);
				Real const gux9 = zr * xr * yr * A12;
				Real const gux10 = xr * (A11+zr2*A12);
				Real const guy5 = yr * (A11+xr2*A12);
				Real const guy6 = xr * (A11+yr2*A12);
				Real const guy7 = gux9;
				Real const guy8 = yr * (3.0*A11+yr2*A12);
				Real const guy9 = zr * (A11+yr2*A12);
				Real const guy10 = yr * (A11+zr2*A12);
				Real const guz5 = zr * (A11+xr2*A12);
				Real const guz6 = gux9;
				Real const guz7 = xr * (A11+zr2*A12);
				Real const guz8 = zr * (A11+yr2*A12);
				Real const guz9 = yr * (A11+zr2*A12);
				Real const guz10 = zr * (3.0*A11+zr2*A12);

				Vector f1_dir( 0.0 );
				f1_dir.x() = p2.x()*gux2 + p2.y()*gux3 + p2.z()*gux4
					+ 0.5 * ( q2*gux1 + quad2.xx()*gux5 + quad2.yy()*gux8 + quad2.zz()*gux10
					+ 2.0*(quad2.xy()*gux6+quad2.xz()*gux7 +quad2.yz()*gux9))
					+ 0.5 * ( q2*gc2 + quad2.xx()*gqxx2 + quad2.yy()*gqyy2 + quad2.zz()*gqzz2
					+ 2.0*(quad2.xy()*gqxy2+quad2.xz()*gqxz2 +quad2.yz()*gqyz2));

				f1_dir.y() = p2.x()*guy2 + p2.y()*guy3 + p2.z()*guy4
					+ 0.5 * ( q2*guy1 + quad2.xx()*guy5 + quad2.yy()*guy8 + quad2.zz()*guy10
					+ 2.0*(quad2.xy()*guy6+quad2.xz()*guy7 +quad2.yz()*guy9))
					+ 0.5 * ( q2*gc3 + quad2.xx()*gqxx3 + quad2.yy()*gqyy3 + quad2.zz()*gqzz3
					+ 2.0*(quad2.xy()*gqxy3+quad2.xz()*gqxz3 +quad2.yz()*gqyz3));

				f1_dir.z() = p2.x()*guz2 + p2.y()*guz3 + p2.z()*guz4
					+ 0.5 * ( q2*guz1 + quad2.xx()*guz5 + quad2.yy()*guz8 + quad2.zz()*guz10
					+ 2.0*(quad2.xy()*guz6+quad2.xz()*guz7 +quad2.yz()*guz9))
					+ 0.5 * ( q2*gc4 + quad2.xx()*gqxx4 + quad2.yy()*gqyy4 + quad2.zz()*gqzz4
					+ 2.0*(quad2.xy()*gqxy4+quad2.xz()*gqxz4 +quad2.yz()*gqyz4));

				Vector f2_dir( 0.0 );
				f2_dir.x() = p1.x()*gux2 + p1.y()*gux3 + p1.z()*gux4
					- 0.5 * ( q1*gux1 + quad1.xx()*gux5 + quad1.yy()*gux8 + quad1.zz()*gux10
					+ 2.0*(quad1.xy()*gux6+quad1.xz()*gux7 +quad1.yz()*gux9))
					- 0.5 * ( q1*gc2 + quad1.xx()*gqxx2 + quad1.yy()*gqyy2 + quad1.zz()*gqzz2
					+ 2.0*(quad1.xy()*gqxy2+quad1.xz()*gqxz2 +quad1.yz()*gqyz2));

				f2_dir.y() = p1.x()*guy2 + p1.y()*guy3 + p1.z()*guy4
					- 0.5 * ( q1*guy1 + quad1.xx()*guy5 + quad1.yy()*guy8 + quad1.zz()*guy10
					+ 2.0*(quad1.xy()*guy6+quad1.xz()*guy7 +quad1.yz()*guy9))
					- 0.5 * ( q1*gc3 + quad1.xx()*gqxx3 + quad1.yy()*gqyy3 + quad1.zz()*gqzz3
					+ 2.0*(quad1.xy()*gqxy3+quad1.xz()*gqxz3 +quad1.yz()*gqyz3));

				f2_dir.z() = p1.x()*guz2 + p1.y()*guz3 + p1.z()*guz4
					- 0.5 * ( q1*guz1 + quad1.xx()*guz5 + quad1.yy()*guz8 + quad1.zz()*guz10
					+ 2.0*(quad1.xy()*guz6+quad1.xz()*guz7 +quad1.yz()*guz9))
					- 0.5 * ( q1*gc4 + quad1.xx()*gqxx4 + quad1.yy()*gqyy4 + quad1.zz()*gqzz4
					+ 2.0*(quad1.xy()*gqxy4+quad1.xz()*gqxz4 +quad1.yz()*gqyz4));

				if ( same_res && ( atm1 == atm2 ) ) {
					f1_dir *= 0.5;
					f2_dir *= 0.5;
				}

				//TR << "Pair rf " << rsd1.xyz(atm1) << " " << rsd2.xyz(atm2) << " " << f1_dir << " " << f2_dir << std::endl;
				//TR << "Pair rf " << rsd1.xyz(atm1) << " " << rsd2.xyz(atm2) << " " << A01 << " " << A10 << " " << A11 << " "
				//    << A12 << " " << A20 << " " << A21 << std::endl;

				mp1.Efield_rf_fixed( atm1 ) += f1_dir;
				mp2.Efield_rf_fixed( atm2 ) += f2_dir;
			}
		}
	}
}

void
MultipoleElecPotential::calculate_res_res_induced_fields_for_polarization(
	Residue const & rsd1,
	MultipoleElecResidueInfo & mp1,
	Residue const & rsd2,
	MultipoleElecResidueInfo & mp2
) const {
	using namespace etable::count_pair;

	Size natoms1 = rsd1.natoms();
	Size natoms2 = rsd2.natoms();

	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );

	// There is no use of count_pair for the purpose of accumulating
	// the total electric field due to induced dipoles.  Everyone
	// contributes.

	Vector Efield( 0.0 );

	Real const gkc( 2.455 );
	Real const dprotein( Ep );
	Real const dwater( Ew );
	Real const fd( 2.0*(dprotein-dwater) / (1.0 + 2.0*dwater) );

	for ( Size atm1 = 1 ; atm1 <= natoms1 ; ++atm1 ) {
		Vector const & p1( mp1.induced_dipole( atm1 ) );
		Real const rkirk1( mp1.rKirkwood( atm1 ) );
		for ( Size atm2 = (same_res ? atm1 : 1 ) ; atm2 <= natoms2 ; ++atm2 ) {

			Vector rvec( rsd2.xyz( atm2 ) - rsd1.xyz( atm1 ) );
			Real dist2( rvec.magnitude_squared() );

			Real dist = sqrt( dist2 );
			Vector const & p2( mp2.induced_dipole( atm2 ) );
			Real const rkirk2( mp2.rKirkwood( atm2 ) );

			if ( rkirk1 == 0.0 || rkirk2 == 0.0 || rsd1.is_virtual( atm1 ) || rsd2.is_virtual( atm2 ) ) continue;

			if ( !same_res || !(atm1 == atm2) ) {
				// Scalar products
				//Real const sc2( p1.dot( p2 ) );
				Real const sc3( p1.dot( rvec ) );
				Real const sc4( p2.dot( rvec ) );

				// gl functions
				//Real const gl1( sc2 );
				//Real const gl2( -sc3*sc4 );

				// Damping for scale factors
				Real scale3( 1.0 );
				Real scale5( 1.0 );
				Real scale7( 1.0 );
				get_damped_scale_factors( mp1.mp_param( atm1 ), mp2.mp_param( atm2 ),
					dist, scale3, scale5, scale7 );

				// inverse distances
				Real rr1 = 1.0 / dist;
				Real rr3 = rr1 / dist2;
				Real rr5 = 3.0 * rr3 / dist2;
				rr3 *= scale3;
				rr5 *= scale5;

				//Real const this_energy( gl1*rr3 + gl2*rr5 );

				Vector const Efield_at_1_due_to_2( -rr3*p2 + sc4*rr5*rvec );
				Vector const Efield_at_2_due_to_1( -rr3*p1 + sc3*rr5*rvec );

				//   TR << "Interaction between " << rsd1.xyz( atm1 ) << " and " << rsd2.xyz( atm2 ) << std::endl;
				//   TR << "Induced field at 1 is " << Efield_at_1_due_to_2 << std::endl;
				//   TR << "Induced field at 2 is " << Efield_at_2_due_to_1 << std::endl;

				// Accumulate
				mp1.Efield_induced( atm1 ) += Efield_at_1_due_to_2;
				mp2.Efield_induced( atm2 ) += Efield_at_2_due_to_1;
			}

			if ( use_gen_kirkwood ) {
				// This is the reaction field due to permanent multipoles
				Real const rs1_rs2( rkirk1 * rkirk2 );
				Real const expterm( std::exp( -dist2/( gkc* rs1_rs2 ) ) );
				Real const expc( expterm / gkc );
				Real const gf2( 1.0 / (dist2 + rs1_rs2*expterm ) );
				Real const gf( std::sqrt( gf2 ) );
				Real const gf3( gf2 * gf );
				Real const gf5( gf3 * gf2 );
				Real A10( -gf3 );
				Real A20( 3.0*gf5 );

				Real const expc1( 1.0 - expc );
				Real A11( expc1 * A20 );

				A10 *= fd; A11 *= fd;

				Real const xr( rvec.x() );
				Real const yr( rvec.y() );
				Real const zr( rvec.z() );
				Real const xr2( xr*xr );
				Real const yr2( yr*yr );
				Real const zr2( zr*zr );

				Real const gux2 = A10 + xr2*A11;
				Real const gux3 = xr * yr * A11;
				Real const gux4 = xr * zr * A11;
				Real const guy2 = gux3;
				Real const guy3 = A10 + yr2*A11;
				Real const guy4 = yr * zr * A11;
				Real const guz2 = gux4;
				Real const guz3 = guy4;
				Real const guz4 = A10 + zr2*A11;

				Vector f1_dir( 0.0 );
				f1_dir.x() = p2.x()*gux2 + p2.y()*gux3 + p2.z()*gux4;
				f1_dir.y() = p2.x()*guy2 + p2.y()*guy3 + p2.z()*guy4;
				f1_dir.z() = p2.x()*guz2 + p2.y()*guz3 + p2.z()*guz4;

				Vector f2_dir( 0.0 );
				f2_dir.x() = p1.x()*gux2 + p1.y()*gux3 + p1.z()*gux4;
				f2_dir.y() = p1.x()*guy2 + p1.y()*guy3 + p1.z()*guy4;
				f2_dir.z() = p1.x()*guz2 + p1.y()*guz3 + p1.z()*guz4;

				if ( same_res && ( atm1 == atm2 ) ) {
					f1_dir *= 0.5;
					f2_dir *= 0.5;
				}

				// Accumulate
				mp1.Efield_induced( atm1 ) += f1_dir;
				mp2.Efield_induced( atm2 ) += f2_dir;
			}
		}
	}

	return;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultipoleElecPotential::calculate_and_store_all_derivs(
	pose::Pose const & pose
) const {
	using namespace etable::count_pair;

	// TR << "Precalculating all of the derivatives" << std::endl;

	MultipoleElecPoseInfo const & multipole_info( static_cast< MultipoleElecPoseInfo const & >( pose.data().get( core::pose::datacache::CacheableDataType::MULTIPOLE_POSE_INFO)));

	// Make sure the cache for the derivs is correctly sized
	cached_atom_derivs_.resize( pose.size() );
	for ( Size ires = 1 ; ires <= pose.size() ; ++ires ) {
		Size const num_atoms( pose.residue( ires ).natoms() );
		cached_atom_derivs_[ ires ].resize( num_atoms );
		for ( Size iat = 1 ; iat <= num_atoms ; ++iat ) {
			cached_atom_derivs_[ ires ][ iat ].f1() = 0.0;
			cached_atom_derivs_[ ires ][ iat ].f2() = 0.0;
		}
	}

	// Now loop over all residue-residue combos
	for ( Size resi = 1 ; resi <= pose.size() ; ++resi ) {
		for ( Size resj = resi ; resj <= pose.size() ; ++resj ) {

			Residue const & rsd1( pose.residue( resi ) );
			Residue const & rsd2( pose.residue( resj ) );

			MultipoleElecResidueInfo const & mp1( multipole_info.residue_info( resi ) );
			MultipoleElecResidueInfo const & mp2( multipole_info.residue_info( resj ) );

			assert( pose.energies().use_nblist() );

			bool const same_res( resi == resj );

			// TR << "Calculating deriv for " << resi << " and " << resj << std::endl;

			CountPairFunctionOP cpfxn( nullptr );
			if ( same_res ) {
				cpfxn = etable::count_pair::CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_34 );
			} else {
				cpfxn = etable::count_pair::CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_34 );
			}

			for ( Size atomi = 1 ; atomi <= rsd1.natoms() ; ++ atomi ) {

				// if ( rsd1.is_virtual( atomi ) ) continue;

				//    if( !ShouldItCount( rsd1, atomi ) ) continue;

				Vector const & xyzi( rsd1.xyz( atomi ) );
				Real const qi( mp1.monopole( atomi ) );
				Vector const & p1( mp1.dipole( atomi ) );
				Vector const & ip1( mp1.induced_dipole( atomi ) );
				Matrix const & quadi( mp1.quadrupole( atomi ) );

				Size const start_atom( same_res ? atomi : 1 );
				for ( Size atomj=start_atom, atomj_end=rsd2.natoms(); atomj<= atomj_end; ++atomj ) {

					// if ( rsd2.is_virtual( atomj ) ) continue;

					//     if( !ShouldItCount( rsd2, atomj ) ) continue;

					Vector const & xyzj( rsd2.xyz(atomj) );
					Vector const & Rij( xyzj - xyzi );
					Real const qj( mp2.monopole( atomj ) );
					Vector const & p2( mp2.dipole( atomj ) );
					Vector const & ip2( mp2.induced_dipole( atomj ) );
					Matrix const & quadj( mp2.quadrupole( atomj ) );
					Real const dist2( Rij.magnitude_squared() );
					Vector const deriv_dr_f1( xyzj.cross_product( xyzi ) );
					Vector const deriv_dr_f2( xyzj - xyzi );
					Real const dist = sqrt( dist2 );

					Vector const p1tot( p1 + 0.5*ip1 );
					Vector const p2tot( p2 + 0.5*ip2 );

					Vector const I1( quadi * Rij );
					Vector const I2( quadj * Rij );
					//Real const I1_dot_p2( I1.dot( p2 ) );
					//Real const I2_dot_p1( I2.dot( p1 ) );
					//Real const I1_dot_ip2( I1.dot( ip2 ) );
					//Real const I2_dot_ip1( I2.dot( ip1 ) );
					Real const S1( I1.dot( Rij ) );
					Real const S2( I2.dot( Rij ) );
					Vector const J1tot( quadi * p2tot );
					Vector const J2tot( quadj * p1tot );
					Vector const J1( quadi * p2 );
					Vector const J2( quadj * p1 );
					Vector const Ji1( quadi * ip2 );
					Vector const Ji2( quadj * ip1 );
					Vector const K1( quadi * I2 );
					Vector const K2( quadj * I1 );

					Real const pi_dot_r( p1.dot_product( Rij ) );
					Real const pj_dot_r( p2.dot_product( Rij ) );

					Vector Xaxis( 1.0, 0.0, 0.0 );
					Vector Yaxis( 0.0, 1.0, 0.0 );
					Vector Zaxis( 0.0, 0.0, 1.0 );

					Vector const Ox( quadj * Xaxis );
					Vector const Oy( quadj * Yaxis );
					Vector const Oz( quadj * Zaxis );
					Vector const Px( quadi * Xaxis );
					Vector const Py( quadi * Yaxis );
					Vector const Pz( quadi * Zaxis );
					Vector const Rx( quadj * Px );
					Vector const Ry( quadj * Py );
					Vector const Rz( quadj * Pz );
					Vector const Ux( quadi * Ox );
					Vector const Uy( quadi * Oy );
					Vector const Uz( quadi * Oz );

					Real const T( 2.0 * ( quadi.xy()*quadj.xy() + quadi.xz()*quadj.xz() + quadi.yz()*quadj.yz() ) +
						quadi.xx()*quadj.xx() + quadi.yy()*quadj.yy() + quadi.zz()*quadj.zz() );


					Real cp_weight( 1.0 );
					Size path_dist( 0 );
					if ( cpfxn->count( atomi, atomj, cp_weight, path_dist ) ) {

						Real polar_weight( cp_weight > 0.0 ? 1.0 : 0.0 );
						if ( same_polarization_group( rsd1, mp1, atomi, rsd2, atomj ) && cp_weight == 0.4 ) {
							//      TR << "ENERGY USING SPECIAL POLAR WEIGHT" << std::endl;
							polar_weight = 0.5;
						}

						Real const inv_dist( 1.0 / dist );
						Real const inv_dist_2( inv_dist * inv_dist );
						Real const inv_dist_3( inv_dist / dist2 );
						Real const inv_dist_4( inv_dist_2 / dist2 );
						Real const inv_dist_5( inv_dist_3 / dist2 );
						Real const inv_dist_6( inv_dist_4 / dist2 );
						Real const inv_dist_7( inv_dist_5 / dist2 );
						Real const inv_dist_8( inv_dist_6 / dist2 );
						Real const inv_dist_9( inv_dist_7 / dist2 );
						Real const inv_dist_11( inv_dist_9 / dist2 );
						Vector const pi_cross_pj( p1.cross_product( p2 ) );
						Vector const pi_cross_r( p1.cross_product( Rij ) );
						Vector const pj_cross_r( p2.cross_product( Rij ) );
						Vector const pi_cross_rj( p1.cross_product( xyzj ) );
						Vector const pj_cross_rj( p2.cross_product( xyzj ) );
						Vector const pi_cross_ri( p1.cross_product( xyzi ) );
						Vector const pj_cross_ri( p2.cross_product( xyzi ) );
						Real const pi_dot_pj( p1.dot_product( p2 ) );

						Vector const ipi_cross_pj( ip1.cross_product( p2 ) );
						Vector const pi_cross_ipj( p1.cross_product( ip2 ) );
						Vector const ipi_cross_r( ip1.cross_product( Rij ) );
						Vector const ipj_cross_r( ip2.cross_product( Rij ) );
						Vector const ipi_cross_rj( ip1.cross_product( xyzj ) );
						Vector const ipj_cross_rj( ip2.cross_product( xyzj ) );
						Vector const ipi_cross_ri( ip1.cross_product( xyzi ) );
						Vector const ipj_cross_ri( ip2.cross_product( xyzi ) );
						Real const ipi_dot_pj( ip1.dot_product( p2 ) );
						Real const pi_dot_ipj( p1.dot_product( ip2 ) );
						Real const ipi_dot_r( ip1.dot_product( Rij ) );
						Real const ipj_dot_r( ip2.dot_product( Rij ) );

						// Damping for scale factors
						Real scale3( 1.0 );
						Real scale5( 1.0 );
						Real scale7( 1.0 );
						Real dscale3_dr( 0.0 );
						Real dscale5_dr( 0.0 );
						Real dscale7_dr( 0.0 );
						get_damped_scale_factors_with_derivs( mp1.mp_param( atomi ), mp2.mp_param( atomj ),
							dist, scale3, scale5, scale7, dscale3_dr, dscale5_dr, dscale7_dr );

						///******** Begin correct derivatives
						///*** You can comment them out, but never, ever alter them
						Real const beta( 15.0*inv_dist_6*pi_dot_r*pj_dot_r - 3.0*inv_dist_4*pi_dot_pj );

						// monopole-monopole interaction derivative
						Vector const m_m_f1( qi * qj * inv_dist_3 * deriv_dr_f1 );
						Vector const m_m_f2( qi * qj * inv_dist_3 * deriv_dr_f2 );

						// monopole-dipole interaction derivative
						Vector const m_d_f1(  qi*inv_dist_3*pj_cross_ri - 3.0*qi*inv_dist_5*pj_dot_r*deriv_dr_f1 );
						Vector const m_d_f2(  qi*inv_dist_3*p2 - 3.0*qi*inv_dist_5*pj_dot_r*deriv_dr_f2 );

						// monopole-induced dipole interaction derivative
						Real const m_id_energy( -qi*ipj_dot_r*inv_dist_4 ); // Includes extra dist factor
						Vector const m_id_f1(  scale3*qi*inv_dist_3*ipj_cross_ri -
							(3.0*scale3*qi*inv_dist_5*ipj_dot_r )*deriv_dr_f1 -
							dscale3_dr*m_id_energy*deriv_dr_f1 );
						Vector const m_id_f2(  scale3*qi*inv_dist_3*ip2 -
							(3.0*scale3*qi*inv_dist_5*ipj_dot_r )*deriv_dr_f2 -
							dscale3_dr*m_id_energy*deriv_dr_f2 );

						// dipole-monopole interaction derivative
						Vector const d_m_f1( -qj*inv_dist_3*pi_cross_rj + 3.0*qj*inv_dist_5*pi_dot_r*deriv_dr_f1 );
						Vector const d_m_f2( -qj*inv_dist_3*p1 + 3.0*qj*inv_dist_5*pi_dot_r*deriv_dr_f2 );

						// induced dipole-monopole interaction derivative
						Real const id_m_energy( qj*ipi_dot_r*inv_dist_4 );  // Includes extra dist factor
						Vector const id_m_f1( -scale3*qj*inv_dist_3*ipi_cross_rj +
							3.0*scale3*qj*inv_dist_5*ipi_dot_r*deriv_dr_f1 -
							dscale3_dr*id_m_energy*deriv_dr_f1 );
						Vector const id_m_f2( -scale3*qj*inv_dist_3*ip1 +
							3.0*scale3*qj*inv_dist_5*ipi_dot_r*deriv_dr_f2 -
							dscale3_dr*id_m_energy*deriv_dr_f2 );

						// dipole-dipole interaction derivative
						Vector const dp_f1( -1.0 * inv_dist_3 * pi_cross_pj - 3.0 * inv_dist_5 * pi_dot_r * pj_cross_r );
						Vector const dp_f2( 0.0 );

						Vector const drvec_f1( 3.0*inv_dist_5*( pi_dot_r * pj_cross_rj + pj_dot_r * pi_cross_rj ) );
						Vector const drvec_f2( 3.0*inv_dist_5*( pi_dot_r * p2 + pj_dot_r * p1 ) );

						Vector const dr_f1( -1.0 * beta * inv_dist * deriv_dr_f1 );
						Vector const dr_f2( -1.0 * beta * inv_dist * deriv_dr_f2 );

						Vector const d_d_f1( dp_f1 + drvec_f1 + dr_f1 );
						Vector const d_d_f2( dp_f2 + drvec_f2 + dr_f2 );


						// dipole-induced dipole interaction derivative
						Real const d_id_beta( 15.0*scale5*inv_dist_6*pi_dot_r*ipj_dot_r - 3.0*scale3*inv_dist_4*pi_dot_ipj );
						Real const d_id_energy_scale3( pi_dot_ipj*inv_dist_4 ); // This includes another dist factor from dr_dalpha
						Real const d_id_energy_scale5( -3.0*pi_dot_r*ipj_dot_r*inv_dist_6 ); // This includes another dist factor from dr_dalpha

						Vector const d_id_dp_f1( -1.0 * scale3 * inv_dist_3 * pi_cross_ipj - 3.0 * scale5 * inv_dist_5 * pi_dot_r * ipj_cross_r );
						Vector const d_id_dp_f2( 0.0 );

						Vector const d_id_drvec_f1( 3.0*scale5 * inv_dist_5*( pi_dot_r * ipj_cross_rj + ipj_dot_r * pi_cross_rj ) );
						Vector const d_id_drvec_f2( 3.0*scale5 * inv_dist_5*( pi_dot_r * ip2 + ipj_dot_r * p1 ) );

						Vector const d_id_dr_f1( -1.0 * d_id_beta * inv_dist * deriv_dr_f1 );
						Vector const d_id_dr_f2( -1.0 * d_id_beta * inv_dist * deriv_dr_f2 );

						Vector const d_id_dscale3_f1( -1.0 * dscale3_dr * d_id_energy_scale3 * deriv_dr_f1 );
						Vector const d_id_dscale3_f2( -1.0 * dscale3_dr * d_id_energy_scale3 * deriv_dr_f2 );

						Vector const d_id_dscale5_f1( -1.0 * dscale5_dr * d_id_energy_scale5 * deriv_dr_f1 );
						Vector const d_id_dscale5_f2( -1.0 * dscale5_dr * d_id_energy_scale5 * deriv_dr_f2 );

						Vector const d_id_f1( d_id_dp_f1 + d_id_drvec_f1 + d_id_dr_f1 + d_id_dscale3_f1 + d_id_dscale5_f1 );
						Vector const d_id_f2( d_id_dp_f2 + d_id_drvec_f2 + d_id_dr_f2 + d_id_dscale3_f2 + d_id_dscale5_f2 );

						// induced dipole-dipole interaction derivative
						Real const id_d_beta( 15.0*scale5*inv_dist_6*ipi_dot_r*pj_dot_r - 3.0*scale3*inv_dist_4*ipi_dot_pj );
						Real const id_d_energy_scale3( ipi_dot_pj*inv_dist_4 );
						Real const id_d_energy_scale5( -3.0*ipi_dot_r*pj_dot_r*inv_dist_6 );

						Vector const id_d_dp_f1( -1.0 * scale3 * inv_dist_3 * ipi_cross_pj - 3.0 * scale5 * inv_dist_5 * ipi_dot_r * pj_cross_r );
						Vector const id_d_dp_f2( 0.0 );

						Vector const id_d_drvec_f1( 3.0*scale5 * inv_dist_5*( ipi_dot_r * pj_cross_rj + pj_dot_r * ipi_cross_rj ) );
						Vector const id_d_drvec_f2( 3.0*scale5 * inv_dist_5*( ipi_dot_r * p2 + pj_dot_r * ip1 ) );

						Vector const id_d_dr_f1( -1.0 * id_d_beta * inv_dist * deriv_dr_f1 );
						Vector const id_d_dr_f2( -1.0 * id_d_beta * inv_dist * deriv_dr_f2 );

						Vector const id_d_dscale3_f1( -1.0 * dscale3_dr * id_d_energy_scale3 * deriv_dr_f1 );
						Vector const id_d_dscale3_f2( -1.0 * dscale3_dr * id_d_energy_scale3 * deriv_dr_f2 );

						Vector const id_d_dscale5_f1( -1.0 * dscale5_dr * id_d_energy_scale5 * deriv_dr_f1 );
						Vector const id_d_dscale5_f2( -1.0 * dscale5_dr * id_d_energy_scale5 * deriv_dr_f2 );

						Vector const id_d_f1( id_d_dp_f1 + id_d_drvec_f1 + id_d_dr_f1 + id_d_dscale3_f1 + id_d_dscale5_f1 );
						Vector const id_d_f2( id_d_dp_f2 + id_d_drvec_f2 + id_d_dr_f2 + id_d_dscale3_f2 + id_d_dscale5_f2 );

						// Start quadrupole derivatives

						///*********** Quadrupole-charge
						Vector const q_m_f1( -6.0*qj*inv_dist_5*I1.cross_product(xyzj) + 15.0*qj*inv_dist_7*I1.dot(Rij)*deriv_dr_f1 );
						Vector const q_m_f2( -6.0*qj*inv_dist_5*I1 + 15.0*qj*inv_dist_7*I1.dot(Rij)*deriv_dr_f2 );
						///*********** Charge-quadrupole
						Vector const m_q_f1( -6.0*qi*inv_dist_5*I2.cross_product(xyzi) + 15.0*qi*inv_dist_7*I2.dot(Rij)*deriv_dr_f1 );
						Vector const m_q_f2( -6.0*qi*inv_dist_5*I2 + 15.0*qi*inv_dist_7*I2.dot(Rij)*deriv_dr_f2 );

						///*********** Quadrupole-dipole term 1
						Vector const q_d_term1_f1( -6.0*inv_dist_5*I1.cross_product(p2) +
							30.0*inv_dist_7*I1.dot(p2)*deriv_dr_f1 -
							6.0*inv_dist_5*J1.cross_product(xyzj) );
						Vector const q_d_term1_f2( -6.0*inv_dist_5*J1 + 30.0*inv_dist_7*I1.dot(p2)*deriv_dr_f2 );
						///*********** Dipole-quadrupole term 1
						Vector const d_q_term1_f1( -6.0*inv_dist_5*I2.cross_product(p1) -
							30.0*inv_dist_7*I2.dot(p1)*deriv_dr_f1 +
							6.0*inv_dist_5*J2.cross_product(xyzi) );
						Vector const d_q_term1_f2( 6.0*inv_dist_5*J2 - 30.0*inv_dist_7*I2.dot(p1)*deriv_dr_f2 );

						///*********** Quadrupole-dipole term 2
						Vector const q_d_term2_f1( 15.0*inv_dist_7*S1*pj_cross_ri +
							30.0*inv_dist_7*pj_dot_r*I1.cross_product( xyzj ) -
							105.0*inv_dist_9*pj_dot_r*S1*deriv_dr_f1 );
						Vector const q_d_term2_f2( 15.0*inv_dist_7*S1*p2 +
							30.0*inv_dist_7*pj_dot_r*I1 -
							105.0*inv_dist_9*S1*pj_dot_r*deriv_dr_f2 );
						///*********** Dipole-quadrupole term 2
						Vector const d_q_term2_f1( -15.0*S2*inv_dist_7*pi_cross_rj -
							30.0*pi_dot_r*inv_dist_7*I2.cross_product( xyzi ) +
							105.0*pi_dot_r*S2*inv_dist_9*deriv_dr_f1 );
						Vector const d_q_term2_f2( -15.0*S2*inv_dist_7*p1 -
							30.0*pi_dot_r*inv_dist_7*I2 +
							105.0*pi_dot_r*S2*inv_dist_9*deriv_dr_f2 );

						///*********** Quadrupole-induced dipole term 1
						Real const q_id_energy_scale5( 6.0*I1.dot(ip2)*inv_dist_6 ); // extra dist factor from dr_dalpha
						Vector const q_id_term1_f1( -6.0*scale5*inv_dist_5*I1.cross_product(ip2) +
							30.0*scale5*inv_dist_7*I1.dot(ip2)*deriv_dr_f1 -
							6.0*scale5*inv_dist_5*Ji1.cross_product(xyzj) -
							dscale5_dr*q_id_energy_scale5*deriv_dr_f1 );
						Vector const q_id_term1_f2( -6.0*scale5*inv_dist_5*Ji1 +
							30.0*scale5*inv_dist_7*I1.dot(ip2)*deriv_dr_f2 -
							dscale5_dr*q_id_energy_scale5*deriv_dr_f2 );
						///*********** Induced dipole-quadrupole term 1
						Real const id_q_energy_scale5( -6.0*I2.dot(ip1)*inv_dist_6 ); // extra dist factor from dr_dalpha
						Vector const id_q_term1_f1( -6.0*scale5*inv_dist_5*I2.cross_product(ip1) -
							30.0*scale5*inv_dist_7*I2.dot(ip1)*deriv_dr_f1 +
							6.0*scale5*inv_dist_5*Ji2.cross_product(xyzi) -
							dscale5_dr*id_q_energy_scale5*deriv_dr_f1 );
						Vector const id_q_term1_f2( 6.0*scale5*inv_dist_5*Ji2 -
							30.0*scale5*inv_dist_7*I2.dot(ip1)*deriv_dr_f2 -
							dscale5_dr*id_q_energy_scale5*deriv_dr_f2 );

						///*********** Quadrupole-induced dipole term 2
						Real const q_id_energy_scale7( -15.0*ipj_dot_r*S1*inv_dist_8 ); // extra dist factor from dr_dalpha
						Vector const q_id_term2_f1( 15.0*scale7*inv_dist_7*S1*ipj_cross_ri +
							30.0*scale7*inv_dist_7*ipj_dot_r*I1.cross_product( xyzj ) -
							105.0*scale7*inv_dist_9*ipj_dot_r*S1*deriv_dr_f1 -
							dscale7_dr*q_id_energy_scale7*deriv_dr_f1 );
						Vector const q_id_term2_f2( 15.0*scale7*inv_dist_7*S1*ip2 +
							30.0*scale7*inv_dist_7*ipj_dot_r*I1 -
							105.0*scale7*inv_dist_9*S1*ipj_dot_r*deriv_dr_f2 -
							dscale7_dr*q_id_energy_scale7*deriv_dr_f2 );
						///*********** Induced dipole-quadrupole term 2
						Real const id_q_energy_scale7( 15.0*ipi_dot_r*S2*inv_dist_8 ); // extra dist factor from dr_dalpha
						Vector const id_q_term2_f1( -15.0*scale7*S2*inv_dist_7*ipi_cross_rj -
							30.0*scale7*ipi_dot_r*inv_dist_7*I2.cross_product( xyzi ) +
							105.0*scale7*ipi_dot_r*S2*inv_dist_9*deriv_dr_f1 -
							dscale7_dr*id_q_energy_scale7*deriv_dr_f1 );
						Vector const id_q_term2_f2( -15.0*scale7*S2*inv_dist_7*ip1 -
							30.0*scale7*ipi_dot_r*inv_dist_7*I2 +
							105.0*scale7*ipi_dot_r*S2*inv_dist_9*deriv_dr_f2 -
							dscale7_dr*id_q_energy_scale7*deriv_dr_f2 );

						///*********** Quadrupole-quadrupole term sc9
						Vector const q_q_term1_f1( 60.0*inv_dist_7*( K1.cross_product( xyzj ) + K2.cross_product( xyzi ) ) -
							420.0*inv_dist_9*I1.dot( I2 )*deriv_dr_f1 +
							60.0*inv_dist_7*I1.cross_product( I2 ) );
						Vector const q_q_term1_f2( 60.0*inv_dist_7*(K1 + K2) -
							420.0*inv_dist_9*I1.dot(I2)*deriv_dr_f2 );

						///*********** Quadrupole-quadrupole term sc5*sc6
						Vector const q_q_term2_f1( -210.0*inv_dist_9*(S1*I2.cross_product(xyzi) + S2*I1.cross_product(xyzj)) +
							945.0*inv_dist_11*S1*S2*deriv_dr_f1 );
						Vector const q_q_term2_f2( -210.0*inv_dist_9*(S1*I2 + S2*I1) +
							945.0*inv_dist_11*S1*S2*deriv_dr_f2 );

						///*********** Quadrupole-quadrupole term sc10
						//Vector const q_q_term3_f1( 30.0*T*inv_dist_7*deriv_dr_f1 +
						//       6.0*inv_dist_5*( Rx.cross_product( Xaxis ) + Ry.cross_product( Yaxis ) + Rz.cross_product( Zaxis ) +
						//       Ox.cross_product( Px ) + Oy.cross_product( Py ) + Oz.cross_product( Pz ) ) );
						Vector const q_q_term3_f1( 30.0*T*inv_dist_7*deriv_dr_f1 +
							12.0*inv_dist_5*( Rx.cross_product( Xaxis ) + Ry.cross_product( Yaxis ) + Rz.cross_product( Zaxis ) ) );
						Vector const q_q_term3_f2( 30.0*T*inv_dist_7*deriv_dr_f2 );

						Vector const f1 = 332.063714 * cp_weight * ( m_m_f1 + m_d_f1 + d_m_f1 + d_d_f1 +
							q_m_f1 + m_q_f1 + q_d_term1_f1 + d_q_term1_f1 + q_d_term2_f1 + d_q_term2_f1 +
							q_q_term1_f1 + q_q_term2_f1 + q_q_term3_f1 ) +
							332.063714 * 0.5 * polar_weight * ( m_id_f1 + id_m_f1 + id_d_f1 +
							d_id_f1 + q_id_term1_f1 + id_q_term1_f1 + q_id_term2_f1 + id_q_term2_f1 );

						Vector const f2 = 332.063714 * cp_weight * ( m_m_f2 + m_d_f2 + d_m_f2 + d_d_f2 +
							q_m_f2 + m_q_f2 + q_d_term1_f2 + d_q_term1_f2 + q_d_term2_f2 + d_q_term2_f2 +
							q_q_term1_f2 + q_q_term2_f2 + q_q_term3_f2 ) +
							332.063714 * 0.5 * polar_weight * ( m_id_f2 + id_m_f2 + id_d_f2 +
							d_id_f2 + q_id_term1_f2 + id_q_term1_f2 + q_id_term2_f2 + id_q_term2_f2 );

#ifdef NOTDEF
						Vector const f1 = 332.063714 * 0.5 * polar_weight * ( m_id_f1 + id_m_f1 + id_d_f1 +
											d_id_f1 + q_id_term1_f1 + id_q_term1_f1 + q_id_term2_f1 + id_q_term2_f1 );

						Vector const f2 = 332.063714 * 0.5 * polar_weight * ( m_id_f2 + id_m_f2 + id_d_f2 +
											d_id_f2 + q_id_term1_f2 + id_q_term1_f2 + q_id_term2_f2 + id_q_term2_f2 );
#endif

						cached_atom_derivs_[ resi ][ atomi ].f1() += f1;
						cached_atom_derivs_[ resi ][ atomi ].f2() += f2;
						cached_atom_derivs_[ resj ][ atomj ].f1() += -1*f1;
						cached_atom_derivs_[ resj ][ atomj ].f2() += -1*f2;

						// Get the atom indices for coordinate determining atoms when frame
						// rotation correction is needed.

						Size const atomi_ref_index( mp1.coord_frame_ref( atomi ).atomno() );
						Size const atomj_ref_index( mp2.coord_frame_ref( atomj ).atomno() );
						Size const atomi_ref_res( mp1.coord_frame_ref( atomi ).rsd() );
						Size const atomj_ref_res( mp2.coord_frame_ref( atomj ).rsd() );

						//      TR << "Ref info 1: res " << atomi_ref_res << " index " << atomi_ref_index << " 2: res " << atomj_ref_res << " index " << atomj_ref_index << std::endl;

						///************ BEGIN correct coordinate correction for simple dipole
						///******* Do not ever change this stuff

						if ( atomi_ref_index != 0 ) {
							// This is the complete dipole coord frame correction derivative up through dipole-dipole terms
							Vector dp_f1_at1( inv_dist_3 * pi_cross_pj - 3.0 * inv_dist_5 * pj_dot_r * pi_cross_r + qj*inv_dist_3*pi_cross_r );
							Vector dq_f1_at1( 0.0 );

							Vector dp_pol_f1_at1( 0.0 );
							Vector dq_pol_f1_at1( 0.0 );

							dp_pol_f1_at1 += ( inv_dist_3 * scale3 * ipi_cross_pj - 3.0 * scale5*inv_dist_5 * pj_dot_r * ipi_cross_r +
								inv_dist_3 * scale3 * pi_cross_ipj - 3.0 * scale5*inv_dist_5 * ipj_dot_r * pi_cross_r +
								scale3*qj*inv_dist_3*ipi_cross_r );
							dq_pol_f1_at1 += ( 0.0 );

							// This is the quadrupole correction through monopole-quadrupole
							// terms q1*sc6 + q2*sc5 in the energy part
							dq_f1_at1 += ( 6.0*qj*inv_dist_5*I1.cross_product(Rij) );

							// This is the quadrupole and dipole corrections for 2.0*(sc7 - sc8 ) above
							dp_f1_at1 += ( 6.0*inv_dist_5*I2.cross_product( p1 ) );
							dq_f1_at1 += ( 6.0*inv_dist_5*( I1.cross( p2 ) + J1.cross_product( Rij ) ) );

							dp_pol_f1_at1 += scale5*( 6.0*inv_dist_5*I2.cross_product( ip1 ) );
							dq_pol_f1_at1 += scale5*( 6.0*inv_dist_5*( I1.cross( ip2 ) + Ji1.cross_product( Rij ) ) );

							// This is the quadrupole and dipole correction for sc3*sc6 - sc4*sc5 above
							dp_f1_at1 += ( 15.0*S2*inv_dist_7*pi_cross_r );
							dq_f1_at1 += ( -30.0*pj_dot_r*inv_dist_7*I1.cross_product( Rij ) );

							dp_pol_f1_at1 += scale7*( 15.0*S2*inv_dist_7*ipi_cross_r );
							dq_pol_f1_at1 += scale7*( -30.0*ipj_dot_r*inv_dist_7*I1.cross_product( Rij ) );

							// This is the quadrupole and dipole correction for sc9
							dq_f1_at1 += ( -60.0*inv_dist_7*( K1.cross_product(Rij) + I1.cross_product(I2) ) );

							// This is the quadrupole correction for sc5*sc6
							dq_f1_at1 += ( 210.0*inv_dist_9*S2*I1.cross_product(Rij) );

							// This is the quadrupole correction for sc10
							dq_f1_at1 += ( 12.0*inv_dist_5*( Ux.cross_product( Xaxis ) +
								Uy.cross_product( Yaxis ) + Uz.cross_product( Zaxis ) ) );

							Vector const coord_frame_f1_at1 = 332.063714 * cp_weight * ( dp_f1_at1 + dq_f1_at1 ) +
								0.5 * 332.063714 * polar_weight * ( dp_pol_f1_at1 + dq_pol_f1_at1 );

							//       Vector const coord_frame_f1_at1 = 0.5 * 332.063714 * polar_weight * scale3*qj*inv_dist_3*ipi_cross_r;

							//       Vector const coord_frame_f1_at1 = 0.5 * 332.063714 * polar_weight * ( dp_pol_f1_at1 + dq_pol_f1_at1 );

							cached_atom_derivs_[ atomi_ref_res ][ atomi_ref_index ].f1() -= coord_frame_f1_at1;
							cached_atom_derivs_[ resi ][ atomi ].f1() += coord_frame_f1_at1;

						}

						if ( atomj_ref_index != 0 ) {
							// This is the complete dipole coord frame correction derivative up through dipole-dipole terms
							Vector dp_f1_at2(   -1.0 * inv_dist_3 * pi_cross_pj - 3.0 * inv_dist_5 * pi_dot_r * pj_cross_r - qi*inv_dist_3*pj_cross_r );
							Vector dq_f1_at2( 0.0 );

							Vector dp_pol_f1_at2( 0.0 );
							Vector dq_pol_f1_at2( 0.0 );

							dp_pol_f1_at2 += ( -scale3*inv_dist_3 * ipi_cross_pj - 3.0 * scale5*inv_dist_5 * pi_dot_r * ipj_cross_r -
								scale3*inv_dist_3 * pi_cross_ipj - 3.0 * scale5*inv_dist_5 * ipi_dot_r * pj_cross_r -
								scale3*qi*inv_dist_3*ipj_cross_r );
							dq_pol_f1_at2 += ( 0.0 );

							// This is the quadrupole correction through monopole-quadrupole
							// terms q1*sc6 + q2*sc5 in the energy part
							dq_f1_at2 += ( 6.0*qi*inv_dist_5*I2.cross_product(Rij) );

							// This is the quadrupole and dipole corrections for 2.0*(sc7 - sc8 ) above
							dp_f1_at2 += ( -6.0*inv_dist_5*I1.cross_product( p2 ) );
							dq_f1_at2 += ( -6.0*inv_dist_5*( I2.cross( p1 ) + J2.cross_product( Rij ) ) );

							dp_pol_f1_at2 += scale5*( -6.0*inv_dist_5*I1.cross_product( ip2 ) );
							dq_pol_f1_at2 += scale5*( -6.0*inv_dist_5*( I2.cross( ip1 ) + Ji2.cross_product( Rij ) ) );

							// This is the quadrupole and dipole correction for sc3*sc6 - sc4*sc5 above
							dp_f1_at2 += ( -15.0*S1*inv_dist_7*pj_cross_r );
							dq_f1_at2 += ( 30.0*pi_dot_r*inv_dist_7*I2.cross_product( Rij ) );

							dp_pol_f1_at2 += scale7*( -15.0*S1*inv_dist_7*ipj_cross_r );
							dq_pol_f1_at2 += scale7*( 30.0*ipi_dot_r*inv_dist_7*I2.cross_product( Rij ) );

							// This is the quadrupole correction for sc9
							dq_f1_at2 += ( -60.0*inv_dist_7*( K2.cross_product(Rij) - I1.cross_product(I2) ) );

							// This is the quadrupole correction for sc5*sc6
							dq_f1_at2 += ( 210.0*inv_dist_9*S1*I2.cross_product(Rij) );

							// This is the quadrupole correction for sc10
							dq_f1_at2 += ( 12.0*inv_dist_5*( Rx.cross_product( Xaxis ) +
								Ry.cross_product( Yaxis ) + Rz.cross_product( Zaxis ) ) );

							Vector const coord_frame_f1_at2 = 332.063714 * cp_weight * ( dp_f1_at2 + dq_f1_at2 ) +
								0.5 * 332.063714 * polar_weight * ( dp_pol_f1_at2 + dq_pol_f1_at2 );

							//       Vector const coord_frame_f1_at2 = 0.5 * 332.063714 * polar_weight * -1.0 * scale3*qi*inv_dist_3*ipj_cross_r;

							//       Vector const coord_frame_f1_at2 = 0.5 * 332.063714 * polar_weight * ( dp_pol_f1_at2 + dq_pol_f1_at2 );

							cached_atom_derivs_[ atomj_ref_res ][ atomj_ref_index ].f1() -= coord_frame_f1_at2;
							cached_atom_derivs_[ resj ][ atomj ].f1() += coord_frame_f1_at2;

						} ///************ END coordinate frame correction
					} // End of count_pair check

					// Add in derivatives for reaction field if using implicit solvent
					// electrostatics.
					// This is troubling.  Each atom has a self energy, but there is no
					// torsion that can decrease it as deduced from energy terms.  The
					// reason is that the self-energy is only dependent on the Kirkwood
					// radius for the atom, and we don't follow through the math to
					// get d(rK)/d_alpha because it would be insanely costly.
					//
					// But I think I'll eventually have to do it anyway.
					//

					if ( use_gen_kirkwood && (!same_res || (atomi != atomj) ) ) {
						//Real const rf_scale( ( same_res and (atm1 == atm2 ) ? 0.5 : 1.0 ) );
						Real const gkc( 2.455 );
						Real const dprotein( Ep );
						Real const dwater( Ew );
						Real const coul( 332.063714 );
						Real const fm( coul * (dprotein-dwater) / dwater );
						Real const fd( coul * 2.0*(dprotein-dwater) / (1.0 + 2.0*dwater) );
						Real const fq( coul * 3.0*(dprotein-dwater) / (2.0 + 3.0*dwater) );
						Real const rkirk1( mp1.rKirkwood( atomi ) );
						Real const rkirk2( mp2.rKirkwood( atomj ) );

						if ( rkirk1 == 0.0 || rkirk2 == 0.0 ) continue;

						Real const rs1_rs2( rkirk1 * rkirk2 );
						Real const expterm( std::exp( -dist2/( gkc* rs1_rs2 ) ) );
						Real const expc( expterm / gkc );
						Real const dexpc( -2.0/(gkc*rs1_rs2) );
						Real const gf2( 1.0 / (dist2 + rs1_rs2*expterm ) );
						Real const gf( std::sqrt( gf2 ) );
						Real const gf3( gf2 * gf );
						Real const gf5( gf3 * gf2 );
						Real const gf7( gf5 * gf2 );
						Real const gf9( gf7 * gf2 );
						Real const gf11( gf9 * gf2 );

						Real const expc1( 1.0 - expc );
						Real const expcdexpc( -expc * dexpc );

						//Real const xr( Rij.x() );
						//Real const yr( Rij.y() );
						//Real const zr( Rij.z() );
						//Real const xr2( xr*xr );
						//Real const yr2( yr*yr );
						//Real const zr2( zr*zr );

						// Reaction field monopole-monopole
						Vector const rf_m_m_f1( fm * qi * qj * gf3 * expc1 * deriv_dr_f1 );
						Vector const rf_m_m_f2( fm * qi * qj * gf3 * expc1 * deriv_dr_f2 );

						// Reaction field monopole-dipole and dipole-monopole
						Real const r_factor( (fm*expc1 + fd)*gf3 );
						Real const d_r_factor( -3.0*(fm*expc1 + fd)*expc1*gf5 + fm*expcdexpc*gf3 );

						Vector const rf_m_d_f1(
							+0.5*(qi*p2tot.dot(Rij) - qj*p1tot.dot(Rij))*d_r_factor*deriv_dr_f1
							-0.5*r_factor*( qi*xyzi.cross( p2tot ) - qj*xyzj.cross( p1tot ) )
						);

						Vector const rf_m_d_f2(
							+0.5*(qi*p2tot.dot(Rij) - qj*p1tot.dot(Rij))*d_r_factor*deriv_dr_f2
							+0.5*r_factor*( qi*p2tot - qj*p1tot )
						);

						Real const dot_prods( p1.dot(p2) + 0.5*ip1.dot(p2) + 0.5*p1.dot(ip2) );
						Real const dot_r_prods( p1.dot(Rij)*p2.dot(Rij) +
							0.5*ip1.dot(Rij)*p2.dot(Rij) +
							0.5*p1.dot(Rij)*ip2.dot(Rij) );
						Real const d_d_dr_part( -5.0*expc1*expc1*gf7 + expcdexpc*gf5 );

						Vector const rf_d_d_f1(
							+3.0*fd*dot_prods*gf5*expc1*deriv_dr_f1 +
							fd*gf3*( p2.cross(p1tot) + 0.5*ip2.cross(p1) ) +
							3.0*fd*dot_r_prods*d_d_dr_part*deriv_dr_f1 -
							3.0*fd*expc1*gf5*( p1tot.dot(Rij)*xyzj.cross(p2) + p2tot.dot(Rij)*xyzj.cross(p1) +
							0.5*p1.dot(Rij)*xyzj.cross(ip2) + 0.5*p2.dot(Rij)*xyzj.cross(ip1) ) -
							3.0*fd*expc1*gf5*( p1tot.dot(Rij)*p2.cross(Rij) + 0.5*p1.dot(Rij)*ip2.cross(Rij) )
						);

						Vector const rf_d_d_f2(
							+3.0*fd*dot_prods*gf5*expc1*deriv_dr_f2 +
							3.0*fd*dot_r_prods*d_d_dr_part*deriv_dr_f2 +
							3.0*fd*expc1*gf5*( p1.dot(Rij)*p2 + p2.dot(Rij)*p1 + 0.5*p1.dot(Rij)*ip2 +
							0.5*ip2.dot(Rij)*p1 + 0.5*ip1.dot(Rij)*p2 + 0.5*p2.dot(Rij)*ip1 )
						);

						Real const R_m_q( 0.5*(3.0*fq*gf5 + 3.0*expc1*expc1*fm*gf5 - fm*expcdexpc*gf3) );
						Real const dR_m_q( 0.5*( -15.0*fq*gf7*expc1 - 15.0*expc1*expc1*fm*gf7*expc1
							+ 6.0*fm*gf5*expc1*expcdexpc + 3.0*fm*expcdexpc*gf5*expc1
							- fm*gf3*expcdexpc*dexpc ) );

						Vector const rf_m_q_f1(
							-(qi*S2 + qj*S1)*dR_m_q*deriv_dr_f1
							+ 2.0*R_m_q*( qj*xyzj.cross(I1) + qi*xyzi.cross(I2) )
						);
						Vector const rf_m_q_f2(
							-(qi*S2 + qj*S1)*dR_m_q*deriv_dr_f2
							- 2.0*R_m_q*( qi*I2 + qj*I1 )
						);

						Real const R1_d_q( 6.0*fd*gf5*expc1 + 6.0*fq*gf5 );
						Real const dR1_d_q( 6.0*( fd*gf5*expcdexpc - 5.0*fd*gf7*expc1*expc1 -5.0*fq*gf7*expc1 ) );

						Real const R2_d_q( 3.0*( fd*gf5*expcdexpc - 5.0*fd*gf7*expc1*expc1 - 5.0*fq*gf7*expc1 ) );
						Real const dR2_d_q( 3.0*( fd*gf5*expcdexpc*dexpc -5.0*fd*expcdexpc*gf7*expc1
							-10.0*fd*gf7*expcdexpc*expc1 + 35.0*fd*gf9*expc1*expc1*expc1
							-5.0*fq*gf7*expcdexpc + 35.0*fq*gf9*expc1*expc1 ) );

						Vector const rf_d_q_f1(
							-0.5*dR1_d_q*( I1.dot(p2tot) - I2.dot(p1tot) )*deriv_dr_f1
							- 0.5*dR2_d_q*( p2tot.dot(Rij)*S1 - p1tot.dot(Rij)*S2 )*deriv_dr_f1
							+ 0.5*R1_d_q*( xyzj.cross(J1tot) - xyzj.cross(J2tot) )
							+ 0.5*R2_d_q*( 2.0*p2tot.dot(Rij)*xyzj.cross(I1) + S1*xyzj.cross(p2tot)
							-2.0*p1tot.dot(Rij)*xyzj.cross(I2) - S2*xyzj.cross(p1tot) )
							+ 0.5*R1_d_q*p2tot.cross(I1) + 0.5*R2_d_q*S1*p2tot.cross(Rij)
							- 0.5*R1_d_q*( J2tot.cross(Rij) + I2.cross(p1tot) )
							- R2_d_q*p1tot.dot(Rij)*I2.cross(Rij) );

						Vector const rf_d_q_f2(
							-0.5*dR1_d_q*( I1.dot(p2tot) - I2.dot(p1tot) )*deriv_dr_f2
							- 0.5*dR2_d_q*( p2tot.dot(Rij)*S1 - p1tot.dot(Rij)*S2 )*deriv_dr_f2
							- 0.5*R1_d_q*( J1tot - J2tot )
							- 0.5*R2_d_q*( 2.0*p2tot.dot(Rij)*I1 + S1*p2tot - 2.0*p1tot.dot(Rij)*I2 - S2*p1tot ) );

						Real const R1_q_q( fq*( 105.0*gf9*expc1*expc1 - 15.0*gf7*expcdexpc ) );
						Real const dR1_q_q( fq*( 210.0*gf9*expc1*expcdexpc - 945.0*expc1*expc1*expc1*gf11
							-15.0*gf7*expcdexpc*dexpc + 105.0*expcdexpc*gf9*expc1 ) );

						Real const R2_q_q( -60.0*fq*gf7*expc1 );
						Real const dR2_q_q( -60.0*fq*( gf7*expcdexpc - 7.0*gf9*expc1*expc1 ) );

						Real const R3_q_q( 6.0*fq*gf5 );
						Real const dR3_q_q( -30.0*fq*gf7*expc1 );

						Vector const rf_q_q_f1(
							-S1*S2*dR1_q_q*deriv_dr_f1
							+ 2.0*R1_q_q*( S2*xyzj.cross( I1 ) + S1*xyzi.cross( I2 ) )
							-dR2_q_q*I1.dot(I2)*deriv_dr_f1
							+R2_q_q*( xyzj.cross( K1 ) + xyzi.cross( K2 ) )
							-R2_q_q*( I1.cross( I2 ) )
							-T*dR3_q_q*deriv_dr_f1
							+ 2.0*R3_q_q*( Rx.cross_product( Xaxis ) + Ry.cross_product( Yaxis ) + Rz.cross_product( Zaxis ) )
						);

						Vector const rf_q_q_f2(
							-S1*S2*dR1_q_q*deriv_dr_f2
							- 2.0*R1_q_q*( S2*I1 + S1*I2 )
							-dR2_q_q*I1.dot(I2)*deriv_dr_f2
							-R2_q_q*( K1 + K2 )
							-T*dR3_q_q*deriv_dr_f2
						);

						Vector const rf_f1( rf_m_m_f1 + rf_m_d_f1 + rf_d_d_f1 + rf_m_q_f1 + rf_d_q_f1 +rf_q_q_f1 );
						Vector const rf_f2( rf_m_m_f2 + rf_m_d_f2 + rf_d_d_f2 + rf_m_q_f2 + rf_d_q_f2 +rf_q_q_f2 );

						cached_atom_derivs_[ resi ][ atomi ].f1() += rf_f1;
						cached_atom_derivs_[ resi ][ atomi ].f2() += rf_f2;
						cached_atom_derivs_[ resj ][ atomj ].f1() += -1.0 * rf_f1;
						cached_atom_derivs_[ resj ][ atomj ].f2() += -1.0 * rf_f2;


						// Add coordinate fixes

						Size const atomi_ref_index( mp1.coord_frame_ref( atomi ).atomno() );
						Size const atomj_ref_index( mp2.coord_frame_ref( atomj ).atomno() );
						Size const atomi_ref_res( mp1.coord_frame_ref( atomi ).rsd() );
						Size const atomj_ref_res( mp2.coord_frame_ref( atomj ).rsd() );

						if ( atomi_ref_index != 0 ) {
							// This is the complete dipole coord frame correction derivative up through dipole-dipole terms
							//     Vector dp_f1_at1( inv_dist_3 * pi_cross_pj - 3.0 * inv_dist_5 * pj_dot_r * pi_cross_r + qj*inv_dist_3*pi_cross_r );
							// monopole-dipole contribution
							Vector dp_f1_at1( -0.5*r_factor*qj*Rij.cross( p1tot ) );
							// dipole-dipole contribution
							dp_f1_at1 += ( fd*gf3*( p1.cross(p2tot) + 0.5*ip1.cross(p2) )
								-3.0*expc1*fd*gf5*( p2tot.dot(Rij)*p1.cross(Rij) + 0.5*(p2.dot(Rij)*ip1.cross(Rij) ) ) );
							Vector dq_f1_at1( 2.0*R_m_q*qj*I1.cross(Rij) );


							dp_f1_at1 -= ( 0.5*( R1_d_q*p1tot.cross(I2) + R2_d_q*S2*p1tot.cross(Rij) ) );
							dq_f1_at1 -= ( -0.5*R1_d_q*( J1tot.cross(Rij) + I1.cross(p2tot) )
								- R2_d_q*p2tot.dot(Rij)*I1.cross(Rij) );

							// quad-quad parts 1 and 2
							dq_f1_at1 += ( 2.0*R1_q_q*S2*I1.cross(Rij) );
							dq_f1_at1 += ( R2_q_q*( K1.cross_product(Rij) + I1.cross_product(I2) ) );

							dq_f1_at1 += ( 2.0*R3_q_q*( Ux.cross_product( Xaxis ) +
								Uy.cross_product( Yaxis ) + Uz.cross_product( Zaxis ) ) );


							Vector const coord_frame_f1_at1 = ( dp_f1_at1 + dq_f1_at1 );

							cached_atom_derivs_[ atomi_ref_res ][ atomi_ref_index ].f1() -= coord_frame_f1_at1;
							cached_atom_derivs_[ resi ][ atomi ].f1() += coord_frame_f1_at1;
						}

						if ( atomj_ref_index != 0 ) {
							// monopole-dipole contribution
							Vector dp_f1_at2( 0.5*r_factor*qi*Rij.cross( p2tot ) );
							// dipole-dipole contribution
							dp_f1_at2 += ( fd*gf3*( p2.cross(p1tot) + 0.5*ip2.cross(p1) )
								-3.0*expc1*fd*gf5*( p1tot.dot(Rij)*p2.cross(Rij) + 0.5*(p1.dot(Rij)*ip2.cross(Rij) ) ) );

							Vector dq_f1_at2( 2.0*R_m_q*qi*I2.cross(Rij) );

							dp_f1_at2 += ( 0.5*( R1_d_q*p2tot.cross(I1) + R2_d_q*S1*p2tot.cross(Rij) ) );
							dq_f1_at2 += ( -0.5*R1_d_q*( J2tot.cross(Rij) + I2.cross(p1tot) )
								- R2_d_q*p1tot.dot(Rij)*I2.cross(Rij) );

							// quad-quad parts 1 and 2
							dq_f1_at2 += ( 2.0*R1_q_q*S1*I2.cross(Rij) );
							dq_f1_at2 += ( R2_q_q*( K2.cross_product(Rij) - I1.cross_product(I2) ) );


							dq_f1_at2 += ( 2.0*R3_q_q*( Rx.cross_product( Xaxis ) +
								Ry.cross_product( Yaxis ) + Rz.cross_product( Zaxis ) ) );

							Vector const coord_frame_f1_at2 = ( dp_f1_at2 + dq_f1_at2 );

							cached_atom_derivs_[ atomj_ref_res ][ atomj_ref_index ].f1() -= coord_frame_f1_at2;
							cached_atom_derivs_[ resj ][ atomj ].f1() += coord_frame_f1_at2;
						} ///************ END coordinate frame correction
					}
				}
			}
			//   for( Size iat = 1 ; iat <= rsd1.natoms() ; ++iat ) {
			//    TR << "deriv for res1 atomno " << iat << " is f1 " << cached_atom_derivs_[1][ iat ] .f1() << " and f2 " << cached_atom_derivs_[1][ iat ].f2() << std::endl;
			//    }
			//
			//    std::exit(1);
		}
	}
}


// This should only get called once per residue.  Currently
// getting called in the intrares_deriv part.
void
MultipoleElecPotential::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & ,
	MultipoleElecResidueInfo const &,
	MultipoleElecResidueInfo const &,
	pose::Pose const & , // provides context
	Real const & factor,
	utility::vector1< DerivVectorPair > & r1_at_derivs,
	utility::vector1< DerivVectorPair > &
) const
{

#ifdef NOTDEF
			TR << "Fetching deriv between " << rsd1.seqpos() << " and " << rsd2.seqpos() << std::endl;
			TR << "Factor is " << factor << std::endl;
#endif

	for ( Size iat = 1 ; iat <= rsd1.natoms() ; ++iat ) {
		r1_at_derivs[ iat ].f1() += factor * cached_atom_derivs_[ rsd1.seqpos() ][ iat ].f1();
		r1_at_derivs[ iat ].f2() += factor * cached_atom_derivs_[ rsd1.seqpos() ][ iat ].f2();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace scoring
} // namespace core


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::MultipoleParameter::MultipoleParameter() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::MultipoleParameter::save( Archive & arc ) const {
	arc( CEREAL_NVP( coord_type_ ) ); // enum core::scoring::MultipoleAxisType
	arc( CEREAL_NVP( atom_type_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( chirality_sign_ ) ); // Real
	arc( CEREAL_NVP( monopole_ ) ); // Real
	arc( CEREAL_NVP( dipole_ ) ); // Vector
	arc( CEREAL_NVP( quadrupole_ ) ); // Matrix
	arc( CEREAL_NVP( polarity_ ) ); // Real
	arc( CEREAL_NVP( thole_ ) ); // Real
	arc( CEREAL_NVP( pdamp_ ) ); // Real
	arc( CEREAL_NVP( group_members_ ) ); // utility::vector1<Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::MultipoleParameter::load( Archive & arc ) {
	arc( coord_type_ ); // enum core::scoring::MultipoleAxisType
	arc( atom_type_ ); // utility::vector1<Size>
	arc( chirality_sign_ ); // Real
	arc( monopole_ ); // Real
	arc( dipole_ ); // Vector
	arc( quadrupole_ ); // Matrix
	arc( polarity_ ); // Real
	arc( thole_ ); // Real
	arc( pdamp_ ); // Real
	arc( group_members_ ); // utility::vector1<Size>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::MultipoleParameter );
CEREAL_REGISTER_TYPE( core::scoring::MultipoleParameter )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::MultipoleElecPoseInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( residue_info_ ) ); // utility::vector1<MultipoleElecResidueInfoOP>
	arc( CEREAL_NVP( placeholder_residue_ ) ); // utility::vector1<ResidueOP>
	arc( CEREAL_NVP( placeholder_info_ ) ); // utility::vector1<MultipoleElecResidueInfoOP>
	arc( CEREAL_NVP( being_packed_ ) ); // utility::vector1<_Bool>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::MultipoleElecPoseInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( residue_info_ ); // utility::vector1<MultipoleElecResidueInfoOP>
	arc( placeholder_residue_ ); // utility::vector1<ResidueOP>
	arc( placeholder_info_ ); // utility::vector1<MultipoleElecResidueInfoOP>
	arc( being_packed_ ); // utility::vector1<_Bool>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::MultipoleElecPoseInfo );
CEREAL_REGISTER_TYPE( core::scoring::MultipoleElecPoseInfo )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::MultipoleElecResidueInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( type_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( coord_frame_ref_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( monopole_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( rKirkwood_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( dipole_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( induced_dipole_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( stored_induced_dipole_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( induced_rf_dipole_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( stored_induced_rf_dipole_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( Efield_fixed_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( Efield_induced_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( Efield_rf_fixed_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( Efield_rf_induced_ ) ); // utility::vector1<Vector>
	arc( CEREAL_NVP( quadrupole_ ) ); // utility::vector1<Matrix>
	arc( CEREAL_NVP( local_coord_matrix_ ) ); // utility::vector1<Matrix>
	arc( CEREAL_NVP( my_group_ ) ); // utility::vector1<utility::vector1<id::AtomID> >
	arc( CEREAL_NVP( my_local_coord_frame_ ) ); // utility::vector1<utility::vector1<id::AtomID> >
	arc( CEREAL_NVP( mp_param_ ) ); // utility::vector1<MultipoleParameterOP>
	arc( CEREAL_NVP( rosetta_res_type_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::MultipoleElecResidueInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( type_ ); // utility::vector1<Size>
	arc( coord_frame_ref_ ); // utility::vector1<Size>
	arc( monopole_ ); // utility::vector1<Real>
	arc( rKirkwood_ ); // utility::vector1<Real>
	arc( dipole_ ); // utility::vector1<Vector>
	arc( induced_dipole_ ); // utility::vector1<Vector>
	arc( stored_induced_dipole_ ); // utility::vector1<Vector>
	arc( induced_rf_dipole_ ); // utility::vector1<Vector>
	arc( stored_induced_rf_dipole_ ); // utility::vector1<Vector>
	arc( Efield_fixed_ ); // utility::vector1<Vector>
	arc( Efield_induced_ ); // utility::vector1<Vector>
	arc( Efield_rf_fixed_ ); // utility::vector1<Vector>
	arc( Efield_rf_induced_ ); // utility::vector1<Vector>
	arc( quadrupole_ ); // utility::vector1<Matrix>
	arc( local_coord_matrix_ ); // utility::vector1<Matrix>
	arc( my_group_ ); // utility::vector1<utility::vector1<id::AtomID> >
	arc( my_local_coord_frame_ ); // utility::vector1<utility::vector1<id::AtomID> >
	arc( mp_param_ ); // utility::vector1<MultipoleParameterOP>
	arc( rosetta_res_type_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::MultipoleElecResidueInfo );
CEREAL_REGISTER_TYPE( core::scoring::MultipoleElecResidueInfo )


/// @brief Default constructor required by cereal to deserialize this class
core::scoring::MultipoleElecRotamerSetInfo::MultipoleElecRotamerSetInfo() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::MultipoleElecRotamerSetInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( residue_info_ ) ); // utility::vector1<MultipoleElecResidueInfoOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::MultipoleElecRotamerSetInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( residue_info_ ); // utility::vector1<MultipoleElecResidueInfoOP>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::MultipoleElecRotamerSetInfo );
CEREAL_REGISTER_TYPE( core::scoring::MultipoleElecRotamerSetInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_MultipoleElecPotential )
#endif // SERIALIZATION


