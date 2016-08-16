// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/BetaTurnDetection.cc
/// @brief determine the presence and type of beta turn at a specific postion in a pose
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

// Unit Headers
#include <protocols/features/BetaTurnDetection.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace features {

using std::string;
using std::map;
using core::Size;
using core::Real;
using core::pose::Pose;
using utility::excn::EXCN_Msg_Exception;
using utility::vector1;

BetaTurnDetection::BetaTurnDetection() :
	utility::pointer::ReferenceCount(), beta_turn_length_( 3 ), beta_turn_distance_cutoff_( 7.0 )
{}

BetaTurnDetection::BetaTurnDetection( BetaTurnDetection const & ) :
	utility::pointer::ReferenceCount(), beta_turn_length_( 3 ), beta_turn_distance_cutoff_( 7.0 )
{}

BetaTurnDetection::~BetaTurnDetection() {}

string
BetaTurnDetection::type_name() const { return "BetaTurnDetection"; }

map< string, string > const & BetaTurnDetection::get_conformation_to_turn_type_map()
{
	// It pisses me off that C++ works this way, but it does. Sergey promises this line will only ever be executed once.
	static map< string, string > * conformation_to_turn_type = 0;

	if ( conformation_to_turn_type == 0 ) {
		conformation_to_turn_type = new map< string, string >;

		// Turn types will be notated thusly: TurnXX[_NUMERAL], where XX is the Ramachandran hash of residues 2 and 3,
		// _NUMERAL will be present if the turn is a classically recognized turn type, a trailing "p" stands for prime.
		// (e.g. a Type I turn will be annotated "TurnAA_I")

		( *conformation_to_turn_type )[ "AA" ] = "TurnAA_I";
		( *conformation_to_turn_type )[ "AB" ] = "TurnAB_VIII";
		( *conformation_to_turn_type )[ "AL" ] = "TurnAL_IX";
		( *conformation_to_turn_type )[ "AE" ] = "TurnAE";

		( *conformation_to_turn_type )[ "BA" ] = "TurnBA";
		( *conformation_to_turn_type )[ "BB" ] = "TurnBB";
		( *conformation_to_turn_type )[ "BL" ] = "TurnBL_II";
		( *conformation_to_turn_type )[ "BE" ] = "TurnBE";

		( *conformation_to_turn_type )[ "LA" ] = "TurnLA_IXp";
		( *conformation_to_turn_type )[ "LB" ] = "TurnLB";
		( *conformation_to_turn_type )[ "LL" ] = "TurnLL_Ip";
		( *conformation_to_turn_type )[ "LE" ] = "TurnLE_VIIIp";

		( *conformation_to_turn_type )[ "EA" ] = "TurnEA_IIp";
		( *conformation_to_turn_type )[ "EB" ] = "TurnEB";
		( *conformation_to_turn_type )[ "EL" ] = "TurnEL";
		( *conformation_to_turn_type )[ "EE" ] = "TurnEE";

		// Well characterized turn types with Cis residues
		( *conformation_to_turn_type )[ "Ba" ] = "TurnCis3_VIa";
		( *conformation_to_turn_type )[ "Bb" ] = "TurnCis3_VIb";

		// Other possible Cis conformations
		( *conformation_to_turn_type )[ "xX" ] = "TurnCis2";
		( *conformation_to_turn_type )[ "xx" ] = "TurnCis2Cis3";
		( *conformation_to_turn_type )[ "Xx" ] = "TurnCis3other";
	}
	return *conformation_to_turn_type;
}

vector1< string > const & BetaTurnDetection::get_valid_ramachandran_hashes()
{
	// It pisses me off that C++ works this way, but it does. Sergey promises this line will only ever be executed once.
	static vector1< string > * valid_ramachandran_hashes = 0;

	if ( valid_ramachandran_hashes == 0 ) {
		valid_ramachandran_hashes = new vector1< string >;
		valid_ramachandran_hashes->resize( number_of_ramachandran_hashes );

		( *valid_ramachandran_hashes )[ A ] = "A";
		( *valid_ramachandran_hashes )[ B ] = "B";
		( *valid_ramachandran_hashes )[ L ] = "L";
		( *valid_ramachandran_hashes )[ E ] = "E";
	}

	return *valid_ramachandran_hashes;
}

bool BetaTurnDetection::all_turn_residues_are_on_the_same_chain( Pose const & pose, Size first_residue ) const
{
	Size chain = pose.residue( first_residue ).chain();
	for ( Size residue_number = first_residue + 1; residue_number <= first_residue + beta_turn_length_; ++residue_number ) {
		if ( pose.residue( first_residue ).chain() != chain ) {
			return false;
		}
	}
	return true;
}


bool BetaTurnDetection::residue_range_is_protein( Pose const & pose, Size range_begin, Size range_end ) const
{
	for ( Size current_residue = range_begin; current_residue <= range_end; ++current_residue ) {
		if ( !pose.residue( current_residue ).is_protein() ) {
			return false;
		}
	}
	return true;
}


bool BetaTurnDetection::beta_turn_present( Pose const & pose, Size first_residue ) const
{
	return ( pose.residue( first_residue ).xyz( "CA" ) -  pose.residue( first_residue + beta_turn_length_ ).xyz( "CA" ) ).norm() <= beta_turn_distance_cutoff_;
}

string const & BetaTurnDetection::beta_turn_type( Pose const & pose, Size first_residue ) const
{
	string rama_hash = determine_ramachandran_hash( pose, first_residue );
	validate_ramachandran_hash( rama_hash );
	return get_conformation_to_turn_type_map().find( rama_hash )->second;
}


string BetaTurnDetection::determine_ramachandran_hash( Pose const & pose, Size first_residue ) const
{
	string rama_hash = "";
	for ( Size residue_number = first_residue + 1; residue_number < first_residue + beta_turn_length_; ++residue_number ) {
		rama_hash += determine_ramachandran_hash_for_residue_with_dihedrals( pose.phi( residue_number ), pose.psi( residue_number ), pose.omega( residue_number ) );
	}
	return rama_hash;
}

/// @brief For the purposes of classifying beta-turns, Ramachandran space has been hashed into four large areas.
/// In most turns the dihedral angles are not close to the boundaries as defined, so this provides a simple
/// way of accurately classifying beta-turns.
///
/// @details The four regions of space are defined as:
/// A: phi <= 0, -100 < psi <= 50
/// B: phi <= 0, psi > 50 OR psi <= -100
/// L: phi > 0, -50 < psi <= 100
/// E: phi > 0, psi > 100 OR psi <= -50
///
/// Note: In the case of a Cis peptide plane, the lowercase letter for the hash will be returned.
///
/// Pictoral representation of Ramachandran hashing used for beta-turn classification:
/// <pre>
///       |----------------------|
///       |          |           |
///       |    B     |     E     |
///       |          |===========| 100
///    50 |==========|           |
///  p    |          |     L     |
///  s  0 |--- A ----------------|
///  i    |          |           |
///       |          |===========| -50
///  -100 |==========|           |
///       |          |     E     |
///       |    B     |          |
///       |----------------------|
///     -180         0          180
///                 phi
/// </pre>
string BetaTurnDetection::determine_ramachandran_hash_for_residue_with_dihedrals( Real phi, Real psi, Real omega ) const
{
	string rama_hash;
	if ( phi <= 0. ) {
		if ( psi > -100. && psi <= 50. ) {
			rama_hash = get_valid_ramachandran_hashes()[ A ];
		} else {
			rama_hash = get_valid_ramachandran_hashes()[ B ];
		}
	} else {
		if ( psi  > -50. && psi <= 100. ) {
			rama_hash = get_valid_ramachandran_hashes()[ L ];
		} else {
			rama_hash = get_valid_ramachandran_hashes()[ E ];
		}
	}

	// Return the lower case letter for the hash for Cis peptide planes
	if ( omega > -90 && omega <= 90 ) {
		transform( rama_hash.begin(), rama_hash.end(), rama_hash.begin(), ::tolower );
	}
	return rama_hash;
}

void BetaTurnDetection::validate_ramachandran_hash( std::string & rama_hash ) const
{
	if ( ! get_conformation_to_turn_type_map().count(  rama_hash ) ) {
		string cis_trans_hash = "";

		for ( string::const_iterator it = rama_hash.begin(); it != rama_hash.end(); ++it ) {
			bool cis_peptide_bond = islower( * it );
			string single_residue_rama_hash( 1, toupper( * it ) );

			if ( ! get_valid_ramachandran_hashes().contains( single_residue_rama_hash ) ) {
				throw EXCN_Msg_Exception( "The Ramachandran hash '" + rama_hash + "' contains '" + string( 1, * it ) + ",' which is not valid. " +
					"Valid Ramachandran hashes are 'A', 'B', 'L' and 'E' for trans peptide bonds, and 'a', 'b', 'l' and 'e' for cis peptide bonds."
				);
			}
			cis_trans_hash += cis_peptide_bond ? "x" : "X";
		}

		if ( ! get_conformation_to_turn_type_map().count(  cis_trans_hash ) ) {
			throw EXCN_Msg_Exception( "The Ramachandran hash '" + rama_hash +
				"' is not recognized as a valid beta-turn type.  " +
				"The attempt to create a generic hash based on the omega dihedral angle resulted in '" +
				cis_trans_hash + ",' which is also not recognized as a valid beta-turn type."
			);
		}

		rama_hash = cis_trans_hash;
	}
}

} // namesapce features
} // namespace protocols
