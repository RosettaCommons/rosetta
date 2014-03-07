// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	protocols/features/BetaTurnDetectionFeatures.cc
/// @brief	report comments stored with each pose
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

// Unit Headers
#include <protocols/features/BetaTurnDetectionFeatures.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

//External

// Platform Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.hh>

// Numeric Headers
#include <numeric/HomogeneousTransform.hh>

// Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>

#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>


// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/foreach.hpp>

// C++ Headers
#include <algorithm>
#include <map>
#include <sstream>
#include <string>

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::transform;
using std::endl;
using std::map;
using basic::database::safely_write_to_database;
using basic::database::safely_prepare_statement;
using core::Size;
using core::SSize;
using core::Real;
using core::pose::Pose;
using core::conformation::Residue;
using basic::datacache::DataMap;
using protocols::filters::Filters_map;
using protocols::moves::Movers_map;
using numeric::HomogeneousTransform;
using numeric::xyzVector;
using utility::sql_database::sessionOP;
using utility::excn::EXCN_Msg_Exception;
using utility::vector1;
using utility::tag::TagCOP;
using cppdb::statement;
using cppdb::result;

BetaTurnDetectionFeatures::BetaTurnDetectionFeatures() :
	FeaturesReporter(), beta_turn_length( 3 ), beta_turn_distance_cutoff( 7.0 )
{}

BetaTurnDetectionFeatures::BetaTurnDetectionFeatures( BetaTurnDetectionFeatures const & ) :
	FeaturesReporter(), beta_turn_length( 3 ), beta_turn_distance_cutoff( 7.0 )
{}

BetaTurnDetectionFeatures::~BetaTurnDetectionFeatures() {}

string
BetaTurnDetectionFeatures::type_name() const { return "BetaTurnDetectionFeatures"; }

void
BetaTurnDetectionFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_beta_turns_table_schema(db_session);
}

void
BetaTurnDetectionFeatures::write_beta_turns_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt());
	Column residue_begin("residue_begin", new DbInteger());
	Column turn_type("turn_type", new DbText());

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(residue_begin);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(residue_begin);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey foreign_key(foreign_key_columns, "residues", reference_columns, true);

	Schema table("beta_turns", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(turn_type);

	table.write(db_session);
}

utility::vector1<std::string>
BetaTurnDetectionFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

/// @details
/// An anchor is a take off and landing for a loop.
/// Every residue in the loop must be relevant in order for the loop to be stored.
Size
BetaTurnDetectionFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
	string beta_turns_stmt_string = "INSERT INTO beta_turns (struct_id, residue_begin, turn_type) VALUES (?,?,?);";
	statement beta_turns_stmt(
		safely_prepare_statement(beta_turns_stmt_string, db_session));

	for(SSize begin=1; begin <= SSize( pose.total_residue() - beta_turn_length ); ++begin){
		Size end = begin + beta_turn_length;

		if ( !check_relevant_residues_range( relevant_residues, begin, end ) || !residue_range_is_protein( pose, begin, end ) || !all_turn_residues_are_on_the_same_chain( pose, begin ) || !beta_turn_present( pose, begin ) )
		{
				continue;
		}

		// Add stuff to database
		beta_turns_stmt.bind(1,struct_id);
		beta_turns_stmt.bind(2,begin);
		beta_turns_stmt.bind(3, beta_turn_type( pose, begin ) );
		basic::database::safely_write_to_database( beta_turns_stmt );

	}
	return 0;
}

map< string, string > const & BetaTurnDetectionFeatures::get_conformation_to_turn_type_map()
{
	// It pisses me off that C++ works this way, but it does. Sergey promises this line will only ever be executed once.
	static map< string, string > * conformation_to_turn_type = 0;
	
	if ( conformation_to_turn_type == 0 )
	{
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

vector1< string > const & BetaTurnDetectionFeatures::get_valid_ramachandran_hashes()
{
	// It pisses me off that C++ works this way, but it does. Sergey promises this line will only ever be executed once.
	static vector1< string > * valid_ramachandran_hashes = 0;
	
	if ( valid_ramachandran_hashes == 0 )
	{
		valid_ramachandran_hashes = new vector1< string >;
		valid_ramachandran_hashes->resize( number_of_ramachandran_hashes );

		( *valid_ramachandran_hashes )[ A ] = "A";
		( *valid_ramachandran_hashes )[ B ] = "B";
		( *valid_ramachandran_hashes )[ L ] = "L";
		( *valid_ramachandran_hashes )[ E ] = "E";		
	}
	
	return *valid_ramachandran_hashes;
}

bool BetaTurnDetectionFeatures::all_turn_residues_are_on_the_same_chain( Pose const & pose, Size first_residue ) const
{
	Size chain = pose.residue( first_residue ).chain();
	for ( Size residue_number = first_residue + 1; residue_number <= first_residue + beta_turn_length; ++residue_number )
	{
		if ( pose.residue( first_residue ).chain() != chain )
		{
			return false;
		}
	}
	return true;
}


bool BetaTurnDetectionFeatures::residue_range_is_protein( Pose const & pose, Size range_begin, Size range_end ) const
{
	for ( Size current_residue = range_begin; current_residue <= range_end; ++current_residue )
	{
		if ( !pose.residue( current_residue ).is_protein() )
		{
				return false;
		}
	}
	return true;
}


bool BetaTurnDetectionFeatures::beta_turn_present( Pose const & pose, Size first_residue ) const
{
	return ( pose.residue( first_residue ).xyz( "CA" ) -  pose.residue( first_residue + beta_turn_length ).xyz( "CA" ) ).norm() <= beta_turn_distance_cutoff;
}

string const & BetaTurnDetectionFeatures::beta_turn_type( Pose const & pose, Size first_residue ) const
{
	string rama_hash = determine_ramachandran_hash( pose, first_residue );
	validate_ramachandran_hash( rama_hash );
	return get_conformation_to_turn_type_map().find( rama_hash )->second;
}


string BetaTurnDetectionFeatures::determine_ramachandran_hash( Pose const & pose, Size first_residue ) const
{
	string rama_hash = "";
	for ( Size residue_number = first_residue + 1; residue_number < first_residue + beta_turn_length; ++residue_number )
	{
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
///
///	      |----------------------|
///	      |	         |           |
///	      |    B     |     E     |
///	      |	         |===========| 100
///    50 |==========|			     |
///	 p    |	         |     L     |
///  s  0 |--- A ----------------|
///	 i    |	         |           |
///	      |	         |===========| -50
///  -100 |==========|			     |
///	      |	         |     E     |
///	      |    B     |	         |
///	      |----------------------|
///	    -180         0          180
///	                phi
string BetaTurnDetectionFeatures::determine_ramachandran_hash_for_residue_with_dihedrals( Real phi, Real psi, Real omega ) const
{	
	string rama_hash;
	if ( phi <= 0. )
	{
		if ( psi > -100. && psi <= 50. )
		{
			rama_hash = get_valid_ramachandran_hashes()[ A ];
		}
		else
		{
			rama_hash = get_valid_ramachandran_hashes()[ B ];
		}
	}
	else
	{
		if ( psi  > -50. && psi <= 100. )
		{
			rama_hash = get_valid_ramachandran_hashes()[ L ];
		}
		else
		{
			rama_hash = get_valid_ramachandran_hashes()[ E ];
		}
	}
	
	// Return the lower case letter for the hash for Cis peptide planes
	if ( omega > -90 && omega <= 90 )
	{
		transform( rama_hash.begin(), rama_hash.end(), rama_hash.begin(), ::tolower );
	}
	return rama_hash;
}

void BetaTurnDetectionFeatures::validate_ramachandran_hash( std::string & rama_hash ) const
{
	if ( ! get_conformation_to_turn_type_map().count(  rama_hash ) )
	{
		string cis_trans_hash = "";
		
		for ( string::const_iterator it = rama_hash.begin(); it != rama_hash.end(); ++it )
		{
			bool cis_peptide_bond = islower( * it );
			string single_residue_rama_hash( 1, toupper( * it ) );
			
			if ( ! get_valid_ramachandran_hashes().contains( single_residue_rama_hash ) )
			{
				throw EXCN_Msg_Exception( "The Ramachandran hash '" + rama_hash + "' contains '" + string( 1, * it ) + ",' which is not valid. " +
					"Valid Ramachandran hashes are 'A', 'B', 'L' and 'E' for trans peptide bonds, and 'a', 'b', 'l' and 'e' for cis peptide bonds."
				);
			}
			cis_trans_hash += cis_peptide_bond ? "x" : "X";
		}
		
		if ( ! get_conformation_to_turn_type_map().count(  cis_trans_hash ) )
		{
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
