// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/BetaTurnDetectionFeatures.cc
/// @brief  report comments stored with each pose
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

// Unit Headers
#include <protocols/features/BetaTurnDetectionFeatures.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

//External
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

// Platform Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <protocols/moves/DataMap.hh>

// Numeric Headers
#include <numeric/HomogeneousTransform.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>

// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <string>
#include <map>
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using std::map;
using basic::database::safely_write_to_database;
using basic::database::safely_prepare_statement;
using core::Size;
using core::SSize;
using core::Real;
using core::pose::Pose;
using core::conformation::Residue;
using protocols::moves::DataMap;
using protocols::filters::Filters_map;
using protocols::moves::Movers_map;
using numeric::HomogeneousTransform;
using numeric::xyzVector;
using utility::sql_database::sessionOP;
using utility::vector1;
using utility::tag::TagPtr;
using cppdb::statement;
using cppdb::result;

bool BetaTurnDetectionFeatures::initialized_( false );
Size const BetaTurnDetectionFeatures::beta_turn_length = 3;
core::Real const BetaTurnDetectionFeatures::beta_turn_distance_cutoff = 7.0;
map< string, string > BetaTurnDetectionFeatures::conformation_to_turn_type_;

BetaTurnDetectionFeatures::BetaTurnDetectionFeatures() :
	FeaturesReporter()
{
	setup_conformation_to_turn_type_map();
}

BetaTurnDetectionFeatures::BetaTurnDetectionFeatures( BetaTurnDetectionFeatures const & ) :
	FeaturesReporter()
{
	setup_conformation_to_turn_type_map();
}

BetaTurnDetectionFeatures::~BetaTurnDetectionFeatures() {}

string
BetaTurnDetectionFeatures::type_name() const { return "BetaTurnDetectionFeatures"; }

string
BetaTurnDetectionFeatures::schema() const {
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);

	if(db_mode == "sqlite3")
	{
		return
			"CREATE TABLE IF NOT EXISTS beta_turns (\n"
			"	struct_id BLOB,\n"
			"	residue_begin INTEGER,\n"
            "   turn_type TEXT,\n"
			"	FOREIGN KEY (struct_id, residue_begin)\n"
			"		REFERENCES residues (struct_id, resNum)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(struct_id, residue_begin));";
    } else if(db_mode == "mysql")
	{
		return
			"CREATE TABLE IF NOT EXISTS beta_turns (\n"
			"	struct_id BINARY(36),\n"
			"	residue_begin INTEGER,\n"
            "   turn_type VARCHAR(255),\n"
			"	FOREIGN KEY (struct_id) REFERENCES structures (struct_id),\n"
			"	PRIMARY KEY(struct_id, residue_begin));";
	} else {
		return "";
	}

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
	boost::uuids::uuid struct_id,
	sessionOP db_session
){
	string beta_turns_stmt_string = "INSERT INTO beta_turns VALUES (?,?,?);";
	statement beta_turns_stmt(
		safely_prepare_statement(beta_turns_stmt_string, db_session));
    
	for(SSize begin=1; begin <= SSize( pose.total_residue() - beta_turn_length ); ++begin){
        Size end = begin + beta_turn_length;
        
        if ( !residue_range_is_relevant( relevant_residues, begin, end ) || !residue_range_is_protein( pose, begin, end ) || !all_turn_residues_are_on_the_same_chain( pose, begin ) || !beta_turn_present( pose, begin ) )
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

void BetaTurnDetectionFeatures::setup_conformation_to_turn_type_map()
{
    if ( initialized_ ) return;
    initialized_ = true;
    
	// These turn types are well characterized
	conformation_to_turn_type_[ "AA" ] = "I";	
	conformation_to_turn_type_[ "LL" ] = "I'";
	conformation_to_turn_type_[ "BL" ] = "II";
	conformation_to_turn_type_[ "EA" ] = "II'";
	conformation_to_turn_type_[ "AB" ] = "VIII";
	
	// Type IV is essentially "other"
	conformation_to_turn_type_[ "XX" ] = "IV";
	
	// These types are the working names for some new turn types
	conformation_to_turn_type_[ "LE" ] = "VIII'";
	conformation_to_turn_type_[ "AL" ] = "IX";
	conformation_to_turn_type_[ "LA" ] = "IX'";
	
	// These turns have a CisProline at position three
	/* My current binning of Ramachandran space is too coarse to differentiate these types, so I will just make a Type VI turn type
	conformation_to_turn_type_[ "BACis" ] = "VIa1";
	conformation_to_turn_type_[ "BACis" ] = "VIa2";
	conformation_to_turn_type_[ "BACis" ] = "VIb";
	*/
	conformation_to_turn_type_[ "BACis" ] = "VI";
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

bool BetaTurnDetectionFeatures::residue_range_is_relevant( vector1< bool > const & relevant_residues, Size range_begin, Size range_end ) const
{
    for ( Size current_residue = range_begin; current_residue <= range_end; ++current_residue )
    {
        if ( current_residue > relevant_residues.size() || !relevant_residues[ current_residue ] )
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
	string rama_hash = "";
	for ( Size residue_number = first_residue + 1; residue_number < first_residue + beta_turn_length; ++residue_number )
	{
		rama_hash += determine_ramachandran_hash( pose.phi( residue_number ), pose.psi( residue_number ), pose.omega( residue_number ) );
	}
	
	if ( rama_hash.size() > 2 )
	{
		if ( pose.residue( first_residue + beta_turn_length - 1 ).name().compare( "PRO" ) )
		{
			rama_hash = "XX";
		}
	}
	
	if (!conformation_to_turn_type_.count( rama_hash ) )
	{
		rama_hash = "XX";
	}
	
	return conformation_to_turn_type_[ rama_hash ];
}

string BetaTurnDetectionFeatures::determine_ramachandran_hash( Real phi, Real psi, Real omega ) const
{
	
	// This method uses extremely crude partitioning of Ramachandran space.  
	// THESE DEFINITIONS ARE TEMPORARY AND WILL BE UPDATED SOON.
	string peptide_bond_isomerization;
	if ( omega > -90 && omega <= 90 )
	{
		peptide_bond_isomerization = "Cis";
	}
	
	
    if ( phi > 0. && phi <= 180. )
	{
        if ( psi > -50. && psi <= 100. )
		{
            return "L" + peptide_bond_isomerization;
		}
		else
		{
			return "E" + peptide_bond_isomerization;
		}
	}
	else 
	{
		if ( psi > -100. && psi <= 50. )
		{
			return "A" + peptide_bond_isomerization;
		}
		else
		{
			return "B" + peptide_bond_isomerization;
		}
	}
	return "X" + peptide_bond_isomerization;
}

} // namesapce features
} // namespace protocols
