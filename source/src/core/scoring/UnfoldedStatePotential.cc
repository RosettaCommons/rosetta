// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/UnfoldedStatePotential.cc
/// @brief  Unfolded state energies based on energies of residues in fragments, definition file
/// @author Ron Jacak (ronj@email.unc.edu)
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <core/scoring/UnfoldedStatePotential.hh>

// Project headers

#include <core/conformation/Residue.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/pose/Pose.hh>

// Numeric headers

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

#include <utility/vector1.hh>

static thread_local basic::Tracer tr( "core.scoring.UnfoldedStatePotential" );


namespace core {
namespace scoring {


UnfoldedStatePotential::UnfoldedStatePotential( std::string const & filename ) {
	read_database_file( filename );
}

UnfoldedStatePotential::~UnfoldedStatePotential() {}


void
UnfoldedStatePotential::read_database_file( std::string const & filename ) {

	using namespace core::chemical;
	using namespace core::scoring;

	if( ! utility::file::file_exists( filename ) ) {
		utility_exit_with_message("Cannot find file '"+filename+"'");
	}

	utility::io::izstream data( filename );
	if ( !data.good() ) {
		utility_exit_with_message("Cannot open file '"+filename+"'");
	}

	// read in all lines in file
	utility::vector1< std::string > lines;
	std::string line;
	while ( getline( data, line ) ) {
		std::istringstream l( line );
		if ( line.size() < 1 || line[0] == '#' ) continue; // skip comment lines
		lines.push_back( line );
	}
	data.close();

	// parse the first line that contains the score types
	std::string const & first_line( lines[1] );
	std::istringstream h( first_line );
	std::string tag;
	utility::vector1< ScoreType > unfolded_score_types;

	h >> tag;
	if ( tag != "AA" ) {
		utility_exit_with_message("Error parsing first line of '"+filename+"'");
	}
	while ( !h.fail() ) {
		h >> tag;
		if ( h.fail() ) break;
		if ( ScoreTypeManager::is_score_type( tag ) )  {
			unfolded_score_types.push_back( ScoreTypeManager::score_type_from_name( tag ) );
		} else {
			utility_exit_with_message("Error, score type '"+tag+"' spcified in '"+filename+"' doesn't exist.");
		}
	}

	// parse the second line that contains the weights
	std::string const & second_line( lines[2] );
	std::istringstream k( second_line );
	Size const ntypes( unfolded_score_types.size() );
	Real weight;

	k >> tag;
	if ( tag != "WEIGHT" ) {
		utility_exit_with_message("Error parsing second line of '"+filename+"'");
	}
	for ( Size i=1; i <= ntypes; ++i ) {
		k >> weight;
		if ( k.fail() ) {
			utility_exit_with_message("Error, number of energies doesn't match number of score types in '"+filename+"'");
		}
		unfolded_potential_file_weights_[ unfolded_score_types[i] ] = weight; // add energy to appropriate score type in temp emap
	}

	// parse the rest of the file
	Size const nlines( lines.size() );

	for ( Size i=3; i <= nlines; ++i ) {
		std::string const & temp_line( lines[i] );
		std::istringstream l( temp_line );
		std::string tlc;
		Real val;
		EnergyMap emap;

		l >> tlc;
		ObjexxFCL::lpad( tlc, 3 ); // RNA codes are rA, etc. (two letters)

		for ( Size i=1; i <= ntypes; ++i ) {
			l >> val;
			if ( l.fail() ) {
				utility_exit_with_message("Error, number of energies doesn't match number of score types in '"+filename+"'");
			}
			emap[ unfolded_score_types[i] ] = val; // add energy to appropriate score type in temp emap
		}

		// add tlc/emap pair to map
		unfolded_energy_[ tlc ] = emap;
	}

}

scoring::EnergyMap
UnfoldedStatePotential::get_unfoled_potential_file_weights() const {
	return unfolded_potential_file_weights_;
}


void
UnfoldedStatePotential::raw_unfolded_state_energymap( std::string const & aa_name3, scoring::EnergyMap & e ) const {

	// the energies stored in the database file aren't probabilities; don't take the log of them, that doesn't make sense!
	e = ( unfolded_energy_.find( aa_name3 ) )->second;

}


void
UnfoldedStatePotential::pose_raw_unfolded_state_energymap( pose::Pose const & pose, scoring::EnergyMap & e ) const {

	EnergyMap unf_total;

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( ! pose.residue(ii).is_protein() && ! pose.residue(ii).is_RNA()  ) {
			continue;
		}

		// EnergyMap's know how to add, so we can take advantage of that feature to come up with a total unfolded state
		// energy (broken down by score type anyway).
		unf_total += ( unfolded_energy_.find( pose.residue(ii).name3() ) )->second;
	}

	e = unf_total;
}

} // namespace scoring
} // namespace core

