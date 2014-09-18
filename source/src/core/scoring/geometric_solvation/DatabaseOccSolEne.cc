// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/scoring/geometric_solvation/DatabaseOccSolEne.cc
/// @brief  Database containing params for OccludedHbondSolEnergy
/// @author John Karanicolas


// Unit Headers
#include <core/scoring/geometric_solvation/DatabaseOccSolEne.hh>

// Package headers


// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/scoring/ScoringManager.fwd.hh>

// AUTO-REMOVED #include <basic/database/open.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/io/izstream.hh>


#include <cmath>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace geometric_solvation {

static thread_local basic::Tracer tr( "core.scoring.DatabaseOccSolEne" );


/// @details ctor, reads data file. Need to configure to allow alternate tables/atom_sets
DatabaseOccSolEne::DatabaseOccSolEne( std::string const & etable_name, Real const & min_occ_energy ):
	min_occ_energy_(min_occ_energy)
{

	if ( etable_name != scoring::FA_STANDARD_DEFAULT ) {
		utility_exit_with_message( "cannot presently use OccludedHbondSolEnergy in any mode but FA_STANDARD_DEFAULT, current mode is " + etable_name );
	}

	// get the relevant atomset
	chemical::AtomTypeSet const & atom_set
		( *chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ) );

	atomic_interaction_cutoff_ = 0.;
	tr << "reading params" << std::endl;
	read_datafile( atom_set, atom_set.directory() + "/occluded_hbond_solvation_params.donors.txt", donor_occ_data_, true );
	read_datafile( atom_set, atom_set.directory() + "/occluded_hbond_solvation_params.acceptors.txt", acc_occ_data_, false );

}


void
DatabaseOccSolEne::read_datafile(
	chemical::AtomTypeSet const & atom_set,
	std::string const & database_name,
	utility::vector1< utility::vector1< utility::vector1< Real > > > & occ_data_,
	bool const process_donors // needed solely to distinguish OH's in reweighting, can go away once reweighting goes into the param files
)
{

	// initialize the values
	Size const n_atom_types( atom_set.n_atomtypes() );
  occ_data_.clear();
	occ_data_.resize( n_atom_types );
	for ( Size j=1; j <= n_atom_types; ++j ) {
		occ_data_[j].resize( n_atom_types );
		for ( Size k=1; k <= n_atom_types; ++k ) {
			// initialize to an unreasonable value, so we can assert lookups are valid
			occ_data_[j][k].resize( OccFitParam_num_params, -2. );
		}
	}

	utility::io::izstream stream;
	stream.open( database_name );
	if ( !stream.good() ) {
		utility_exit_with_message( "Unable to open " + database_name );
	}

	// skip header line
	std::string line;
	// no header for now...
	//	getline( stream, line );

	// process lines containing data
	// Note: file format is: polar_atom_name  occ_atom_name  5 params (amp, 2x mu, 2x sigma)
	while ( getline( stream, line ) ) {
		std::istringstream l(line);
		// get the atom types from the names - will fail if name is not recognized
		std::string polar_atom_type_name, occ_atom_type_name;
		l >> polar_atom_type_name;
		Size const polar_atom_type_index( atom_set.atom_type_index( polar_atom_type_name ) );
		l >> occ_atom_type_name;
		Size const occ_atom_type_index( atom_set.atom_type_index( occ_atom_type_name ) );

		// this term reweights to the exact model to account for lack of actual pairwise additivity
		Real amp_reweighting(0.387829);
		if ( process_donors ) {
			if ( polar_atom_type_name == "NH2O" ) {
				amp_reweighting *= 0.8534;
			} else if ( polar_atom_type_name == "Narg" ) {
		  	amp_reweighting *= 0.9365;
		  } else if ( polar_atom_type_name == "Nbb" ) {
		    amp_reweighting *= 0.7724;
			} else if ( polar_atom_type_name == "Nlys" ) {
				amp_reweighting *= 0.8122;
			} else if ( polar_atom_type_name == "Ntrp" ) {
				amp_reweighting *= 0.9530;
			} else if ( polar_atom_type_name == "OH" ) {
				amp_reweighting *= 0.8968;
			} else {
				utility_exit_with_message( "unrecognized donor polar atom type " + polar_atom_type_name );
			}
		} else {
			if ( polar_atom_type_name == "Nhis" ) {
				amp_reweighting *= 1.0727;
			} else if ( polar_atom_type_name == "OCbb" ) {
				amp_reweighting *= 1.1669;
			} else if ( polar_atom_type_name == "OH" ) {
				amp_reweighting *= 1.0333;
			} else if ( polar_atom_type_name == "ONH2" ) {
				amp_reweighting *= 0.9672;
			} else if ( polar_atom_type_name == "OOC" ) {
				amp_reweighting *= 0.8558;
			} else {
				utility_exit_with_message( "unrecognized donor polar atom type " + polar_atom_type_name );
			}
		}

		// read data
		Real amp, dist_mu, twice_dist_sigma_sq, cos_angle_mu, twice_cos_angle_sigma_sq;
		l >> amp; l >> dist_mu; l >> twice_dist_sigma_sq; l >> cos_angle_mu; l >> twice_cos_angle_sigma_sq;
		if ( l.fail() ) utility_exit_with_message( "bad format in " + database_name);

		// store data
		occ_data_[ polar_atom_type_index ][ occ_atom_type_index ][ OccFitParam_amp ] = amp * amp_reweighting;
		occ_data_[ polar_atom_type_index ][ occ_atom_type_index ][ OccFitParam_dist_mu ] = dist_mu;
		occ_data_[ polar_atom_type_index ][ occ_atom_type_index ][ OccFitParam_twice_dist_sigma_sq ] = twice_dist_sigma_sq;
		occ_data_[ polar_atom_type_index ][ occ_atom_type_index ][ OccFitParam_cos_angle_mu ] = cos_angle_mu;
		occ_data_[ polar_atom_type_index ][ occ_atom_type_index ][ OccFitParam_twice_cos_angle_sigma_sq ] = twice_cos_angle_sigma_sq;

		// set jumpout conditions (and atomic_interaction_cutoff_)
		// at what difference from dist_mu will the distance lead to an energy of min_occ_energy_ (assuming perfect cos_angle)?
		Real dist_away_from_mu = compute_jumpout_diff( amp, twice_dist_sigma_sq );
		Real max_dist_to_compute = dist_mu + dist_away_from_mu;
		occ_data_[ polar_atom_type_index ][ occ_atom_type_index ][ OccFitParam_max_sq_dist ] = max_dist_to_compute * max_dist_to_compute;
		if ( max_dist_to_compute > atomic_interaction_cutoff_ ) atomic_interaction_cutoff_ = max_dist_to_compute;

		// at what difference from cos_angle_mu will the cosine lead to an energy of min_occ_energy_ (assuming perfect dist)?
		Real cos_angle_away_from_mu = compute_jumpout_diff( amp, twice_cos_angle_sigma_sq );
		Real min_cos_angle_to_compute = cos_angle_mu - cos_angle_away_from_mu;
		if ( min_cos_angle_to_compute < -1. ) min_cos_angle_to_compute = -1.;
		occ_data_[ polar_atom_type_index ][ occ_atom_type_index ][ OccFitParam_min_cos_angle ] = min_cos_angle_to_compute;

	}

}


Real
DatabaseOccSolEne::compute_jumpout_diff( Real const & amp, Real const & twice_sigma_sq )
{

	// At what difference from mu will this gaussian (passed in) give a value less than min_occ_energy_ ?
	assert ( amp > min_occ_energy_ );
	assert ( twice_sigma_sq > 0. );
	return sqrt ( - twice_sigma_sq * log( min_occ_energy_ / amp ) );

}


} // geometric_solvation
} // namespace scoring
} // namespace core

