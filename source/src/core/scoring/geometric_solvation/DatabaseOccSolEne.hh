// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/geometric_solvation/DatabaseOccSolEne.hh
/// @brief  Database containing params for OccludedHbondSolEnergy
/// @author John Karanicolas


#ifndef INCLUDED_core_scoring_geometric_solvation_DatabaseOccSolEne_hh
#define INCLUDED_core_scoring_geometric_solvation_DatabaseOccSolEne_hh

// Unit Headers
//#include <core/scoring/AtomVDW.fwd.hh>

// Package headers
#include <core/chemical/AtomTypeSet.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <string>


namespace core {
namespace scoring {
namespace geometric_solvation {

enum OccFitParam {
	OccFitParam_amp = 1,
	OccFitParam_dist_mu,
	OccFitParam_twice_dist_sigma_sq,
	OccFitParam_cos_angle_mu,
	OccFitParam_twice_cos_angle_sigma_sq,
	OccFitParam_max_sq_dist,
	OccFitParam_min_cos_angle,
	OccFitParam_num_params = OccFitParam_min_cos_angle
};


class DatabaseOccSolEne : public utility::pointer::ReferenceCount {

public:

	/// @brief ctor, reads data file
	DatabaseOccSolEne( std::string const & etable_name, Real const & min_occ_energy );

	Real const &
	operator()( bool const polar_atom_donates, Size const polar_atom_type_index, Size const occ_atom_type_index, OccFitParam des_param ) const
	{
		// note: min cos_angle can be negative
		if ( polar_atom_donates ) {
			debug_assert ( donor_occ_data_[ polar_atom_type_index][occ_atom_type_index][des_param] > -1.1 );
			return donor_occ_data_[ polar_atom_type_index][occ_atom_type_index][des_param];
		}
		debug_assert ( acc_occ_data_[ polar_atom_type_index][occ_atom_type_index][des_param] > -1.1 );
		return acc_occ_data_[ polar_atom_type_index][occ_atom_type_index][des_param];
	}

	Real const & atomic_interaction_cutoff() const { return atomic_interaction_cutoff_; }

private:

	void read_datafile( chemical::AtomTypeSet const & atom_set, std::string const & database_name,
		utility::vector1< utility::vector1< utility::vector1< Real > > > & occ_data_, bool const process_donors );

	Real compute_jumpout_diff( Real const & amp, Real const & twice_sigma_sq );

private:

	// hold all data here, indexed by 1) polar_atom_type_index, 2) occ_atom_type_index, 3) parameter type
	// jk note: we don't want to change the first index to donor/acceptor hybridization, since then we'd lose atomic charge info from the fits
	// jk note: we could switch the second index to element type / radius, it's just that at present we can only look this up via a string compare
	// jk note: at the moment we have two tables, because we're indexing on polar_atom_type_index and Ser/Thr/Tyr can be donor or acceptor
	utility::vector1< utility::vector1< utility::vector1< Real > > > donor_occ_data_;
	utility::vector1< utility::vector1< utility::vector1< Real > > > acc_occ_data_;

	// energies below this can be neglected, used in determining jumpout geometry conditions
	Real const min_occ_energy_;

	// max distance between interacting atoms (computed from min_occ_energy)
	Real atomic_interaction_cutoff_;

};

} // geometric_solvation
} // scoring
} // core

#endif // INCLUDED_core_scoring_geometric_solvation_DatabaseOccSolEne_HH

