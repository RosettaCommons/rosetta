// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/geometric_solvation/DatabaseOccSolEne.hh
/// @brief  Database containing params for OccludedHbondSolEnergy (i.e., for the pwSHO model)
/// @author John Karanicolas
/// @author Andrea Bazzoli


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
			debug_assert ( donor_occ_data_[ polar_atom_type_index][occ_atom_type_index][des_param] > NO_PARAMETERS_ + 0.01 );
			return donor_occ_data_[ polar_atom_type_index][occ_atom_type_index][des_param];
		}
		debug_assert ( acc_occ_data_[ polar_atom_type_index][occ_atom_type_index][des_param] > NO_PARAMETERS_ + 0.01 );
		return acc_occ_data_[ polar_atom_type_index][occ_atom_type_index][des_param];
	}

	Real const & atomic_interaction_cutoff() const { return atomic_interaction_cutoff_; }

	/// @brief tells, for a given atom type, if there are pwSHO parameters for donors having that atom type
	bool don_type_has_data( Size atom_type_index ) const { return don_type_has_data_[ atom_type_index ]; }

	/// @brief tells, for a given atom type, if there are pwSHO parameters for acceptors having that atom type
	bool acc_type_has_data( Size atom_type_index ) const { return acc_type_has_data_[ atom_type_index ]; }

	/// @brief tells, for a given atom type, if there are pwSHO parameters for donor-occluding atoms having
	/// that type
	bool don_occ_type_has_data( Size atom_type_index ) const { return don_occ_type_has_data_[ atom_type_index ]; }

	/// @brief tells, for a given atom type, if there are pwSHO parameters for acceptor-occluding atoms having
	/// that type
	bool acc_occ_type_has_data( Size atom_type_index ) const { return acc_occ_type_has_data_[ atom_type_index ]; }

	/// @brief maps each atom type to one for which there are pwSHO parameters when atoms of the new type are
	/// evaluated as donors
	Size don_type_mapping( Size atom_type_index ) const { return don_type_mapping_[ atom_type_index ]; }

	/// @brief maps each atom type to one for which there are pwSHO parameters when atoms of the new type are
	/// evaluated as acceptors
	Size acc_type_mapping( Size atom_type_index ) const { return acc_type_mapping_[ atom_type_index ]; }

	/// @brief maps each atom type to one for which there are pwSHO parameters when atoms of the new type are
	/// evaluated as occluding solvent H-bonds to donors
	Size don_occ_type_mapping( Size atom_type_index ) const { return don_occ_type_mapping_[ atom_type_index ]; }

	/// @brief maps each atom type to one for which there are pwSHO parameters when atoms of the new type are
	/// evaluated as occluding solvent H-bonds to acceptors
	Size acc_occ_type_mapping( Size atom_type_index ) const { return acc_occ_type_mapping_[ atom_type_index ]; }


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

	// tell, for each atom type, if there are pwSHO parameters for donors, acceptors, atoms occluding solvent H-bonds to
	// donors, and atoms occluding solvent H-bonds to acceptors, respectively, having that atom type.
	utility::vector1< bool > don_type_has_data_;
	utility::vector1< bool > acc_type_has_data_;
	utility::vector1< bool > don_occ_type_has_data_;
	utility::vector1< bool > acc_occ_type_has_data_;

	// map each atom type to one for which the pwSHO has parameters when atoms of the new type are evaluated as donors,
	// acceptors, atoms occluding solvent H-bonds to donors, and atoms occluding solvent H-bonds to acceptors, respectively.
	utility::vector1< Size > don_type_mapping_;
	utility::vector1< Size > acc_type_mapping_;
	utility::vector1< Size > don_occ_type_mapping_;
	utility::vector1< Size > acc_occ_type_mapping_;

	/// @brief maps each atom type to itself
	void init_to_self_mapping( utility::vector1< Size > & v );

	/// @brief maps donor types not having pwSHO parameters to donor types having them
	void init_don_mapping(chemical::AtomTypeSet const & atom_set);

	/// @brief maps acceptor types not having pwSHO parameters to acceptor types having them
	void init_acc_mapping(chemical::AtomTypeSet const & atom_set);

	/// @brief maps occluding-atom types not having pwSHO parameters to occluding-atom types having them
	void init_atom_occ_mapping(
		chemical::AtomTypeSet const & atom_set,
		utility::vector1< Size > & atom_occ_type_mapping,
		utility::vector1< bool > const & atom_occ_type_has_data);

	// energies below this can be neglected, used in determining jumpout geometry conditions
	Real const min_occ_energy_;

	// max distance between interacting atoms (computed from min_occ_energy)
	Real atomic_interaction_cutoff_;

	// dummy value indicating the absence of pwSHO parameters for a given atom type
	Real const NO_PARAMETERS_;
};

} // geometric_solvation
} // scoring
} // core

#endif // INCLUDED_core_scoring_geometric_solvation_DatabaseOccSolEne_HH

