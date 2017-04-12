// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/hbonds/hbonds.hh
///
/// @brief Set of utilities for detection of hydrogen bonds
/// @author Lots of people...

#ifndef INCLUDED_core_scoring_hbonds_hbonds_hh
#define INCLUDED_core_scoring_hbonds_hbonds_hh

// Package Headers
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#include <core/pose/Pose.fwd.hh>

// Unit Headers
#include <utility/vector1.hh>
#include <boost/unordered_map.hpp>

namespace core {
namespace scoring {
namespace hbonds {

//fpd ss-len-dep weight parameters
struct SSWeightParameters {
	SSWeightParameters() : ssdep_(false), l_(0.5), h_(2.0), len_l_(4), len_h_(17) {}
	bool ssdep_;
	Real l_, h_;
	Size len_l_, len_h_;
};

// The idea is that for each of the atoms involved in the hbond, there needs
// to be a coorespondence to the entry in the HBondDerivs struct that holds
// the derivative vector for that atom. In particular, the logic for the abase
// and abase2 atoms is complicated, and so to encapsulate that logic, we
// need this small struct and also the HBDerivAssigner.
struct AssignmentScaleAndDerivVectID {
	AssignmentScaleAndDerivVectID() : scale_( 1.0 ), dvect_id_( which_hb_unassigned ) {}
	Real scale_;
	which_atom_in_hbond dvect_id_;
};

class HBDerivAssigner
{
public:
	HBDerivAssigner(
		HBondOptions const & hbondoptions,
		HBEvalTuple hb_eval,
		conformation::Residue const & don_rsd,
		Size don_h_atm,
		conformation::Residue const & acc_rsd,
		Size acc_atm
	);

	Size h_ind() const;
	Size d_ind() const;
	Size a_ind() const;
	Size abase_ind() const;
	Size abase_prime_ind() const; // returns zero if no appropriate abase'
	Size abase2_ind() const;      // returns zero if no appropriate abase2

	Size ind( which_atom_in_hbond which );

	/// @brief For a given role in the hbond, return the scale factor and the entry
	/// in the HBondDerivs struct so that the derivative can be properly assigned
	/// to that atom.
	AssignmentScaleAndDerivVectID
	assignment( which_atom_in_hbond which );

private:
	chemical::Hybridization acc_hybrid_;
	Size h_ind_;
	Size d_ind_;
	Size a_ind_;
	Size abase_ind_;
	Size abase_prime_ind_;
	Size abase2_ind_;
	bool measure_sp3acc_BAH_from_hvy_;
};

bool
calculate_intra_res_hbonds( conformation::Residue const & rsd,
	HBondOptions const & options );

void
fill_intra_res_hbond_set(
	pose::Pose const & pose,
	bool const calculate_derivative,
	HBondSet & hbond_set,
	bool const exclude_bb  = false,
	bool const exclude_bsc = false,
	bool const exclude_scb = false,
	bool const exclude_sc  = false
);

void
fill_hbond_set(
	pose::Pose const & pose,
	bool const calculate_derivative,
	HBondSet & hbond_set,
	bool const exclude_bb  = false,
	bool const exclude_bsc = false,
	bool const exclude_scb = false,
	bool const exclude_sc  = false
);

void
fill_hbond_set(
	pose::Pose const & pose,
	bool const calculate_derivative,
	HBondSet & hbond_set,
	SSWeightParameters const &,
	bool const exclude_bb  = false,
	bool const exclude_bsc = false,
	bool const exclude_scb = false,
	bool const exclude_sc  = false
);


core::Real
get_ssdep_weight(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	SSWeightParameters const & ssdep
);

/// @brief Fill HBondSet using the distance between the acceptor and
///hydrogen atoms as the definitional cutoff. Do not exclude any
///contacts and do not evaluate derivatives.
void
fill_hbond_set_by_AHdist_threshold(
	pose::Pose const & pose,
	Real const AHdist_threshold,
	HBondSet & hbond_set);

/// @brief Increment the appropriate places in the input energy map for a particular hydrogen bond
void
increment_hbond_energy(
	HBEvalType const & hbe_type,
	EnergyMap & emap,
	Real hbE
);

/// @brief Increment the appropriate places in the input energy map for a particular hydrogen bond
void
increment_npd_hbond_energy(
	HBEvalType const & hbe_type,
	EnergyMap & emap,
	Real hbE,
	bool intra_res
);

void
get_hbond_energies(
	HBondSet const & hbond_set,
	EnergyMap & emap
);

/*void
get_hbond_energies(
HBondSet const & hbond_set,
EnergyMap & emap);*/

void
identify_hbonds_1way(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_don_bb,
	bool const exclude_don_bsc,
	bool const exclude_acc_scb,
	bool const exclude_acc_sc,
	// output
	HBondSet & hbond_set,
	Real ssdep_weight_factor = 1.0
);

void
identify_hbonds_1way(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_don_bb,
	bool const exclude_don_bsc,
	bool const exclude_acc_scb,
	bool const exclude_acc_sc,
	HBondOptions const & options,
	// output
	EnergyMap & emap,
	Real ssdep_weight_factor = 1.0
);

void
identify_hbonds_1way(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_don_bb,
	bool const exclude_don_bsc,
	bool const exclude_acc_scb,
	bool const exclude_acc_sc,
	HBondOptions const & options,
	// output
	EnergyMap & emap,
	boost::unordered_map<core::Size, core::Size> & num_hbonds,
	Real ssdep_weight_factor = 1.0
);

/// @brief Returns the energy for the hydrogen bond between a given don/acceptor
/// pair
Real
hb_energy(
	HBondDatabase const & database,
	HBondOptions const & hbondoptions,
	HBondSet const & hbset,
	conformation::Residue const & acc,
	Size acc_atm,
	conformation::Residue const & don,
	Size don_atm
);

void
identify_hbonds_1way_AHdist(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	Real const AHdist_threshold,
	HBondSet & hbond_set
);

void
identify_intra_res_hbonds(
	HBondDatabase const & database,
	conformation::Residue const & rsd,
	Size const rsd_nb,
	bool const evaluate_derivative,
	HBondSet & hbond_set,
	bool const exclude_bb  = false,
	bool const exclude_bsc = false,
	bool const exclude_scb = false,
	bool const exclude_sc  = false
);


void
identify_intra_res_hbonds(
	HBondDatabase const & database,
	conformation::Residue const & rsd,
	Size const rsd_nb,
	HBondOptions const & options,
	EnergyMap & emap);

Real
get_environment_dependent_weight(
	HBEvalTuple const & hbe_type,
	int const don_nb,
	int const acc_nb,
	HBondOptions const & options);

Real
hb_eval_type_weight(
	HBEvalType const & hbe_type,
	EnergyMap const & emap,
	bool const intra_res = false,
	bool const put_intra_into_total = true
);

Real
npd_hb_eval_type_weight(
	HBEvalType const & hbe_type,
	EnergyMap const & emap,
	bool const intra_res = false,
	bool const put_intra_into_total = true
);


bool
nonzero_hbond_weight( ScoreFunction const & scorefxn );

//Membranes
// TODO: clean up this code duplication (will not make sense with the new framework):
// Possible Clean up approach:
// (1) Create base class HBondEnvironment
// (2) Derive -> HBondNeighborCountEnvironment which has the number of neighbors for the donor and acceptor residues
// (3) Derive -> HBondMembraneEnvironment which stores the MembraneEmbed or some subset of it
// (4) extend the interface to hbond_set to get a HBondEnvironment.  This will require adding additional member data to the hbond_set
// (5) modify the identify_hbonds_1way interface to take a reference to an HBondEnvironment
// (6) use the passed in HBondEnvironment to compute the environmental weight in identify_hbonds_1way
//  Note: care will have to be taken to not create and destroy lots of copies of HBondEnvironment objects

// Once this has been done, please remove all traces of the identify_hbonds_1way_membrane functions.

/// @brief Identify Membrane Hydrogen Bonds (Env)
/// @details Corrects for strength of hydrogen bonding in the membrane
/// (depth-dependent). This version of the method switches between
/// the previous membrane code (MembraneEmbed cached to the pose) and updated
/// RosettaMP Framework (2015).
/// @note This code still preserves duplicated hbond interface for membranes
/// TODO: REMOVE THE CODE DUPLICATION!!!!
void
identify_hbonds_1way_membrane(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_don_bb,
	bool const exclude_don_bsc,
	bool const exclude_acc_scb,
	bool const exclude_acc_sc,
	HBondSet & hbond_set,
	pose::Pose const & pose
);

void
identify_hbonds_1way_membrane(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_don_bb,
	bool const exclude_don_bsc,
	bool const exclude_acc_scb,
	bool const exclude_acc_sc,
	HBondOptions const & options,
	EnergyMap & emap,
	pose::Pose const & pose
);

Real
get_membrane_depth_dependent_weight(
	pose::Pose const & pose,
	int const don_nb,
	int const acc_nb,
	Vector const & Hxyz,
	Vector const & Axyz
);

Real
get_membrane_depth_dependent_weight(
	Vector const & normal,
	Vector const & center,
	Real const & thickness,
	Real const & steepness,
	int const don_nb,
	int const acc_nb,
	Vector const & Hxyz, // proton
	Vector const & Axyz  // acceptor
);

} // hbonds
} // scoring
} // core

#endif // INCLUDED_core_scoring_hbonds_hbonds_hh
