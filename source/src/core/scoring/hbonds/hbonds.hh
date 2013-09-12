// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_hbonds_hbonds_hh
#define INCLUDED_core_scoring_hbonds_hbonds_hh

#include <core/scoring/hbonds/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace hbonds {

void
fill_intra_res_hbond_set(
	pose::Pose const & pose,
	bool const calculate_derivative,
	HBondSet & hbond_set
);

void
fill_hbond_set(
	pose::Pose const & pose,
	bool const calculate_derivative,
	HBondSet & hbond_set,
	bool const exclude_bb  = false,
	bool const exclude_bsc = false,
	bool const exclude_scb = false,
	bool const exclude_sc  = false);

///@breif Fill HBondSet using the distance between the acceptor and
///hydrogen atoms as the definitional cutoff. Do not exclude any
///contacts and do not evaluate derivatives.
void
fill_hbond_set_by_AHdist_threshold(
	pose::Pose const & pose,
	Real const AHdist_threshold,
	HBondSet & hbond_set);

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
	HBondSet & hbond_set
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
	EnergyMap & emap
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
  bool const evaluate_derivative,
	HBondSet & hbond_set);


void
identify_intra_res_hbonds( 
	HBondDatabase const & database,
	conformation::Residue const & rsd, 
	bool const evaluate_derivative, 
	HBondOptions const & options,
	EnergyMap & emap);

Real
get_environment_dependent_weight(
	HBEvalTuple const & hbe_type,
	int const don_nb,
	int const acc_nb,
	HBondOptions const & options);

// TODO: clean up this code duplication:

//The membrane version of identify_hbonds_1way use the MembraneEmbed
//object cached with the pose to compute a membrane-aware
//environmental dependence hydrogen bonding potential.
// See for example http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694939/

//One way to clean this up is to do the following
//  create a base class HBondEnvironment
//  Derive -> HBondNeighborCountEnvironment which has the number of neighbors for the donor and acceptor residues
//  Derive -> HBondMembraneEnvironment which stores the MembraneEmbed or some subset of it
//  extend the interface to hbond_set to get a HBondEnvironment.  This will require adding additional member data to the hbond_set
//  modify the identify_hbonds_1way interface to take a reference to an HBondEnvironment
//  use the passed in HBondEnvironment to compute the environmental weight in identify_hbonds_1way
//  Note: care will have to be taken to not create and destroy lots of copies of HBondEnvironment objects

// Once this has been done, please remove all traces of the identify_hbonds_1way_membrane functions.
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
	// output
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
	// output
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
	Vector const & Hxyz,
	Vector const & Axyz
);

Real
hb_eval_type_weight(
	HBEvalType const & hbe_type,
	EnergyMap const & emap,
	bool const intra_res
);


bool
nonzero_hbond_weight( ScoreFunction const & scorefxn );

/*void
get_atom_hbond_derivative(
	id::AtomID const & atom,
	HBondSet const & hbond_set,
	EnergyMap const & weights,
	Vector & f1,
	Vector & f2
);*/

}
}
}

#endif
