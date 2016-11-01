// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_hbonds_hbonds_geom_hh
#define INCLUDED_core_scoring_hbonds_hbonds_geom_hh

// Package headers
#include <core/scoring/DerivVectorPair.fwd.hh>

#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/hbtrie/HBAtom.fwd.hh>

// Project headers
#include <core/chemical/types.hh>
#include <core/conformation/Residue.fwd.hh>


#include <numeric/constants.hh>


namespace core {
namespace scoring {
namespace hbonds {

void
fade_energy(
	Real & energy,
	Real & dE_dr = DUMMY_DERIV,
	Real & dE_dxD = DUMMY_DERIV,
	Real & dE_dxH = DUMMY_DERIV,
	Real & dE_dBAH = DUMMY_DERIV,
	Real & dE_dchi = DUMMY_DERIV
);


HBAccChemType
get_hb_acc_chem_type(
	Size const aatm,
	conformation::Residue const & acc_rsd
);

HBDonChemType
get_hb_don_chem_type(
	Size const datm,
	conformation::Residue const & don_rsd
);

HBEvalTuple
hbond_evaluation_type(
	hbtrie::HBAtom const & datm,
	Size const don_rsd,
	hbtrie::HBAtom const & aatm,
	Size const acc_rsd
);

HBSeqSep
get_seq_sep(
	HBDonChemType const & don_chem_type,
	HBAccChemType const & acc_chem_type,
	int const & sep
);

// hbond evaluation type -- determines what scoring function to use
extern Real DUMMY_DERIV;
extern bool DUMMY_BOOL;
extern HBGeoDimType DUMMY_HBGEODIMTYPE;
extern HBondDerivs DUMMY_DERIVS;
extern HBondDerivs const ZERO_DERIV2D;


HBEvalTuple
hbond_evaluation_type(
	Size const datm,
	conformation::Residue const & don_rsd,
	Size const aatm,
	conformation::Residue const & acc_rsd
);

// @brief Inner-most function for hydrogen bond energy evaluation.  This function
// evaluates the polynomials describing the energy for a hydrogen bond in terms of
// its geometry.  AHdis is the distance between the acceptor and the hydrogen, xD is
// the cosine of (180 - the angle between the acceptor, the hydrogen and the donor heavy atom),
// xH is the cosine of (180 - the angle between the acceptor-base, the acceptor, and the hydrogen),
// and chi is the dihedral defined by the acceptor-base 2, the acceptor base, the acceptor, and
// the hydrogen.  The chi term only contributes to the energy if hbondoptions.sp2_chi_penalty()
// returns true.  If the last 6 parameters are not specified, then the derivatives are not evaluated.
void
hbond_compute_energy(
	HBondDatabase const & database,
	HBondOptions const & hbondoptions,
	HBEvalTuple const & hbt,   // used internally & by geometric solvation
	Real const AHdis, // acceptor proton distance
	Real const xD,    // -cos(180-theta), where theta is defined by Tanja K.
	Real const xH,      // cos(180-phi), where phi is defined by Tanja K.
	Real const chi,     // AB2-AB-A-H dihdral angle for sp2 hybridized acceptors
	Real & energy,      // main return value #1: sum of the dAH term, the xD term and the xH term.
	bool & apply_chi_torsion_penalty = DUMMY_BOOL, // did this hbond get the chi torsion penalty?
	HBGeoDimType & AHD_geometric_dimension = DUMMY_HBGEODIMTYPE, // measure in angle in cosine space?
	Real & dE_dr = DUMMY_DERIV,
	Real & dE_dxD = DUMMY_DERIV,
	Real & dE_dxH = DUMMY_DERIV,
	Real & dchipen_dBAH = DUMMY_DERIV,
	Real & dchipen_dchi = DUMMY_DERIV
);

/// @brief Evaluate the hydrogen bond energy and derivatives after having first calculated
/// the HD and BA *u*nit vectors
void
hb_energy_deriv_u(
	HBondDatabase const & database,
	HBondOptions const & hbondoptions,
	HBEvalTuple const & hbt, // hbond evaluation tuple
	Vector const & Hxyz, // proton
	Vector const & Dxyz, // donor -- only needed for derivative evaluation
	Vector const & HDunit, // proton-to-donor unit vector
	Vector const & Axyz, // acceptor
	Vector const & Bxyz, // pseudo acceptor-base coordinate -- only needed for derivative evaluation
	Vector const & BAunit, // unit vector towards the acceptor base
	Vector const & B2xyz, // coordinate of acceptor-base 2
	Real & energy,
	bool const calculate_derivative = false,
	HBondDerivs & deriv = DUMMY_DERIVS
);

/// @brief Evaluate the hydrogen bond energy and derivatives after having first calculated
/// the HD and BA *u*nit vectors; deriv type must have been chosen (why does this exist?)
void
hb_energy_deriv_u2(
	HBondDatabase const & database,
	HBondOptions const & hbondoptions,
	HBEvalTuple const & hbt, // hbond evaluation type
	HBDerivType const deriv_type,
	Vector const & Hxyz, // proton
	Vector const & Dxyz, // donor
	Vector const & HDunit, // donor
	Vector const & Axyz, // acceptor
	Vector const & Bxyz, // pseudo acceptor-base coordinate
	Vector const & BAunit, // unit vector towards base
	Vector const & B2xyz, // coordinate of acceptor-base 2
	Real & energy,
	HBondDerivs & deriv = DUMMY_DERIVS
);


//hbond evaluation type -- determines what scoring function to use
void
hb_energy_deriv(
	HBondDatabase const & database,
	HBondOptions const & hbondoptions,
	HBEvalTuple const & hbt, // hbond evaluation type
	Vector const & Dxyz, // donor coords
	Vector const & Hxyz, // proton
	Vector const & Axyz, // acceptor
	Vector const & Bxyz, // acceptor base
	Vector const & B2xyz, // 2nd acceptor base for ring & SP2 acceptors
	Real & energy,
	bool const calculate_derivative = false,
	HBondDerivs & deriv = DUMMY_DERIVS // f1/f2 for four atoms
);

void
hb_energy_deriv(
	HBondDatabase const & database,
	HBondOptions const & hbondoptions,
	HBEvalTuple const & hbt, // hbond evaluation type
	Vector const & Dxyz, // donor coords
	Vector const & Hxyz, // proton
	Vector const & Axyz, // acceptor
	Vector const & Bxyz, // acceptor base
	Vector const & B2xyz, // 2nd acceptor base for ring & SP2 acceptors
	Real & energy,
	HBDerivType const deriv_type,
	HBondDerivs & deriv = DUMMY_DERIVS // f1/f2 for four atoms
);

Vector
create_acc_orientation_vector(
	HBondOptions const & hbondoptions,
	conformation::Residue const & residue,
	int atom_id
);

// hack?
extern Real DUMMY_DERIV;


Vector
create_don_orientation_vector(
	conformation::Residue const & residue,
	int atom_id
);

void
make_hbBasetoAcc_unitvector(
	HBondOptions const & hbondoptions,
	chemical::Hybridization const & acc_hybrid,
	Vector const & Axyz,
	Vector const & Bxyz,
	Vector const & B2xyz,
	Vector & PBxyz, /// the coordinate for the pseudo-acceptor-base atom, used in derivative evaluation
	Vector & BAunit
);

void
assign_abase_derivs(
	HBondOptions const & hbondoptions,
	conformation::Residue const & acc_rsd,
	Size acc_atom,
	chemical::Hybridization const & acc_hybrid,
	DerivVectorPair const & abase_deriv,
	Real weighted_energy,
	utility::vector1< DerivVectorPair > & acc_atom_derivs
);

void
assign_abase_derivs(
	HBondOptions const & hbondoptions,
	conformation::Residue const & acc_rsd,
	Size acc_atom,
	HBEvalTuple const & hbe_type,
	DerivVectorPair const & abase_deriv,
	Real weighted_energy,
	utility::vector1< DerivVectorPair > & acc_atom_derivs
);


}
}
}

#endif
