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

HBAccChemType
get_hb_acc_chem_type(
	int const aatm,
	conformation::Residue const & acc_rsd
);

HBDonChemType
get_hb_don_chem_type(
	int const datm,
	conformation::Residue const & don_rsd
);

HBEvalTuple
hbond_evaluation_type(
	hbtrie::HBAtom const & datm,
	int const & don_rsd,
	hbtrie::HBAtom const & aatm,
	int const & acc_rsd
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
	int const datm,
	conformation::Residue const & don_rsd,
	int const aatm,
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
	HBEvalTuple hbt,   // used internally & by geometric solvation
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

/// @brief Fade the energy smoothly to zero over the energy range [-0.1, 0.1]
/// @detail Because of the additive functional form, in order to make
///derivative continuous at the boundary of definition, we fade the
///energy function smoothly to zero.
///
/// Check that f(x) = -0.025 + 0.5x - 2.5x^2 satisfies
///     f(-.1) = -0.025 +   0.5*(-.1) - 2.5*(-.1)^2 = -.1
///     f( .1) = -0.025 +    0.5*(.1) -  2.5*(.1)^2 = 0
///    f'(-.1) =  0.5   - 2.5*2*(-.1)               = 1
///     f'(.1) =  0.5   -  2.5*2*(.1)               = 0
inline
void
fade_energy(
	Real & energy,
	Real & dE_dr = DUMMY_DERIV,
	Real & dE_dxD = DUMMY_DERIV,
	Real & dE_dxH = DUMMY_DERIV,
	Real & dE_dBAH = DUMMY_DERIV,
	Real & dE_dchi = DUMMY_DERIV
) {
	Real const input_energy( energy );
	if ( input_energy > 0.1L ) {
		energy = 0;
		if ( &dE_dxH != &DUMMY_DERIV ) {
			dE_dr  = 0;
			dE_dxD = 0;
			dE_dxH = 0;
			dE_dBAH = 0;
			dE_dchi = 0;
		}
	} else if ( input_energy > -0.1L ) {
		energy = -0.025 + 0.5*energy - 2.5*energy*energy;
		if ( &dE_dxH != &DUMMY_DERIV ) {
			dE_dr  *= 5*(0.1-input_energy);
			dE_dxD *= 5*(0.1-input_energy);
			dE_dxH *= 5*(0.1-input_energy);
			dE_dBAH *= 5*(0.1-input_energy);
			dE_dchi *= 5*(0.1-input_energy);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
///
/// @brief
///Evaluate the Base-Acceptor-Hydrogen angle and Base-Acceptor
///torsion portion of the hydrogen bond energy and derivatives
///
/// @details
/// Formula #11
///
/// F (chi=0 or chi=pi)               | G (chi=pi/2 or chi=3*pi/2)  |
///------\                   /--------|-----\_              _/------|-  m - 0.5
///|      \                 /         |       \_          _/        |-  1
///m       \               /          |         \_      _/          |
///|        \     /-\     /        ---|           \----/            |-  d - 0.5
///|_    |   \___/   \___/         _d_|                             |_  -0.5
///      |<-l->|                      |                             |
///      |     |<-BAH=2pi/3
///      |
///      |<-BAh=2pi/3 - l
////
///
///BAH := Base-Acceptor-Hydrogen interior Angle
///       BAH=pi when linear and BAH=pi/2 when perpendicular
///
///chi :  Torsion angle defined by ABase2-Base-Acceptor-Hydrogen
///       The Sp2 orbials are in the ABase2-Base-Acceptor plane
///       For Backbone acceptors ABase2=C-alpha
///
///  d := distance from minimum value of -0.5 at BAH=120 to BAH=180 in f
///       defined by HBondOptions::sp2_BAH180_rise() which is set by
///       -corrections:score:hb_sp2_BAH180_rise flag and
///       defaults to 0.75
///
///  m := distance from minimum to maximum values of f
///       must rise high enough so that
///           E_fade_max = maxBAH_CHI + minAHD + minAHdist
///                (0.1) = m + (minBAH_CHI) + (-0.5) + (-0.5)
///                   m  = 1.6
///  l := period/2 of the BAH=120 to BAH=60 piece of F
///       emperically fit to be 0.357
///
///  F := d/2 * cos(3(pi-BAH) + d/2 - 0.5                        BAH > 2pi/3
///       m/2 * cos(pi - (2pi/3 - BAH)/l) + m/2 - 0.5    2pi/3 > BAH > pi(2/3 - l)
///       m-0.5                                    pi(2/3 - l) > BAH
///
///  G := d - 0.5                                                BAH > 2pi/3
///       (m-d)/2 * cos(pi - (2pi/3 - BAH)/l) + (m-d)/2 + d + 0.5
///                                                      2pi/3 > BAH > pi(2/3 - l)
///       m-0.5                                    pi(2/3 - l) > BAH
///
///
///  H := inteprolate smoothly betwen F and G going around chi
///       (cos(2*chi) + 1)/2
///
///  E := Energy for BAH/CHI term
///       H*F + (1-H)*G
///
/// dE/dchi := dH/dchi*f - dH/dchi*g
///          = -sin(2*chi)*F + sin(2*chi)*G
///
/// dE/dBAH := H*dF/dBAH + (1-H)*dG/dBAH
///
/// dF/dBAH := 3 * d/2 * sin(3(pi-BAH))                           BAH > 2pi/3
///            m/2 * -1/l * sin(pi - (2pi/3 - BAH)/l)     2pi/3 > BAH > pi(2/3 - l)
///            0                                    pi(2/3 - l) > BAH
///
/// dG/dBAH := 0                                                  BAH > 2pi/3
///            (m-d)/2 * -1/l * sin(pi - (2pi/3 - BAH)/)  2pi/3 > BAH > pi(2/3 - l)
///            0                                    pi(2/3 - l) > BAH
inline
void
bah_chi_compute_energy_sp2(
	Real const d,
	Real const m,
	Real const l,
	Real const xH,
	Real const chi,
	Real const acc_don_scale,
	Real & energy,
	Real & dE_dBAH,
	Real & dE_dchi
) {
	using std::cos;
	using std::sin;
	using numeric::constants::d::pi;

	Real const PI_minus_BAH( acos(xH) );
	Real const BAH = pi - ( PI_minus_BAH );

	Real const  H((cos(2*chi) + 1) * 0.5);
	Real F(0), G(0);

	if ( BAH >= pi * 2/3 ) {
		F = d/2 * cos(3 * PI_minus_BAH) + d/2 - 0.5;
		G = d - 0.5;
	} else if ( BAH >= pi * (2/3 - l) ) {
		Real const outer_rise(cos(pi - (pi*2/3 -  BAH)/l));
		F = m/2 * outer_rise + m/2 - 0.5;
		G = (m - d)/2 * outer_rise + (m - d)/2 + d - 0.5;
	} else {
		F = m-0.5;
		G = m-0.5;
	}

	energy += acc_don_scale * ( H*F + (1-H)*G );

	if ( &dE_dchi != &DUMMY_DERIV ) {
		Real const dH_dchi(-1 * sin(2*chi));
		Real dF_dBAH(0), dG_dBAH(0);
		if ( BAH >= pi * 2/3 ) {
			dF_dBAH = 3 * d/2 * sin(3 * PI_minus_BAH);
		} else if ( BAH >= pi * (2/3 - l) ) {
			Real const d_outer_rise_dBAH( -1/l * sin(pi - (2*pi/3 - BAH)/l) );
			dF_dBAH = m/2 * d_outer_rise_dBAH;
			dG_dBAH = (m - d)/2 * d_outer_rise_dBAH;
		}
		dE_dchi = acc_don_scale * ( F*dH_dchi - G*dH_dchi );
		dE_dBAH = acc_don_scale * ( H*dF_dBAH + (1-H)*dG_dBAH );
	}
}


inline
void
bah_chi_compute_energy_sp3(
	Real const /*xH*/,
	Real const chi,
	Real const acc_don_scale,
	Real & energy,
	Real & dE_dBAH,
	Real & dE_dchi
) {
	// just add in a penalty directly to the energy sum; the chi-penalty
	// is only multiplied in for the sp2 term.
	Real const max_penalty = 0.25;
	Real cos2ChiShifted = max_penalty * ( 1 + std::cos(chi)) / 2;
	energy += acc_don_scale * cos2ChiShifted;

	if ( &dE_dBAH != &DUMMY_DERIV ) {
		dE_dchi = -1 * max_penalty * std::sin(chi)/2 * acc_don_scale;
	}
}


/// @brief Evaluate the hydrogen bond energy and derivatives after having first calculated
/// the HD and BA *u*nit vectors
void
hb_energy_deriv_u(
	HBondDatabase const & database,
	HBondOptions const & hbondoptions,
	HBEvalTuple const hbt, // hbond evaluation tuple
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
	HBEvalTuple const hbt, // hbond evaluation type
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
	HBEvalTuple const hbt, // hbond evaluation type
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
	HBEvalTuple const hbt, // hbond evaluation type
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
	HBEvalTuple const hbe_type,
	DerivVectorPair const & abase_deriv,
	Real weighted_energy,
	utility::vector1< DerivVectorPair > & acc_atom_derivs
);


}
}
}

#endif
