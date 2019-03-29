// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/AromaticBackboneRestraintEnergy.cc
/// @brief  Aramid backbone restraint energy method class implementation
/// @author Andrew Watkins (amw579@stanford.edu)

// Unit Headers
#include <core/scoring/methods/AromaticBackboneRestraintEnergy.hh>
#include <core/scoring/methods/AromaticBackboneRestraintEnergyCreator.hh>

// Package Headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/DerivVectorPair.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>

// options
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Utility headers
#include <numeric/conversions.hh>
#include <utility/vector1.hh>
#include <numeric/deriv/dihedral_deriv.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>

//C++ header
#include <cstdio>

namespace core {
namespace scoring {
namespace methods {

static basic::Tracer TR( "core.scoring.methods.AromaticBackboneRestraintEnergy" );

/// @details This must return a fresh instance of the AromaticBackboneRestraintEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
AromaticBackboneRestraintEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< AromaticBackboneRestraintEnergy >();
}

ScoreTypes
AromaticBackboneRestraintEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( aromatic_restraint );
	return sts;
}

/// @brief Constructor.
///
AromaticBackboneRestraintEnergy::AromaticBackboneRestraintEnergy() :
	parent( utility::pointer::make_shared< AromaticBackboneRestraintEnergyCreator >() )//,
	//std_dev_sq_(basic::options::option[ basic::options::OptionKeys::score::ring_close_shadow_constraint ])
{
	//std_dev_sq_ = std_dev_sq_*std_dev_sq_; //Since we'll only ever use the square of this value, let's store the square of the value, calculated once, rather than recalculating it a zillion times.
}

/// @brief Copy constructor.
///
AromaticBackboneRestraintEnergy::AromaticBackboneRestraintEnergy( AromaticBackboneRestraintEnergy const & /*src*/ ):
	parent( utility::pointer::make_shared< AromaticBackboneRestraintEnergyCreator >() )
{}

/// @brief Clone -- creates a copy and returns an owning pointer to the copy.
///
EnergyMethodOP
AromaticBackboneRestraintEnergy::clone() const
{
	return utility::pointer::make_shared< AromaticBackboneRestraintEnergy >();
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

void
AromaticBackboneRestraintEnergy::residue_energy(
	conformation::Residue const &rsd,
	pose::Pose const & ,//pose,
	EnergyMap &emap
) const
{
	if ( rsd.is_virtual_residue() ) return; //Skip virtual residues.

	if ( rsd.is_ortho_aramid() ) {
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 2 ] - 0.0 ) * ( rsd.mainchain_torsions()[ 2 ] - 0.0 );
	} else if ( rsd.is_pre_methylene_ortho_aramid() ) {
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 3 ] - 0.0 ) * ( rsd.mainchain_torsions()[ 3 ] - 0.0 );
	} else if ( rsd.is_post_methylene_ortho_aramid() ) {
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 2 ] - 0.0 ) * ( rsd.mainchain_torsions()[ 2 ] - 0.0 );
	} else if ( rsd.is_pre_methylene_post_methylene_ortho_aramid() ) {
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 3 ] - 0.0 ) * ( rsd.mainchain_torsions()[ 3 ] - 0.0 );

	} else if ( rsd.is_meta_aramid() ) {
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 2 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 2 ] - 180.0 );
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 3 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 3 ] - 180.0 );
	} else if ( rsd.is_pre_methylene_meta_aramid() ) {
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 3 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 3 ] - 180.0 );
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 4 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 4 ] - 180.0 );
	} else if ( rsd.is_post_methylene_meta_aramid() ) {
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 2 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 2 ] - 180.0 );
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 3 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 3 ] - 180.0 );
	} else if ( rsd.is_pre_methylene_post_methylene_meta_aramid() ) {
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 3 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 3 ] - 180.0 );
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 4 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 4 ] - 180.0 );

	} else if ( rsd.is_para_aramid() ) {
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 2 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 2 ] - 180.0 );
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 3 ] - 0.0 ) * ( rsd.mainchain_torsions()[ 3 ] - 0.0 );
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 4 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 3 ] - 180.0 );
	} else if ( rsd.is_pre_methylene_para_aramid() ) {
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 3 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 3 ] - 180.0 );
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 4 ] - 0.0 ) * ( rsd.mainchain_torsions()[ 4 ] - 0.0 );
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 5 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 5 ] - 180.0 );
	} else if ( rsd.is_post_methylene_para_aramid() ) {
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 2 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 2 ] - 180.0 );
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 3 ] - 0.0 ) * ( rsd.mainchain_torsions()[ 3 ] - 0.0 );
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 4 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 3 ] - 180.0 );
	} else if ( rsd.is_pre_methylene_post_methylene_para_aramid() ) {
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 3 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 3 ] - 180.0 );
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 4 ] - 0.0 ) * ( rsd.mainchain_torsions()[ 4 ] - 0.0 );
		emap[ aromatic_restraint ] += ( rsd.mainchain_torsions()[ 5 ] - 180.0 ) * ( rsd.mainchain_torsions()[ 5 ] - 180.0 );

	} else if ( rsd.is_aramid() ) {
		TR.Error << "Case not accounted for: " << rsd.name() << std::endl;
	}
}


void
eval(
	Real const phi0,
	Real const weight,
	Size const rt1,
	Size const rt2,
	Size const rt3,
	Size const rt4,
	conformation::Residue const & rsd,
	utility::vector1< DerivVectorPair > & atom_derivs
) {
	Real phi = 0, dE_dphi;

	auto f1 = Vector( 0.0 );
	auto f2 = Vector( 0.0 );
	numeric::deriv::dihedral_p1_cosine_deriv(
		rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), rsd.xyz( rt4 ), phi, f1, f2 );
	Real del_phi = numeric::conversions::degrees( basic::subtract_radian_angles(phi, phi0) );
	dE_dphi = weight * 2 * del_phi;

	atom_derivs[ rt1 ].f1() += dE_dphi * f1;
	atom_derivs[ rt1 ].f2() += dE_dphi * f2;

	f1 = f2 = Vector(0.0);
	numeric::deriv::dihedral_p2_cosine_deriv(
		rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), rsd.xyz( rt4 ), phi, f1, f2 );
	atom_derivs[ rt2 ].f1() += dE_dphi * f1;
	atom_derivs[ rt2 ].f2() += dE_dphi * f2;

	f1 = f2 = Vector(0.0);
	numeric::deriv::dihedral_p2_cosine_deriv(
		rsd.xyz( rt4 ), rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), phi, f1, f2 );

	atom_derivs[ rt3 ].f1() += dE_dphi * f1;
	atom_derivs[ rt3 ].f2() += dE_dphi * f2;

	f1 = f2 = Vector(0.0);
	numeric::deriv::dihedral_p1_cosine_deriv(
		rsd.xyz( rt4 ), rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), phi, f1, f2 );
	atom_derivs[ rt4 ].f1() += dE_dphi * f1;
	atom_derivs[ rt4 ].f2() += dE_dphi * f2;
}

/// @brief Evaluate the derivatives for all atoms in this residue.
///
void
AromaticBackboneRestraintEnergy::eval_residue_derivatives(
	conformation::Residue const &rsd,
	ResSingleMinimizationData const & /*min_data*/,
	pose::Pose const & /*pose*/,
	EnergyMap const &weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const {

	if ( rsd.is_virtual_residue() ) return; //Skip virtual residues.

	auto weight = weights[ aromatic_restraint ];

	Real phi0 = 0;
	Size rt1, rt2, rt3, rt4;
	auto const & mca = rsd.mainchain_atoms();
	if ( rsd.is_ortho_aramid() ) {
		// eval for mca 2
		phi0 = 0;
		rt1 = mca[ 1 ];
		rt2 = mca[ 2 ];
		rt3 = mca[ 3 ];
		rt4 = mca[ 4 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );

	} else if ( rsd.is_pre_methylene_ortho_aramid() ) {
		// eval for mca 3
		phi0 = 0;
		rt1 = mca[ 2 ];
		rt2 = mca[ 3 ];
		rt3 = mca[ 4 ];
		rt4 = mca[ 5 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );

	} else if ( rsd.is_post_methylene_ortho_aramid() ) {
		// eval for mca 2
		phi0 = 0;
		rt1 = mca[ 1 ];
		rt2 = mca[ 2 ];
		rt3 = mca[ 3 ];
		rt4 = mca[ 4 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );

	} else if ( rsd.is_pre_methylene_post_methylene_ortho_aramid() ) {
		// eval for mca 3
		phi0 = 0;
		rt1 = mca[ 2 ];
		rt2 = mca[ 3 ];
		rt3 = mca[ 4 ];
		rt4 = mca[ 5 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
	} else if ( rsd.is_meta_aramid() ) {
		// eval for mca 2
		phi0 = 180;
		rt1 = mca[ 1 ];
		rt2 = mca[ 2 ];
		rt3 = mca[ 3 ];
		rt4 = mca[ 4 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
		rt1 = mca[ 2 ];
		rt2 = mca[ 3 ];
		rt3 = mca[ 4 ];
		rt4 = mca[ 5 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );

	} else if ( rsd.is_pre_methylene_meta_aramid() ) {
		// eval for mca 3
		phi0 = 180;
		rt1 = mca[ 2 ];
		rt2 = mca[ 3 ];
		rt3 = mca[ 4 ];
		rt4 = mca[ 5 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
		rt1 = mca[ 3 ];
		rt2 = mca[ 4 ];
		rt3 = mca[ 5 ];
		rt4 = mca[ 6 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );

	} else if ( rsd.is_post_methylene_meta_aramid() ) {
		// eval for mca 2
		phi0 = 180;
		rt1 = mca[ 1 ];
		rt2 = mca[ 2 ];
		rt3 = mca[ 3 ];
		rt4 = mca[ 4 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
		rt1 = mca[ 2 ];
		rt2 = mca[ 3 ];
		rt3 = mca[ 4 ];
		rt4 = mca[ 5 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );

	} else if ( rsd.is_pre_methylene_post_methylene_meta_aramid() ) {
		// eval for mca 3
		phi0 = 180;
		rt1 = mca[ 2 ];
		rt2 = mca[ 3 ];
		rt3 = mca[ 4 ];
		rt4 = mca[ 5 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
		rt1 = mca[ 3 ];
		rt2 = mca[ 4 ];
		rt3 = mca[ 5 ];
		rt4 = mca[ 6 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
	}  else if ( rsd.is_para_aramid() ) {
		// eval for mca 2
		phi0 = 180;
		rt1 = mca[ 1 ];
		rt2 = mca[ 2 ];
		rt3 = mca[ 3 ];
		rt4 = mca[ 4 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
		phi0 = 0;
		rt1 = mca[ 2 ];
		rt2 = mca[ 3 ];
		rt3 = mca[ 4 ];
		rt4 = mca[ 5 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
		phi0 = 180;
		rt1 = mca[ 3 ];
		rt2 = mca[ 4 ];
		rt3 = mca[ 5 ];
		rt4 = mca[ 6 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );

	} else if ( rsd.is_pre_methylene_para_aramid() ) {
		// eval for mca 3
		phi0 = 180;
		rt1 = mca[ 2 ];
		rt2 = mca[ 3 ];
		rt3 = mca[ 4 ];
		rt4 = mca[ 5 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
		phi0 = 0;
		rt1 = mca[ 3 ];
		rt2 = mca[ 4 ];
		rt3 = mca[ 5 ];
		rt4 = mca[ 6 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
		phi0 = 180;
		rt1 = mca[ 4 ];
		rt2 = mca[ 5 ];
		rt3 = mca[ 6 ];
		rt4 = mca[ 7 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );

	} else if ( rsd.is_post_methylene_para_aramid() ) {
		// eval for mca 2
		phi0 = 180;
		rt1 = mca[ 1 ];
		rt2 = mca[ 2 ];
		rt3 = mca[ 3 ];
		rt4 = mca[ 4 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
		phi0 = 0;
		rt1 = mca[ 2 ];
		rt2 = mca[ 3 ];
		rt3 = mca[ 4 ];
		rt4 = mca[ 5 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
		phi0 = 180;
		rt1 = mca[ 3 ];
		rt2 = mca[ 4 ];
		rt3 = mca[ 5 ];
		rt4 = mca[ 6 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );

	} else if ( rsd.is_pre_methylene_post_methylene_para_aramid() ) {
		// eval for mca 3
		phi0 = 180;
		rt1 = mca[ 2 ];
		rt2 = mca[ 3 ];
		rt3 = mca[ 4 ];
		rt4 = mca[ 5 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
		phi0 = 0;
		rt1 = mca[ 3 ];
		rt2 = mca[ 4 ];
		rt3 = mca[ 5 ];
		rt4 = mca[ 6 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
		phi0 = 180;
		rt1 = mca[ 4 ];
		rt2 = mca[ 5 ];
		rt3 = mca[ 6 ];
		rt4 = mca[ 7 ];
		eval( phi0, weight, rt1, rt2, rt3, rt4, rsd, atom_derivs );
	} else if ( rsd.is_aramid() ) {
		TR.Error << "Case not accounted for: " << rsd.name() << std::endl;
	}
	return;
}

/// @brief AromaticBackboneRestraint Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
AromaticBackboneRestraintEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
)
const
{}
core::Size
AromaticBackboneRestraintEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

