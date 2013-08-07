// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/GeometricSolEnergy.fwd.hh
/// @brief  Geometric solvation energy.
/// @author Rhiju Das

// Unit Headers
#include <core/scoring/geometric_solvation/GeometricSolEnergy.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyCreator.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.hh>

// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/TenANeighborGraph.hh>

// Project headers

// Utility headers
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/prof.hh>

#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/StringVectorOption.hh>
#include <ObjexxFCL/FArray3D.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>

static basic::Tracer tr( "core.scoring.geometric_solvation.GeometricSolEnergy" );

//////////////////////////////////////////////////////////////////
// All the good stuff is now in GeometricSolEnergyEvaluator, which is
// shared with (faster) ContextIndependentGeometricSolEnergy.
//////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace geometric_solvation {

/// @details This must return a fresh instance of the GeometricSolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
GeometricSolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new GeometricSolEnergy( options );
}

ScoreTypes
GeometricSolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( geom_sol );
	sts.push_back( geom_sol_intra_RNA );
	return sts;
}


///@brief copy c-tor
GeometricSolEnergy::GeometricSolEnergy( methods::EnergyMethodOptions const & opts
) :
	parent( new GeometricSolEnergyCreator ),
	options_( new methods::EnergyMethodOptions( opts ) ),
	evaluator_( new GeometricSolEnergyEvaluator( opts ) )
{
}

/// copy ctor
GeometricSolEnergy::GeometricSolEnergy( GeometricSolEnergy const & src ):
	ContextDependentTwoBodyEnergy( src ),
	options_( new methods::EnergyMethodOptions( *src.options_ ) ),
	evaluator_( src.evaluator_ )
{}

/// clone
methods::EnergyMethodOP
GeometricSolEnergy::clone() const
{
	return new GeometricSolEnergy( *this );
}

///
void
GeometricSolEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	using core::scoring::EnergiesCacheableDataType::HBOND_SET;

	// We need the H-bond set -- well, at least the backbone/backbone h-bonds
	// when computing geometric solvation scores.
	// Since this is probably being computed elsewhere, might make sense
	// to have a "calculated" flag.
	// But, anyway, the geometric sol calcs take way longer than this.
	pose.update_residue_neighbors();

	hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( options_->hbond_options() ) );
	hbond_set->setup_for_residue_pair_energies( pose );
	pose.energies().data().set( HBOND_SET, hbond_set );
}

/////////////////////////////////////////////////////////////////////////////
void
GeometricSolEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	using core::scoring::EnergiesCacheableDataType::HBOND_SET;

// same setup as for HBondEnergy.cc. Note that this is probably repeating some work
// that has already occurred in hbonds calculation. Probably could have a "calculated"
// flag to save the computation ... but the geometric sol. calculation takes so
// much longer that this initial hbond loop isn't too bad.
 	pose.update_residue_neighbors();

	hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( options_->hbond_options() ) );
	hbonds::fill_hbond_set( pose, true /*calc derivs*/, *hbond_set );
	hbond_set->resize_bb_donor_acceptor_arrays( pose.total_residue() );
	pose.energies().data().set( HBOND_SET, hbond_set );
}


//void
//GeometricSolEnergy::setup_for_derivatives( pose::Pose & /*pose*/, ScoreFunction const & ) const
//{
//}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// Everything in here.
void
GeometricSolEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap ) const
{
	evaluator_->residue_pair_energy( rsd1, rsd2, pose, scorefxn, emap );
}

void
GeometricSolEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const & scorefxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	evaluator_->eval_atom_derivative( atom_id, pose, domain_map, scorefxn, weights, F1, F2 );
}

Distance
GeometricSolEnergy::atomic_interaction_cutoff() const
{
	return evaluator_->atomic_interaction_cutoff();
}

///////////////////////////////////////////////////////////////////////////////////////
bool
GeometricSolEnergy::defines_intrares_energy( EnergyMap const & weights ) const
{
	//Change to this on Feb 06, 2012. Ensure that the function returns false if weights[geom_sol_intra_RNA] == 0.0
	bool condition_1 = (weights[geom_sol_intra_RNA] > 0.0001) ? true : false;
	return condition_1;
}

///////////////////////////////////////////////////////////////////////////////////////
void
GeometricSolEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap ) const
{
	evaluator_->eval_intrares_energy( rsd, pose, scorefxn, emap );
}

///////////////////////////////////////////////////////////////////////////////////////
///@brief GeometricSolEnergy is context sensitive
void
GeometricSolEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & context_graphs_required
) const
{
	context_graphs_required[ ten_A_neighbor_graph ] = true;
}

core::Size
GeometricSolEnergy::version() const
{
	return 2; // Initial versioning
}


} // hbonds
} // scoring
} // core

