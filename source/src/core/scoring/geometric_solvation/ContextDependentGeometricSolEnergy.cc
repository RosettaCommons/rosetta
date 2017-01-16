// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ContextDependentGeometricSolEnergy.fwd.hh
/// @brief  Geometric solvation energy, with environment dependence allowed. Note ContextIndependent version is faster.
/// @author Rhiju Das

// Unit Headers
#include <core/scoring/geometric_solvation/ContextDependentGeometricSolEnergy.hh>
#include <core/scoring/geometric_solvation/ContextDependentGeometricSolEnergyCreator.hh>
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
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/types.hh>


// Project headers

// Utility headers
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/StringVectorOption.hh>
#include <ObjexxFCL/FArray3D.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>

static THREAD_LOCAL basic::Tracer tr( "core.scoring.geometric_solvation.GeometricSolEnergy" );

//////////////////////////////////////////////////////////////////
// All the good stuff is now in GeometricSolEnergyEvaluator, which is
// shared with (faster) ContextIndependentGeometricSolEnergy.
//
// See notes in GeometricSolEnergyEvaluator for what this does.
//
// If you make changes in here, please make sure to also update
//  ContextIndependentGeometricSolEnergy, ideally but putting a shared
//  function in GeometricSolEnergyEvaluator to prevent code copying.
//
//////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace geometric_solvation {

/// @details This must return a fresh instance of the GeometricSolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ContextDependentGeometricSolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new ContextDependentGeometricSolEnergy( options ) );
}

ScoreTypes
ContextDependentGeometricSolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( geom_sol );
	sts.push_back( geom_sol_intra_RNA );
	return sts;
}


/// @brief copy c-tor
ContextDependentGeometricSolEnergy::ContextDependentGeometricSolEnergy( methods::EnergyMethodOptions const & opts
) :
	parent( methods::EnergyMethodCreatorOP( new ContextDependentGeometricSolEnergyCreator ) ),
	options_( opts ),
	evaluator_( GeometricSolEnergyEvaluatorOP( new GeometricSolEnergyEvaluator( opts ) ) ),
	precalculated_bb_bb_energy_(0.0f),
	using_extended_method_( false )
{
}

/// copy ctor
ContextDependentGeometricSolEnergy::ContextDependentGeometricSolEnergy( ContextDependentGeometricSolEnergy const & src ):
	ContextDependentTwoBodyEnergy( src ),
	options_( src.options_ ),
	evaluator_( src.evaluator_ ),
	precalculated_bb_bb_energy_(src.precalculated_bb_bb_energy_),
	using_extended_method_( src.using_extended_method_ )
{}

/// clone
methods::EnergyMethodOP
ContextDependentGeometricSolEnergy::clone() const
{
	return methods::EnergyMethodOP( new ContextDependentGeometricSolEnergy( *this ) );
}

void
ContextDependentGeometricSolEnergy::setup_for_packing(
	pose::Pose  & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const & designing_residues
) const
{
	// Need an HBOND_SET if we don't already have one from regular HBONDS
	// Note: The logic here is kinda hacky.
	using EnergiesCacheableDataType::HBOND_SET;
	if ( ! pose.energies().data().has( HBOND_SET ) ) {
		hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( options_.hbond_options() ) );
		hbond_set->setup_for_residue_pair_energies( pose );
		pose.energies().data().set( HBOND_SET, hbond_set );
	}

	bool might_be_designing = std::any_of(
		designing_residues.begin(), designing_residues.end(),
		[]( bool const b ){ return b; } );
	precalculated_bb_bb_energy_ = 0.0f;

	//if nothing can be designed no reason to precalculate backbone/backbone geom_solv
	if ( !might_be_designing ) return;
	precalculated_bb_bb_energy_ = evaluator_->precalculate_bb_bb_energy_for_design( pose );
}


void
ContextDependentGeometricSolEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();

	// Need an HBOND_SET if we don't already have one from regular HBONDS
	// Note: The logic here is kinda hacky.
	using EnergiesCacheableDataType::HBOND_SET;
	if ( ! pose.energies().data().has( HBOND_SET ) ) {
		hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( options_.hbond_options() ) );
		hbond_set->setup_for_residue_pair_energies( pose );
		pose.energies().data().set( HBOND_SET, hbond_set );
	}
}

///////////////////////////////////////////////////////////////////////////////
bool
ContextDependentGeometricSolEnergy::defines_score_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	bool res_moving_wrt_eachother
) const
{
	return evaluator_->defines_score_for_residue_pair( rsd1, rsd2, res_moving_wrt_eachother);
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
ContextDependentGeometricSolEnergy::get_count_pair_function(
	Size const res1,
	Size const res2,
	pose::Pose const & pose,
	ScoreFunction const &
) const
{
	return evaluator_->get_count_pair_function( res1, res2, pose );
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
ContextDependentGeometricSolEnergy::get_count_pair_function(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const
{
	return evaluator_->get_count_pair_function( rsd1, rsd2 );
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
ContextDependentGeometricSolEnergy::get_intrares_countpair(
	conformation::Residue const & res,
	pose::Pose const &,
	ScoreFunction const &
) const
{
	return evaluator_->get_intrares_countpair( res );
}

///////////////////////////////////////////////////////////////////////////////
bool
ContextDependentGeometricSolEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}


/////////////
bool
ContextDependentGeometricSolEnergy::requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & ) const
{
	return false;
}

////////////////////
void
ContextDependentGeometricSolEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	if ( !defines_intrares_energy( weights ) ) return;
	evaluator_->eval_intrares_derivatives( rsd, pose, weights[ geom_sol_intra_RNA ], atom_derivs, true /*just_RNA*/);
	if ( options_.put_intra_into_total() ) evaluator_->eval_intrares_derivatives( rsd, pose, weights[ geom_sol ], atom_derivs, false /*just_RNA*/ );
}

////////////////////////////////////////////////////
void
ContextDependentGeometricSolEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & ires,
	conformation::Residue const & jres,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	evaluator_->eval_residue_pair_derivatives( ires, jres, min_data, pose, weights[ geom_sol ], r1_atom_derivs, r2_atom_derivs );
}

////////////////////////////////////////////////////////////////////////////////
void
ContextDependentGeometricSolEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	using_extended_method_ = true;
	emap[ geom_sol ] += evaluator_->residue_pair_energy_ext( rsd1, rsd2, min_data, pose );
}

////////////////////////////////////////////////////////////////////////////////
void
ContextDependentGeometricSolEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData & pair_data
) const
{
	evaluator_->setup_for_minimizing_for_residue_pair( rsd1, rsd2, pair_data );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
/// Everything in here.
void
ContextDependentGeometricSolEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap ) const
{
	using_extended_method_ = false;

	//if the backbone/backbone energy has already been calculated in setup_for_packing
	//  this is done only if we are doing fix backbone design and the backbone/backbone
	//  energy cannot change (Joseph Yesselman 9/11/13)
	// note from rhiju -- this looks dangerous
	if ( precalculated_bb_bb_energy_ > 0.0f ) {
		emap[ geom_sol ] += evaluator_->geometric_sol_one_way_sc(rsd1, rsd2, pose) +
			evaluator_->geometric_sol_one_way_sc(rsd2, rsd1, pose);
	} else {
		evaluator_->residue_pair_energy( rsd1, rsd2, pose, scorefxn, emap );
	}
}

/////////////////////////////
bool
ContextDependentGeometricSolEnergy::minimize_in_whole_structure_context( pose::Pose const & ) const
{
	return false; //pose.energies().use_nblist_auto_update(); //???
}


void
ContextDependentGeometricSolEnergy::finalize_total_energy(
	pose::Pose &,
	ScoreFunction const &,
	EnergyMap & totals ) const
{
	if ( !using_extended_method_ ) {
		totals[ geom_sol ] += precalculated_bb_bb_energy_;
	}
}


Distance
ContextDependentGeometricSolEnergy::atomic_interaction_cutoff() const
{
	return evaluator_->atomic_interaction_cutoff();
}

///////////////////////////////////////////////////////////////////////////////////////
bool
ContextDependentGeometricSolEnergy::defines_intrares_energy( EnergyMap const &  ) const
{
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////
void
ContextDependentGeometricSolEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap ) const
{
	EnergyMap emap_local;
	evaluator_->eval_intrares_energy( rsd, pose, scorefxn, emap );
	emap[ geom_sol_intra_RNA ] += emap_local[ geom_sol_intra_RNA ];
	emap[ geom_sol ]           += emap_local[ geom_sol ];
}


///////////////////////////////////////////////////////////////////////////////////////
/// @brief GeometricSolEnergy is context sensitive
void
ContextDependentGeometricSolEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & context_graphs_required
) const
{
	context_graphs_required[ ten_A_neighbor_graph ] = true;
}

core::Size
ContextDependentGeometricSolEnergy::version() const
{
	return 3; // Initial versioning
}


} // hbonds
} // scoring
} // core

