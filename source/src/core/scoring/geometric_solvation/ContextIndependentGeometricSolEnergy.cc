// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ContextIndependentGeometricSolEnergy.cc
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)
/// @brief  Similar to the standard version of GeometricSolEnergy.cc BUT without the CONTEXT_DEPENDENT stuff.
/// @brief  Significantly speed up when used in src/protocol/swa/rna/StepwiseRNA_Sampler.cc

// Unit Headers
#include <core/scoring/geometric_solvation/ContextIndependentGeometricSolEnergy.hh>
#include <core/scoring/geometric_solvation/ContextIndependentGeometricSolEnergyCreator.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.hh>

// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/kinematics/MinimizerMapBase.hh>

#include <core/id/AtomID.hh>

// Project headers
#include <core/conformation/Residue.hh>

// Utility headers
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>


static THREAD_LOCAL basic::Tracer TR( "core.scoring.geometric_solvation.ContextIndependentGeometricSolEnergy" );

///////////////////////////////////////////////////////////////////
// See notes in GeometricSolEnergyEvaluator for what this does.
//
// If you make changes in here, please make sure to also update
//  ContextDependentGeometricSolEnergy, ideally but putting a shared
//  function in GeometricSolEnergyEvaluator to prevent code copying.
//
///////////////////////////////////////////////////////////////////
namespace core {
namespace scoring {
namespace geometric_solvation {

/// @details This must return a fresh instance of the GeometricSolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ContextIndependentGeometricSolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new ContextIndependentGeometricSolEnergy( options ) );
}

ScoreTypes
ContextIndependentGeometricSolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( geom_sol_fast );
	sts.push_back( geom_sol_fast_intra_RNA );
	return sts;
}

/// @brief copy c-tor
ContextIndependentGeometricSolEnergy::ContextIndependentGeometricSolEnergy( methods::EnergyMethodOptions const & opts) :
	parent( methods::EnergyMethodCreatorOP( new ContextIndependentGeometricSolEnergyCreator ) ),
	options_( opts ),
	evaluator_( GeometricSolEnergyEvaluatorOP( new GeometricSolEnergyEvaluator( opts ) ) ),
	precalculated_bb_bb_energy_(0.0f),
	using_extended_method_(false)
{
	if ( options_.hbond_options().use_hb_env_dep() ) {
		utility_exit_with_message( "Environment dependent hydrogen bonds are not compatible with geom_sol_fast. You can turn off hydrogen-bond environment dependence (e.g., with NO_HB_ENV_DEP in score file), and protein & RNA structure prediction won't get worse! (Or if you insist on using environment dependence, use geom_sol instead of geom_sol_fast.)" );
	}
}

/// copy ctor
ContextIndependentGeometricSolEnergy::ContextIndependentGeometricSolEnergy( ContextIndependentGeometricSolEnergy const & src ):
	ContextIndependentTwoBodyEnergy( src ),
	options_( src.options_ ),
	evaluator_( src.evaluator_ ),
	precalculated_bb_bb_energy_(src.precalculated_bb_bb_energy_),
	using_extended_method_(src.using_extended_method_)
{
}

/// clone
methods::EnergyMethodOP
ContextIndependentGeometricSolEnergy::clone() const
{
	return methods::EnergyMethodOP( new ContextIndependentGeometricSolEnergy( *this ) );
}


void
ContextIndependentGeometricSolEnergy::precalculate_bb_bb_energy_for_design(
	pose::Pose const &
) const {
}

void
ContextIndependentGeometricSolEnergy::setup_for_packing(
	pose::Pose  & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const & designing_residues
) const
{
	bool might_be_designing = std::any_of(
		designing_residues.begin(), designing_residues.end(),
		[]( bool const b ){ return b; } );

	precalculated_bb_bb_energy_ = 0.0f;

	//if nothing can be designed no reason to precalculate backbone/backbone geom_solv
	if ( !might_be_designing ) return;
	precalculated_bb_bb_energy_ = evaluator_->precalculate_bb_bb_energy_for_design( pose );
}


void
ContextIndependentGeometricSolEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & /*scorefxn*/ ) const
{
	pose.update_residue_neighbors();
}

///////////////////////////////////////////////////////////////////////////////
bool
ContextIndependentGeometricSolEnergy::defines_score_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	bool res_moving_wrt_eachother
) const
{
	return evaluator_->defines_score_for_residue_pair( rsd1, rsd2, res_moving_wrt_eachother);
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
ContextIndependentGeometricSolEnergy::get_count_pair_function(
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
ContextIndependentGeometricSolEnergy::get_count_pair_function(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const
{
	return evaluator_->get_count_pair_function( rsd1, rsd2 );
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
ContextIndependentGeometricSolEnergy::get_intrares_countpair(
	conformation::Residue const & res,
	pose::Pose const &,
	ScoreFunction const &
) const
{
	return evaluator_->get_intrares_countpair( res );
}

///////////////////////////////////////////////////////////////////////////////
bool
ContextIndependentGeometricSolEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}

///////////////////////////
void
ContextIndependentGeometricSolEnergy::setup_for_minimizing_for_residue(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData &
) const
{}


////////////////////////////////////////////////////////////////////////////////
void
ContextIndependentGeometricSolEnergy::setup_for_minimizing_for_residue_pair(
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

/////////////
bool
ContextIndependentGeometricSolEnergy::requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & ) const
{
	return false;
}

////////////////////
void
ContextIndependentGeometricSolEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	if ( !defines_intrares_energy( weights ) ) return;
	evaluator_->eval_intrares_derivatives( rsd, pose, weights[ geom_sol_fast_intra_RNA ], atom_derivs, true /*just_RNA*/ );
	if ( options_.put_intra_into_total() ) evaluator_->eval_intrares_derivatives( rsd, pose, weights[ geom_sol_fast ], atom_derivs, false /*just_RNA*/ );
}

////////////////////////////////////////////////////
void
ContextIndependentGeometricSolEnergy::eval_residue_pair_derivatives(
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
	evaluator_->eval_residue_pair_derivatives( ires, jres, min_data, pose, weights[ geom_sol_fast], r1_atom_derivs, r2_atom_derivs );
}

////////////////////////////////////////////////////////////////////////////////
void
ContextIndependentGeometricSolEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	using_extended_method_ = true;
	emap[ geom_sol_fast ] += evaluator_->residue_pair_energy_ext( rsd1, rsd2, min_data, pose );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
/// Everything in here.
void
ContextIndependentGeometricSolEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap
) const
{
	using_extended_method_ = false;

	//if the backbone/backbone energy has already been calculated in setup_for_packing
	//  this is done only if we are doing fix backbone design and the backbone/backbone
	//  energy cannot change (Joseph Yesselman 9/11/13)
	// note from rhiju -- this looks dangerous
	if ( precalculated_bb_bb_energy_ > 0.0f ) {
		emap[ geom_sol_fast ] += evaluator_->geometric_sol_one_way_sc(rsd1, rsd2, pose) +
			evaluator_->geometric_sol_one_way_sc(rsd2, rsd1, pose);
	} else {
		EnergyMap emap_local;
		evaluator_->residue_pair_energy( rsd1, rsd2, pose, scorefxn, emap_local );
		emap[ geom_sol_fast ] += emap_local[ geom_sol ];
	}
}

bool
ContextIndependentGeometricSolEnergy::minimize_in_whole_structure_context( pose::Pose const & ) const
{
	return false;
}

void
ContextIndependentGeometricSolEnergy::finalize_total_energy(
	pose::Pose & /*pose*/,
	ScoreFunction const &,
	EnergyMap & totals
) const
{
	// note from rhiju -- this looks dangerous -- check and see how it is done elsewhere.
	if ( !using_extended_method_ ) {
		totals[ geom_sol_fast ] += precalculated_bb_bb_energy_;
		return;
	}
}

Distance
ContextIndependentGeometricSolEnergy::atomic_interaction_cutoff() const
{
	return evaluator_->atomic_interaction_cutoff();
}

///////////////////////////////////////////////////////////////////////////////////////
bool
ContextIndependentGeometricSolEnergy::defines_intrares_energy( EnergyMap const &  ) const
{
	//Change to this on Feb 06, 2012. Ensure that the function returns false if weights[geom_sol_intra_RNA] == 0.0
	// return ( weights[geom_sol_fast_intra_RNA] > 0.0 || options_.include_intra() );
	return true; // always calculate.
}


///////////////////////////////////////////////////////////////////////////////////////
void
ContextIndependentGeometricSolEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap
) const{

	EnergyMap emap_local;
	evaluator_->eval_intrares_energy( rsd, pose, scorefxn, emap_local );
	emap[ geom_sol_fast_intra_RNA ] += emap_local[ geom_sol_intra_RNA ];
	emap[ geom_sol_fast ]           += emap_local[ geom_sol ];
}

/// @brief ContextIndependentGeometricSolEnergy is not context sensitive, of course.
void
ContextIndependentGeometricSolEnergy::indicate_required_context_graphs(utility::vector1< bool > &) const
{}

core::Size
ContextIndependentGeometricSolEnergy::version() const
{
	return 3;
}


} // geometric_solvation
} // scoring
} // core

