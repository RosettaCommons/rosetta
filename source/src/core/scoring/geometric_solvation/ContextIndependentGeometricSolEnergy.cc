// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ContextIndependentGeometricSolEnergy.cc
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)
/// @brief  Similar to the standard version of GeometricSolEnergy.cc BUT without the CONTEXT_DEPENDENT stuff. ALOT OF CODE DUPLICATION!
/// @brief  Significantly speed up when used in src/protocol/swa/rna/StepwiseRNA_Sampler.cc!

// Unit Headers
#include <core/scoring/geometric_solvation/ContextIndependentGeometricSolEnergy.hh>
#include <core/scoring/geometric_solvation/ContextIndependentGeometricSolEnergyCreator.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.hh>

// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunction.hh>

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


static basic::Tracer tr("core.scoring.geometric_solvation.ContextIndependentGeometricSolEnergy" );

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// Added on July. 22, 2011, Parin Sripakdeevong (sripakpa@stanford.edu).
//
// This copies huge amounts of code from GeometricSolEnergy.cc.
// Should instead make a GeometricSolPotential.cc, which holds *all* the core functions,
// and then GeometricSolEnergy and ContextIndependentGeometricSolEnergy can both call those core functions.
///////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace geometric_solvation {

/// @details This must return a fresh instance of the GeometricSolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ContextIndependentGeometricSolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new ContextIndependentGeometricSolEnergy( options );
}

ScoreTypes
ContextIndependentGeometricSolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( geom_sol_fast );
	sts.push_back( geom_sol_fast_intra_RNA );
	return sts;
}

///@brief copy c-tor
ContextIndependentGeometricSolEnergy::ContextIndependentGeometricSolEnergy( methods::EnergyMethodOptions const & opts) :
	parent( new ContextIndependentGeometricSolEnergyCreator ),
	options_( new methods::EnergyMethodOptions( opts ) ),
	evaluator_( new GeometricSolEnergyEvaluator( opts ) )
{
	if ( options_->hbond_options().use_hb_env_dep() ) {
		utility_exit_with_message( "Environment dependent hydrogen bonds are not compatible with geom_sol_fast. You can turn off hydrogen-bond environment dependence (e.g., with NO_HB_ENV_DEP in score file), and protein & RNA structure prediction won't get worse! (Or if you insist on using environment dependence, use geom_sol instead of geom_sol_fast.)" );
	}
}

/// copy ctor
ContextIndependentGeometricSolEnergy::ContextIndependentGeometricSolEnergy( ContextIndependentGeometricSolEnergy const & src ):
	ContextIndependentTwoBodyEnergy( src ),
	options_( new methods::EnergyMethodOptions( *src.options_ ) ),
	evaluator_( src.evaluator_ )
{
}

/// clone
methods::EnergyMethodOP
ContextIndependentGeometricSolEnergy::clone() const
{
	return new ContextIndependentGeometricSolEnergy( *this );
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

	//	A hack. Hopefully not a big slow down.
	EnergyMap emap_local;
	evaluator_->residue_pair_energy( rsd1, rsd2, pose, scorefxn, emap_local );
	emap[ geom_sol_fast ] += emap_local[ geom_sol ];

}


void
ContextIndependentGeometricSolEnergy::eval_atom_derivative(
		 id::AtomID const & atom_id,
		 pose::Pose const & pose,
		 kinematics::DomainMap const & domain_map,
		 ScoreFunction const & scorefxn,
		 EnergyMap const & weights,
		 Vector & F1,
		 Vector & F2
) const
{

	//	A hack. Hopefully not a big slow down.
	EnergyMap weights_local;
	weights_local[ geom_sol ]           = weights[ geom_sol_fast ];
	weights_local[ geom_sol_intra_RNA ] = weights[ geom_sol_fast_intra_RNA ];
	evaluator_->eval_atom_derivative( atom_id, pose, domain_map, scorefxn, weights_local, F1, F2 );

}

Distance
ContextIndependentGeometricSolEnergy::atomic_interaction_cutoff() const
{
	return evaluator_->atomic_interaction_cutoff();
}

///////////////////////////////////////////////////////////////////////////////////////
bool
ContextIndependentGeometricSolEnergy::defines_intrares_energy( EnergyMap const & weights ) const
{
	//Change to this on Feb 06, 2012. Ensure that the function returns false if weights[geom_sol_intra_RNA] == 0.0
	bool condition_1 = (weights[geom_sol_fast_intra_RNA] > 0.0001) ? true : false;
	return condition_1;
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

}

///@brief ContextIndependentGeometricSolEnergy is not context sensitive, of course.
void
ContextIndependentGeometricSolEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & context_graphs_required
) const
{
}


core::Size
ContextIndependentGeometricSolEnergy::version() const
{
	return 2; // Initial versioning
}


} // hbonds
} // scoring
} // core

