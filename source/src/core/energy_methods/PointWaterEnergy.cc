// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/energy_methods/PointWaterEnergy.hh
/// @brief  Statistical point water energy function
/// @author Frank DiMaio


// Unit headers
#include <core/energy_methods/PointWaterEnergy.hh>
#include <core/energy_methods/PointWaterEnergyCreator.hh>
#include <core/scoring/PointWaterPotential.hh>

// Package headers
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/scoring/DerivVectorPair.hh>

// Project headers
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/chemical/AA.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/deriv/distance_deriv.hh>

// Basic headers
#include <basic/Tracer.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

static basic::Tracer TR("core.scoring.methods.PointWaterEnergy");

namespace core {
namespace energy_methods {



/// @details This must return a fresh instance of the PointWaterEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
PointWaterEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &opts
) const {
	return utility::pointer::make_shared< PointWaterEnergy >( opts );
}

core::scoring::ScoreTypes
PointWaterEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( pointwater );
	return sts;
}


////////////////////////////////////////////////////////////////////////////
PointWaterEnergy::PointWaterEnergy( core::scoring::methods::EnergyMethodOptions const &opts ):
	parent( utility::pointer::make_shared< PointWaterEnergyCreator >() ),
	potential_( core::scoring::ScoringManager::get_instance()->get_PointWaterPotential() )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	pwater_ref_wt_ = opts.ordered_pt_wat_penalty();
	pwater_water_bonus_ = option[corrections::water::pointwat_wat_bonus]();
	pwater_water_bonus_width_ = option[corrections::water::pointwat_wat_bonus_width]();
}


////////////////////////////////////////////////////////////////////////////
PointWaterEnergy::PointWaterEnergy( PointWaterEnergy const & src ):
	parent( src ),
	potential_( src.potential_ ),
	pwater_ref_wt_( src.pwater_ref_wt_ )
{ }


core::scoring::methods::EnergyMethodOP
PointWaterEnergy::clone() const
{
	return utility::pointer::make_shared< PointWaterEnergy >( *this );
}

void
PointWaterEnergy::setup_for_minimizing(
	pose::Pose & /*pose*/,
	core::scoring::ScoreFunction const & /*sfxn*/,
	kinematics::MinimizerMapBase const & /*min_map*/
) const
{ }


///
void
PointWaterEnergy::setup_for_derivatives( pose::Pose & /*pose*/, core::scoring::ScoreFunction const & /*sfxn*/ ) const
{ }

///
void
PointWaterEnergy::setup_for_scoring( pose::Pose & /*pose*/, core::scoring::ScoreFunction const & /*sfxn*/ ) const
{ }


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
PointWaterEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & ,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	bool water_water =
		(rsd1.aa() == core::chemical::aa_h2o && rsd2.aa() == core::chemical::aa_h2o );
	bool water_prot = !water_water && (rsd1.aa() == core::chemical::aa_h2o || rsd2.aa() == core::chemical::aa_h2o);

	core::Real score = 0;
	if ( water_prot ) {
		core::conformation::Residue const &r_protein = (rsd1.aa() == core::chemical::aa_h2o) ? rsd2 : rsd1;
		core::conformation::Residue const &r_water   = (rsd1.aa() == core::chemical::aa_h2o) ? rsd1 : rsd2;

		if ( r_water.name() == "PWAT" || r_water.name() == "BB_PWAT" || r_water.name() == "PWAT1" ) {
			score = potential_.eval_pointwater_score( r_protein.aa(), r_protein, r_water.atom(1).xyz() );
		}
	} else if ( water_water && ((rsd1.name() == "PWAT" && rsd2.name() == "PWAT") || (rsd1.name() == "BB_PWAT" && rsd2.name() == "BB_PWAT") ||
			(rsd1.name() == "PWAT1" && rsd2.name() == "PWAT1") ) ) {
		// bonus for good water-water distances
		core::Vector wat1 = rsd1.atom(1).xyz();
		core::Vector wat2 = rsd2.atom(1).xyz();
		core::Real distance = (wat1-wat2).length();
		score = pwater_water_bonus_*(-exp( -(distance-2.7)*(distance-2.7)*pwater_water_bonus_width_ ));
	}

	emap[ core::scoring::pointwater ] += score;
}

void
PointWaterEnergy::eval_intrares_energy(
	conformation::Residue const &rsd1,
	pose::Pose const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap &emap
) const {
	if ( rsd1.name() == "PWAT" || rsd1.name() == "BB_PWAT" || rsd1.name() == "PWAT1" ) {
		emap[ core::scoring::pointwater ] += pwater_ref_wt_;
	}
}

bool
PointWaterEnergy::defines_score_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	bool res_moving_wrt_eachother
) const
{
	bool water_prot =
		(rsd1.aa() == core::chemical::aa_h2o && rsd2.is_protein() ) ||
		(rsd2.aa() == core::chemical::aa_h2o && rsd1.is_protein() );
	bool water_water =
		(rsd1.aa() == core::chemical::aa_h2o && rsd2.aa() == core::chemical::aa_h2o );

	if ( !water_prot && !water_water ) {
		return false;
	}

	return res_moving_wrt_eachother;
}


void
PointWaterEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	core::scoring::ResSingleMinimizationData const &,
	core::scoring::ResSingleMinimizationData const &,
	core::scoring::ResPairMinimizationData const & ,
	pose::Pose const & , // provides context
	core::scoring::EnergyMap const & weights,
	utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
	utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
) const
{
	bool water_prot =
		(rsd1.aa() == core::chemical::aa_h2o && rsd2.is_protein() ) ||
		(rsd2.aa() == core::chemical::aa_h2o && rsd1.is_protein() );
	bool water_water =
		(rsd1.aa() == core::chemical::aa_h2o && rsd2.aa() == core::chemical::aa_h2o );

	if ( water_prot ) {
		core::conformation::Residue const &r_protein = (rsd1.aa() == core::chemical::aa_h2o) ? rsd2 : rsd1;
		core::conformation::Residue const &r_water   = (rsd1.aa() == core::chemical::aa_h2o) ? rsd1 : rsd2;
		utility::vector1< core::scoring::DerivVectorPair > &deriv_protein = (rsd1.aa() == core::chemical::aa_h2o) ? r2_atom_derivs : r1_atom_derivs;
		utility::vector1< core::scoring::DerivVectorPair > &deriv_water   = (rsd1.aa() == core::chemical::aa_h2o) ? r1_atom_derivs : r2_atom_derivs;

		if ( r_water.name() == "PWAT" || r_water.name() == "BB_PWAT" || r_water.name() == "PWAT1" ) {
			potential_.eval_pointwater_derivs( r_protein.aa(), r_protein, r_water.atom(1).xyz(), deriv_protein, deriv_water, weights[ core::scoring::pointwater ] );
		}
	} else if ( water_water && ((rsd1.name() == "PWAT" && rsd2.name() == "PWAT") || (rsd1.name() == "BB_PWAT" && rsd2.name() == "BB_PWAT") ||
			(rsd1.name() == "PWAT1" && rsd2.name() == "PWAT1")) ) {
		// bonus for good water-water distances
		core::Vector wat1 = rsd1.atom(1).xyz();
		core::Vector wat2 = rsd2.atom(1).xyz();
		core::Real distance = (wat1-wat2).length();
		core::Real dEdd = weights[ core::scoring::pointwater ]*2*pwater_water_bonus_width_*pwater_water_bonus_*exp( -(distance-2.7)*(distance-2.7)*pwater_water_bonus_width_ )*(distance-2.7);

		Vector f1(0.0), f2(0.0);
		numeric::deriv::distance_f1_f2_deriv( wat1, wat2, distance, f1, f2 );
		r1_atom_derivs[ 1 ].f1() += dEdd * f1;
		r1_atom_derivs[ 1 ].f2() += dEdd * f2;
		r2_atom_derivs[ 1 ].f1() -= dEdd * f1;
		r2_atom_derivs[ 1 ].f2() -= dEdd * f2;
	}

}


/// @brief PointWaterEnergy distance cutoff
core::Size
PointWaterEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace energy_methods
} // namespace core
