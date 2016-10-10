// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/methods/carbohydrates/SugarBackboneEnergy.hh
/// @brief   Method definitions for SugarBackboneEnergy and SugarBackboneEnergyCreator.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit Headers
#include <core/scoring/methods/carbohydrates/SugarBackboneEnergy.hh>
#include <core/scoring/methods/carbohydrates/SugarBackboneEnergyCreator.hh>

// Package Headers
#include <core/scoring/carbohydrates/CHIEnergyFunction.hh>
#include <core/scoring/carbohydrates/OmegaPreferencesFunction.hh>
#include <core/scoring/carbohydrates/util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>

// Project Headers
#include <core/id/TorsionID.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>

// Numeric Headers
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/angle.functions.hh>

// Basic Header
#include <basic/Tracer.hh>


// Construct tracer.
static THREAD_LOCAL basic::Tracer TR( "core.scoring.methods.carbohydrates.SugarBackboneEnergy" );


namespace core {
namespace scoring {
namespace methods {
namespace carbohydrates {

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
SugarBackboneEnergy::SugarBackboneEnergy() :
		ContextIndependentOneBodyEnergy( EnergyMethodCreatorOP( new SugarBackboneEnergyCreator ) ),
		E_cef_( ScoringManager::get_instance()->get_CHIEnergyFunction() ),
		E_opf_( ScoringManager::get_instance()->get_OmegaPreferencesFunction() )
{}


// General EnergyMethod Methods ///////////////////////////////////////////////
EnergyMethodOP
SugarBackboneEnergy::clone() const
{
	return EnergyMethodOP( new SugarBackboneEnergy );
}


// OneBodyEnergy Methods //////////////////////////////////////////////////////
// Evaluate the one-body carbohydrate backbone energies for a particular residue.
void
SugarBackboneEnergy::residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap ) const
{
	using namespace numeric;
	using namespace chemical::carbohydrates;
	using namespace scoring::carbohydrates;
	using namespace pose::carbohydrates;

	// This is a carbohydrate-only scoring method.
	if ( ! rsd.is_carbohydrate() ) { return; }

	// Phi, psi, and omega are meaningless for reducing-end sugars.
	if ( rsd.is_lower_terminus() ) { return; }

	// Ignore REPLONLY variants.
	if ( rsd.has_variant_type( chemical::REPLONLY ) ) { return; }

	// Get parent residue, CarbohydrateInfo, and exocyclic state for convenience.
	conformation::Residue const & prev_rsd( pose.residue( find_seqpos_of_saccharides_parent_residue( rsd ) ) );
	CarbohydrateInfoCOP info( rsd.carbohydrate_info() );
	bool const is_exocyclic_bond( has_exocyclic_glycosidic_linkage( rsd, prev_rsd ) );

	// Get the angles.
	// (Convert the psi and omega to between 0 and 360 because that's what the functions expect.)
	uint const seqpos( rsd.seqpos() );
	Angle phi( pose.phi( seqpos ) );
	Angle psi( nonnegative_principal_angle_degrees( pose.psi( seqpos ) ) );
	Angle omega( nonnegative_principal_angle_degrees( pose.omega( seqpos ) ) );  // 1st omega, technically

	Energy score( 0.0 );


	// Calculate phi component. ///////////////////////////////////////////////
	// L-Sugars use the mirror image of the score functions.
	if ( info->is_L_sugar() ) {
		phi = -phi;
	}

	// Score is 0 if linkage type is LINKAGE_NA.
	score += E_cef_( get_CHI_energy_function_linkage_type_for_phi_for_residue_in_pose( pose, seqpos ), phi );


	// Calculate psi component. ///////////////////////////////////////////////
	// Psis TO L-Sugars use the mirror image of the score functions, if NOT exocyclic.
	// Psis FROM L-Sugars use the mirror image of the score functions, IF exocyclic
	if ( prev_rsd.is_carbohydrate() ) {
		if ( ! is_exocyclic_bond ) {
			if ( prev_rsd.carbohydrate_info()->is_L_sugar() ) {
				psi = 360 - psi;
			}
		}
	} else if ( is_exocyclic_bond ) {
		if ( info->is_L_sugar() ) {
			psi = 360 - psi;
		}
	}

	// Score is 0 if linkage type is LINKAGE_NA.
	score += E_cef_( get_CHI_energy_function_linkage_type_for_psi_for_residue_in_pose( pose, seqpos ), psi );


	// Calculate omega component. /////////////////////////////////////////////
	// Omegas TO L-Sugars use the mirror image of the score functions.
	if ( prev_rsd.is_carbohydrate() ) {
		if ( ! is_exocyclic_bond ) {
			if ( prev_rsd.carbohydrate_info()->is_L_sugar() ) {
				omega = 360 - psi;
			}
		}
	}

	// Score is 0 if linkage type is PREFERENCE_NA.
	score += E_opf_( get_omega_preference_for_residue_in_pose( pose, seqpos ), omega );


	emap[ sugar_bb ] += score;
}

// Evaluate the DoF derivative for a particular residue.
core::Real
SugarBackboneEnergy::eval_residue_dof_derivative(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & /* min_data */,
		id::DOF_ID const & /* dof_id */,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose,
		ScoreFunction const & /* sf */,
		EnergyMap const & weights ) const
{
	using namespace std;
	using namespace numeric;
	using namespace numeric::constants::d;
	using namespace id;
	using namespace chemical::carbohydrates;
	using namespace pose::carbohydrates;
	using namespace scoring::carbohydrates;

	TR.Debug << "Evaluating torsion: " << torsion_id << endl;

	Real deriv( 0.0 );

	// Ignore REPLONLY variants.
	if ( rsd.has_variant_type( chemical::REPLONLY ) ) { return deriv; }

	// This scoring method only considers glycosidic torsions, which may have either BB, CHI, or BRANCH TorsionIDs.
	if ( ! torsion_id.valid() ) {
		TR.Debug << "Torsion not valid: " << torsion_id << endl;
		return deriv;
	}

	if ( is_glycosidic_phi_torsion( pose, torsion_id ) ) {
		// Rosetta defines this torsion angle differently than IUPAC.  In addition, this torsion belongs to the next
		// residue.  We need to figure out what that residue is, use pose.phi() to get the value, and see if it's
		// alpha or beta.

		// First, what is the next residue?
		uint const next_rsd( get_downstream_residue_that_this_torsion_moves( pose, torsion_id ) );

		// Now, get the next residue's info and its phi.
		CarbohydrateInfoCOP info( pose.residue( next_rsd ).carbohydrate_info() );
		Angle phi( pose.phi( next_rsd ) );
		TR.Debug << "Phi: " << phi << endl;

		if ( info->is_L_sugar() ) {
			// L-Sugars use the mirror image of the score functions. The phi functions run from -180 to 180.
			phi = -phi;
		}

		// Finally, we can evaluate.
		deriv = E_cef_.evaluate_derivative(
				get_CHI_energy_function_linkage_type_for_phi_for_residue_in_pose( pose, next_rsd ), phi );

		if ( info->is_L_sugar() ) {
			deriv = -deriv;
		}
	} else if ( is_glycosidic_psi_torsion( pose, torsion_id ) ) {
		// Rosetta defines this torsion angle the same way as IUPAC; however, it is the psi of the following residue.
		// This should not matter, though, because the CHI energy function only cares about whether the linkage is
		// axial or equatorial.  We can simply convert to a number in degrees between 0 and 360 and call the function.

		// First, make sure this is a sugar.
		if ( ! rsd.is_carbohydrate() ) { return deriv; }

		// Second, get the psi and convert it to between 0 and 360 (because that's what the function expects).
		Angle psi( nonnegative_principal_angle_degrees( pose.torsion( torsion_id ) ) );
		TR.Debug << "Psi: " << psi << endl;

		// Third, we need to deal with L-sugars.
		CarbohydrateInfoCOP info( rsd.carbohydrate_info() );
		uint const next_rsd( get_downstream_residue_that_this_torsion_moves( pose, torsion_id ) );
		CarbohydrateInfoCOP next_info( pose.residue( next_rsd ).carbohydrate_info() );
		bool const is_exocyclic_bond( has_exocyclic_glycosidic_linkage( pose, next_rsd ) );

		// Psis TO L-Sugars use the mirror image of the score functions, if NOT exocyclic.
		// Psis FROM L-Sugars use the mirror image of the score functions, IF exocyclic
		if ( ! is_exocyclic_bond ) {
			if ( info->is_L_sugar() ) {
				psi = 360 - psi;
			}
		} else {
			if ( next_info->is_L_sugar() ) {
				psi = 360 - psi;
			}
		}

		// Finally, evaluate.
		deriv = E_cef_.evaluate_derivative(
				get_CHI_energy_function_linkage_type_for_psi_for_residue_in_pose( pose, next_rsd ), psi );

		// ...and deal with L-sugars again.
		if ( ! is_exocyclic_bond ) {
			if ( info->is_L_sugar() ) {
				deriv = -deriv;
			}
		} else {
			if ( next_info->is_L_sugar() ) {
				deriv = -deriv;
			}
		}
	} else if ( is_glycosidic_omega_torsion( pose, torsion_id ) ) {
		// Rosetta defines this torsion angle the same way as IUPAC; however, it is the omega of the following residue.
		// This should not matter, though, because the omega preference function only cares about whether the O4 of the
		// reducing end sugar (for aldohexopyranoses is axial or equatorial.  We can simply convert to a number in
		// degrees between 0 and 360 and call the function.

		// First, make sure this is a sugar.
		if ( ! rsd.is_carbohydrate() ) { return deriv; }

		CarbohydrateInfoCOP info( rsd.carbohydrate_info() );

		// Second, get the omega and convert it to between 0 and 360 (because that's what the function expects).
		Angle omega( nonnegative_principal_angle_degrees( pose.torsion( torsion_id ) ) );
		TR.Debug << "Omega: " << omega << endl;

		if ( info->is_L_sugar() ) {
			// L-Sugars use the mirror image of the score functions.
			omega = 360 - omega;
		}

		// Third, what is the next residue?
		uint const next_rsd( get_downstream_residue_that_this_torsion_moves( pose, torsion_id ) );

		deriv = E_opf_.evaluate_derivative( get_omega_preference_for_residue_in_pose( pose, next_rsd ), omega );

		if ( info->is_L_sugar() ) {
			deriv = -deriv;
		}
	}
	return weights[ sugar_bb ] * deriv * 180/pi;  // Convert back into radians.
}


// Creator methods ////////////////////////////////////////////////////////////
// Return an up-casted owning pointer (EnergyMethodOP) to the energy method.
EnergyMethodOP
SugarBackboneEnergyCreator::create_energy_method( EnergyMethodOptions const & ) const
{
	return EnergyMethodOP( new SugarBackboneEnergy );
}

// Return the set of ScoreTypes for which this EnergyMethod is responsible.
ScoreTypes
SugarBackboneEnergyCreator::score_types_for_method() const
{
	ScoreTypes types;
	types.push_back( sugar_bb );
	return types;
}

}  // namespace carbohydrates
}  // namespace methods
}  // namespace scoring
}  // namespace core
