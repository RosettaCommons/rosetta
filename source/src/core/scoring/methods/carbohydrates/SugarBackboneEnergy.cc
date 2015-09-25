// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/methods/carbohydrates/SugarBackboneEnergy.hh
/// @brief   Method definitions for SugarBackboneEnergy and SugarBackboneEnergyCreator.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit headers
#include <core/scoring/methods/carbohydrates/SugarBackboneEnergy.hh>
#include <core/scoring/methods/carbohydrates/SugarBackboneEnergyCreator.hh>
#include <core/scoring/carbohydrates/CHIEnergyFunction.hh>
#include <core/scoring/carbohydrates/LinkageType.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/angle.functions.hh>

// Basic header
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
	E_( ScoringManager::get_instance()->get_CHIEnergyFunction() )
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
	using namespace scoring::carbohydrates;

	// This is a carbohydrate-only scoring method.
	if ( ! rsd.is_carbohydrate() ) { return; }

	// Ignore REPLONLY variants.
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) { return; }

	uint const seqpos( rsd.seqpos() );
	Angle phi( pose.phi( seqpos ) );
	//Angle psi( pose.psi( seqpos ) );
	chemical::carbohydrates::CarbohydrateInfoCOP info( rsd.carbohydrate_info() );

	Energy score( 0.0 );

	// Calculate phi component.
	if ( info->is_L_sugar() ) {
		// L-Sugars use the mirror image of the score functions.
		phi = 360 - phi;
	}
	if ( info->is_alpha_sugar() ) {
		score += E_( ALPHA_LINKS, phi );
	} else if ( info->is_beta_sugar() ) {
		score += E_( BETA_LINKS, phi );
	}  // ...else it's a linear sugar, and this scoring method does not apply.

	// TODO: Calculate psi component.

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
	using namespace id;
	using namespace chemical::rings;
	using namespace chemical::carbohydrates;
	using namespace pose::carbohydrates;
	using namespace scoring::carbohydrates;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Evaluating torsion: " << torsion_id << endl;
	}

	Real deriv( 0.0 );

	// This is a carbohydrate-only scoring method.
	if ( ! rsd.is_carbohydrate() ) { return deriv; }

	// Ignore REPLONLY variants.
	if ( rsd.has_variant_type( chemical::REPLONLY ) ) { return deriv; }

	// This scoring method only considers glycosidic torsions, which may have either BB, CHI, or BRANCH TorsionIDs.
	if ( ! torsion_id.valid() ) {
		if ( TR.Debug.visible() ) {
			TR.Debug << "Torsion not valid: " << torsion_id << endl;
		}
		return deriv;
	}

	if ( is_phi_torsion( pose, torsion_id ) ) {
		// Rosetta defines this torsion angle differently than IUPAC.  In addition, this torsion belongs to the next
		// residue.  We need to figure out what that residue is, use pose.phi() to get the value, and see if it's
		// alpha or beta.

		// First, what is the next residue?
		uint next_rsd( 0 );
		if ( torsion_id.type() == BB ) {
			// If this is a main-chain torsion, we can be confident that the next residue on this chain MUST be n+1.
			next_rsd = torsion_id.rsd() + 1;
		} else if ( torsion_id.type() == BRANCH ) {
			Size const n_mainchain_connections( rsd.n_polymeric_residue_connections() );
			next_rsd = rsd.residue_connection_partner( n_mainchain_connections + torsion_id.torsion() );
		} else {
			TR.Error << "Torsion " << torsion_id << " cannot be a phi torsion!" << endl;
			return deriv;
		}

		// Now, get the next residue's phi and its info.
		CarbohydrateInfoCOP info( pose.residue( next_rsd ).carbohydrate_info() );
		Angle phi( pose.phi( next_rsd ) );
		if ( TR.Debug.visible() ) {
			TR.Debug << "Phi: " << phi << endl;
		}
		if ( info->is_L_sugar() ) {
			// L-Sugars use the mirror image of the score functions. The phi functions run from -180 to 180.
			//phi = principal_angle_degrees( 360 - phi );
			phi = -phi;
		}

		// Finally, we can evaluate.
		if ( info->is_alpha_sugar() ) {
			deriv = E_.evaluate_derivative( ALPHA_LINKS, phi );
		} else if ( info->is_beta_sugar() ) {
			deriv = E_.evaluate_derivative( BETA_LINKS, phi );
		}  // ...else it's a linear sugar, and this scoring method does not apply.
	} else if ( is_psi_torsion( pose, torsion_id ) ) {
		// Rosetta defines this torsion angle the same way as IUPAC; however, it is the psi of the following residue.
		// This should not matter, though, because the CHI energy function only cares about whether the linkage is
		// axial or equatorial.  We can simply convert to a number in degrees between 0 and 360 and call the function.

		// First, if this is not a pyranose, do nothing, because we do not have statistics for this residue.
		CarbohydrateInfoCOP info( rsd.carbohydrate_info() );
		if ( ! info->is_pyranose() ) {
			return deriv;
		}

		// Second, what is our connecting atom?
		core::uint connect_atom;
		if ( torsion_id.type() == BB ) {
			// If this is a main-chain torsion, the connect atom is the UPPER_CONNECT.
			connect_atom = rsd.upper_connect_atom();
		} else if ( torsion_id.type() == CHI ) {
			// If this is a side-chain torsion, the CONNECT atom will be the 3rd atom of the CHI definition.
			connect_atom = rsd.chi_atoms()[ torsion_id.torsion() ][ 3 ];
		} else {
			TR.Error << "Torsion " << torsion_id << " cannot be a psi torsion!" << endl;
			return deriv;
		}

		// TODO: If I ever add a score for exocyclic psis, I should check if this atom is exocyclic and act
		// accordingly.  For now, the code will just return 0.0, because an exocyclic atom will be NEITHER. ~Labonte

		// Next, get the psi and convert it to between 0 and 360 (because that's what the function expects).
		Angle psi( nonnegative_principal_angle_degrees( pose.torsion( torsion_id ) ) );
		if ( TR.Debug.visible() ) {
			TR.Debug << "Psi: " << psi << endl;
		}
		if ( info->is_L_sugar() ) {
			// L-Sugars use the mirror image of the score functions.
			psi = -psi;
		}

		// Now, get the ring atoms. We can assume that the ring we care about is the 1st ring, since this is a sugar.
		utility::vector1< core::uint > const ring_atoms( rsd.type().ring_atoms( 1 ) );

		// Finally, check if it's axial or equatorial and call the appropriate function.
		switch ( is_atom_axial_or_equatorial_to_ring( rsd, connect_atom, ring_atoms ) ) {
		case AXIAL :
			if ( torsion_id.torsion() % 2 == 0 ) {  // even
				deriv = E_.evaluate_derivative( _2AX_3EQ_4AX_LINKS, psi );
			} else /* odd */ {
				deriv = E_.evaluate_derivative( _2EQ_3AX_4EQ_LINKS, psi );
			}
			break;
		case EQUATORIAL :
			if ( torsion_id.torsion() % 2 == 0 ) {  // even
				deriv = E_.evaluate_derivative( _2EQ_3AX_4EQ_LINKS, psi );
			} else /* odd */ {
				deriv = E_.evaluate_derivative( _2AX_3EQ_4AX_LINKS, psi );
			}
			break;
		case NEITHER :
			break;
		}
	}
	return weights[ sugar_bb ] * deriv;
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
