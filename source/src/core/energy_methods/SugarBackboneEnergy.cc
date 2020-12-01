// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/energy_methods/SugarBackboneEnergy.hh
/// @brief   Method definitions for SugarBackboneEnergy and SugarBackboneEnergyCreator.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit Headers
#include <core/energy_methods/SugarBackboneEnergy.hh>
#include <core/energy_methods/SugarBackboneEnergyCreator.hh>

// Package Headers
#include <core/scoring/carbohydrates/CHIEnergyFunction.hh>
#include <core/scoring/carbohydrates/OmegaPreferencesFunction.hh>
#include <core/scoring/carbohydrates/util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>

// Project Headers
#include <core/id/AtomID.hh>
#include <core/id/PartialAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>

// Numeric Headers
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/angle.functions.hh>

// Basic Header
#include <basic/Tracer.hh>


// Construct tracer.
static basic::Tracer TR( "core.scoring.methods.carbohydrates.SugarBackboneEnergy" );


namespace core {
namespace energy_methods {

using namespace core::conformation::carbohydrates;

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
SugarBackboneEnergy::SugarBackboneEnergy() :
	core::scoring::methods::ContextIndependentOneBodyEnergy( utility::pointer::make_shared< SugarBackboneEnergyCreator >() ),
	E_cef_( core::scoring::ScoringManager::get_instance()->get_CHIEnergyFunction() ),
	E_opf_( core::scoring::ScoringManager::get_instance()->get_OmegaPreferencesFunction() )
{}


// General EnergyMethod Methods ///////////////////////////////////////////////
core::scoring::methods::EnergyMethodOP
SugarBackboneEnergy::clone() const
{
	return utility::pointer::make_shared< SugarBackboneEnergy >();
}


// OneBodyEnergy Methods //////////////////////////////////////////////////////
/// @details Evaluate the one-body carbohydrate backbone energies for a particular residue.
/// @note This is not a one body energy. It reaches into the pose to get the coordinates of
/// a second residue. This is a two body energy.
void
SugarBackboneEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	core::scoring::EnergyMap & emap ) const
{
	using namespace numeric;
	using namespace chemical::carbohydrates;
	using namespace pose::carbohydrates;

	// This is a carbohydrate-only scoring method.
	if ( ! rsd.is_carbohydrate() ) { return; }

	// Phi, psi, and omega are meaningless for reducing-end sugars.
	if ( rsd.is_lower_terminus() ) { return; }

	// Ignore REPLONLY variants.
	if ( rsd.has_variant_type( chemical::REPLONLY ) ) { return; }

	// Ignore VIRTUAL residues.
	if ( rsd.is_virtual_residue() ) { return; }

	// Get parent residue, CarbohydrateInfo, and exocyclic state for convenience.
	CarbohydrateInfoCOP info( rsd.carbohydrate_info() );

	bool is_exocyclic_bond;
	core::Size prev_rsd_num;

	if ( pose.glycan_tree_set() ) {
		prev_rsd_num =  pose.glycan_tree_set()->get_parent( rsd.seqpos() );
		is_exocyclic_bond = pose.glycan_tree_set()->has_exocyclic_glycosidic_linkage( rsd.seqpos() );
	} else {
		prev_rsd_num =  find_seqpos_of_saccharides_parent_residue( rsd );
		is_exocyclic_bond = has_exocyclic_glycosidic_linkage( pose.conformation(), rsd.seqpos() );
	}

	if ( prev_rsd_num == 0 ) return;

	conformation::Residue const & prev_rsd( pose.residue( prev_rsd_num));
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
	score += E_cef_( core::scoring::carbohydrates::get_CHI_energy_function_linkage_type_for_phi_for_residue_in_pose( pose, seqpos ), phi );


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
	score += E_cef_( core::scoring::carbohydrates::get_CHI_energy_function_linkage_type_for_psi_for_residue_in_pose( pose, seqpos ), psi );


	// Calculate omega component. /////////////////////////////////////////////
	// Omegas TO L-Sugars use the mirror image of the score functions.
	if ( prev_rsd.is_carbohydrate() ) {
		if ( ! is_exocyclic_bond ) {
			if ( prev_rsd.carbohydrate_info()->is_L_sugar() ) {
				omega = 360 - omega;
			}
		}
	}

	// Score is 0 if linkage type is PREFERENCE_NA.
	score += E_opf_( core::scoring::carbohydrates::get_omega_preference_for_residue_in_pose( pose, seqpos ), omega );


	emap[ core::scoring::sugar_bb ] += score;
}

utility::vector1< id::PartialAtomID >
SugarBackboneEnergy::atoms_with_dof_derivatives(
	conformation::Residue const & rsd,
	pose::Pose const & pose
) const
{
	utility::vector1< id::PartialAtomID > atoms;

	// TO DO!
	// The dihedral angles for a particular residue are
	// defined mostly by the coordinates of the atom in
	// the parent residue, therefore this term should
	// be defined as a two-body energy.

	// This is a carbohydrate-only scoring method.
	if ( ! rsd.is_carbohydrate() ) { return atoms; }

	// Phi, psi, and omega are meaningless for reducing-end sugars.
	if ( rsd.is_lower_terminus() ) { return atoms; }

	// Ignore REPLONLY variants.
	if ( rsd.has_variant_type( chemical::REPLONLY ) ) { return atoms; }

	// Ignore VIRTUAL residues.
	if ( rsd.is_virtual_residue() ) { return atoms; }

	// Get parent residue, CarbohydrateInfo, and exocyclic state for convenience.
	chemical::carbohydrates::CarbohydrateInfoCOP info( rsd.carbohydrate_info() );

	core::Size prev_rsd_num;

	if ( pose.glycan_tree_set() ) {
		prev_rsd_num =  pose.glycan_tree_set()->get_parent( rsd.seqpos() );
	} else {
		prev_rsd_num =  find_seqpos_of_saccharides_parent_residue( rsd );
	}

	if ( prev_rsd_num == 0 ) return atoms;

	std::pair< conformation::ResidueCOP, conformation::ResidueCOP > rsd_and_parent =
		conformation::carbohydrates::get_glycosidic_bond_residues(
		pose.conformation(), rsd.seqpos() );
	if ( rsd_and_parent.first->seqpos() == rsd_and_parent.second->seqpos() ) {
		return atoms;
	}

	utility::vector1< id::AtomID > phi_ats = conformation::carbohydrates::get_reference_atoms_for_phi(
		pose.conformation(), rsd.seqpos() );
	//std::set< id::PartialAtomID > atoms;
	atoms.reserve(5);
	atoms.push_back( id::PartialAtomID( phi_ats[1].atomno(), phi_ats[1].rsd() ));
	atoms.push_back( id::PartialAtomID( phi_ats[2].atomno(), phi_ats[2].rsd() ));

	// Find how the parent connects to res
	auto find_conn = [&]() {
		for ( Size ii = 1; ii <= rsd.connect_map_size(); ++ii ) {
			if ( rsd.connect_map(ii).resid() == rsd_and_parent.second->seqpos() ) {
				return rsd.connect_map(ii).connid();
			}
		}
		return Size(0);
	};
	Size const conn_to_phi_res = find_conn();
	Size const parent_res = rsd_and_parent.second->seqpos();
	if ( conn_to_phi_res != 0 ) {
		atoms.push_back( id::PartialAtomID( conn_to_phi_res, parent_res, 0 )); // phi, psi, omega
		atoms.push_back( id::PartialAtomID( conn_to_phi_res, parent_res, 1 )); // phi, psi, omega
		atoms.push_back( id::PartialAtomID( conn_to_phi_res, parent_res, 2 )); // psi, omega
		atoms.push_back( id::PartialAtomID( conn_to_phi_res, parent_res, 3 )); // omega
	}
	return atoms;
}


// Evaluate the DoF derivative for a particular residue.
core::Real
SugarBackboneEnergy::eval_residue_dof_derivative(
	conformation::Residue const & rsd,
	core::scoring::ResSingleMinimizationData const & /* min_data */,
	id::DOF_ID const & /* dof_id */,
	id::TorsionID const & torsion_id,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const & /* sf */,
	core::scoring::EnergyMap const & weights ) const
{
	using namespace std;
	using namespace numeric;
	using namespace numeric::constants::d;
	using namespace id;
	using namespace chemical::carbohydrates;
	using namespace pose::carbohydrates;

	TR.Debug << "Evaluating torsion: " << torsion_id << endl;

	Real deriv( 0.0 );

	// Ignore REPLONLY variants.
	if ( rsd.has_variant_type( chemical::REPLONLY ) ) { return deriv; }

	// Ignore VIRTUAL residues.
	if ( rsd.is_virtual_residue() ) { return deriv; }

	// Ignore LIGAND residues, as they cannot have a sugar backbone energy.
	if ( rsd.is_ligand() ) { return deriv; }

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

		// If the connecting residue defining phi is virtual, phi scoring is undefined.
		if ( pose.residue(next_rsd).is_virtual_residue() ) return deriv;
		if ( next_rsd == 0 ) return deriv;

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
			core::scoring::carbohydrates::get_CHI_energy_function_linkage_type_for_phi_for_residue_in_pose( pose, next_rsd ), phi );

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

		// If the connecting residue defining psi/omega is virtual, scoring is undefined.
		if ( pose.residue(next_rsd).is_virtual_residue() ) return deriv;
		if ( next_rsd == 0 ) return deriv;

		CarbohydrateInfoCOP next_info( pose.residue( next_rsd ).carbohydrate_info() );

		bool const is_exocyclic_bond( pose.glycan_tree_set()->has_exocyclic_glycosidic_linkage( next_rsd ));

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
			core::scoring::carbohydrates::get_CHI_energy_function_linkage_type_for_psi_for_residue_in_pose( pose, next_rsd ), psi );

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
		if ( next_rsd == 0 ) return deriv;

		// If the connecting residue defining phi is virtual, omega scoring is undefined.
		if ( pose.residue(next_rsd).is_virtual_residue() ) return deriv;

		deriv = E_opf_.evaluate_derivative( core::scoring::carbohydrates::get_omega_preference_for_residue_in_pose( pose, next_rsd ), omega );

		if ( info->is_L_sugar() ) {
			deriv = -deriv;
		}
	}
	return weights[ core::scoring::sugar_bb ] * deriv * 180/pi;  // Convert back into radians.
}


// Creator methods ////////////////////////////////////////////////////////////
// Return an up-casted owning pointer (core::scoring::methods::EnergyMethodOP) to the energy method.
core::scoring::methods::EnergyMethodOP
SugarBackboneEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const & ) const
{
	return utility::pointer::make_shared< SugarBackboneEnergy >();
}

// Return the set of ScoreTypes for which this EnergyMethod is responsible.
core::scoring::ScoreTypes
SugarBackboneEnergyCreator::score_types_for_method() const
{
	using namespace core::scoring;
	ScoreTypes types;
	types.push_back( sugar_bb );
	return types;
}

}  // namespace energy_methods
}  // namespace core
