// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.cc
/// @brief  Atom pair energy functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)


// Unit headers
#include <core/scoring/methods/ChainbreakEnergy.hh>
#include <core/scoring/methods/ChainbreakEnergyCreator.hh>

// Package headers


// Project headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/EnergyMap.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>

//Auto Headers
#include <core/id/AtomID.hh>


// Numeric headers


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the ChainbreakEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ChainbreakEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new ChainbreakEnergy;
}

ScoreTypes
ChainbreakEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( chainbreak );
	return sts;
}


ChainbreakEnergy::ChainbreakEnergy() :
	parent( new ChainbreakEnergyCreator )
{}

/// called at the end of energy evaluation
/// In this case (ChainbreakEnergy), all the calculation is done here

void
ChainbreakEnergy::finalize_total_energy(
 pose::Pose & pose,
 ScoreFunction const &,
 EnergyMap & totals
) const
{
	using conformation::Residue;
	Real total_dev(0.0);
	for ( int n=1; n<= pose.fold_tree().num_cutpoint(); ++n ) {
		int const cutpoint( pose.fold_tree().cutpoint( n ) );
		Residue const & lower_rsd( pose.residue( cutpoint ) );
		if ( !lower_rsd.has_variant_type( chemical::CUTPOINT_LOWER ) ) continue;

		Residue const & upper_rsd( pose.residue( cutpoint+1 ) );
		assert( upper_rsd.has_variant_type( chemical::CUTPOINT_UPPER ) );
		Size const nbb( lower_rsd.mainchain_atoms().size() );
		total_dev +=
			( upper_rsd.atom( upper_rsd.mainchain_atoms()[  1] ).xyz().distance_squared( lower_rsd.atom( "OVL1" ).xyz() ) +
				upper_rsd.atom( upper_rsd.mainchain_atoms()[  2] ).xyz().distance_squared( lower_rsd.atom( "OVL2" ).xyz() ) +
				lower_rsd.atom( lower_rsd.mainchain_atoms()[nbb] ).xyz().distance_squared( upper_rsd.atom( "OVU1" ).xyz() ) );
	}
	assert( std::abs( totals[ chainbreak ] ) < 1e-3 );
	totals[ chainbreak ] = total_dev;
}



	/// called during gradient-based minimization inside dfunc
	/**
		 F1 and F2 are not zeroed -- contributions from this atom are
		 just summed in
	**/

void
ChainbreakEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const &, // sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	using conformation::Residue;
	using chemical::CUTPOINT_LOWER;
	using chemical::CUTPOINT_UPPER;

	if ( pose.fold_tree().is_cutpoint( id.rsd() ) && pose.residue(id.rsd() ).has_variant_type( CUTPOINT_LOWER ) ) {
		Residue const & lower_rsd( pose.residue( id.rsd()     ) );
		Residue const & upper_rsd( pose.residue( id.rsd() + 1 ) );
		Vector const & xyz_moving( pose.xyz( id ) );
		bool match( true );
		Vector xyz_fixed;
		Size const nbb( lower_rsd.mainchain_atoms().size() );
		if ( id.atomno() == lower_rsd.mainchain_atoms()[nbb] ) {
			xyz_fixed = upper_rsd.atom( "OVU1" ).xyz();
		} else if ( lower_rsd.atom_name( id.atomno() ) == "OVL1" ) {
			xyz_fixed = upper_rsd.atom( upper_rsd.mainchain_atoms()[1] ).xyz();
		} else if ( lower_rsd.atom_name( id.atomno() ) == "OVL2" ) {
			xyz_fixed = upper_rsd.atom( upper_rsd.mainchain_atoms()[2] ).xyz();
		} else {
			match = false;
		}

		if ( match ) {
			// deriv = 2 * r
			// factor = deriv / r = 2
			F1 += weights[ chainbreak ] * 2 * cross( xyz_moving, xyz_fixed );
			F2 += weights[ chainbreak ] * 2 * ( xyz_moving - xyz_fixed );
		}
	}

	if ( id.rsd() > 1 && pose.fold_tree().is_cutpoint( id.rsd()-1 ) &&
			 pose.residue(id.rsd() ).has_variant_type( CUTPOINT_UPPER )  ) {
		Residue const & lower_rsd( pose.residue( id.rsd() - 1 ) );
		Residue const & upper_rsd( pose.residue( id.rsd()     ) );
		Vector const & xyz_moving( pose.xyz( id ) );
		bool match( true );
		Vector xyz_fixed;
		Size const nbb( lower_rsd.mainchain_atoms().size() );
		if ( id.atomno() == upper_rsd.mainchain_atoms()[1] ) {
			xyz_fixed = lower_rsd.atom( "OVL1" ).xyz();
		} else if ( id.atomno() == upper_rsd.mainchain_atoms()[2] ) {
			xyz_fixed = lower_rsd.atom( "OVL2" ).xyz();
		} else if ( upper_rsd.atom_name( id.atomno() ) == "OVU1" ) {
			xyz_fixed = lower_rsd.atom( lower_rsd.mainchain_atoms()[nbb] ).xyz();
		} else {
			match = false;
		}

		if ( match ) {
			// deriv = 2 * r
			// factor = deriv / r = 2
			F1 += weights[ chainbreak ] * 2 * cross( xyz_moving, xyz_fixed );
			F2 += weights[ chainbreak ] * 2 * ( xyz_moving - xyz_fixed );
		}
	}
}

/// @brief Chainbreak Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
ChainbreakEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const
{}
core::Size
ChainbreakEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace methods
} // namespace scoring
} // namespace core
