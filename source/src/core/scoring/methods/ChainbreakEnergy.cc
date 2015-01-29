// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/methods/ChainbreakEnergy.cc
/// @brief  Method definitions for ChainbreakEnergyCreator and ChainbreakEnergy classes
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)


// Unit headers
#include <core/scoring/methods/ChainbreakEnergy.hh>
#include <core/scoring/methods/ChainbreakEnergyCreator.hh>

// Project headers
#include <core/id/AtomID.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/EnergyMap.hh>

// Utility header
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

/// @details This must return a fresh instance of the ChainbreakEnergy class, never an instance already in use.
methods::EnergyMethodOP
ChainbreakEnergyCreator::create_energy_method( methods::EnergyMethodOptions const & ) const
{
	return methods::EnergyMethodOP( new ChainbreakEnergy );
}

ScoreTypes
ChainbreakEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( chainbreak );
	return sts;
}


ChainbreakEnergy::ChainbreakEnergy() : parent( methods::EnergyMethodCreatorOP( new ChainbreakEnergyCreator ) )
{}


// Called at the end of the energy evaluation.
/// @details In this case (ChainbreakEnergy), all the calculation is done here.
void
ChainbreakEnergy::finalize_total_energy( pose::Pose & pose, ScoreFunction const &, EnergyMap & totals ) const
{
	Size max_res = pose.n_residue();
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		using namespace core::conformation::symmetry;
		SymmetricConformation const & symm_conf(
				dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
		max_res = symm_info->num_independent_residues() - 1;
	}

	using conformation::Residue;
	using namespace core::chemical;
	DistanceSquared total_dev( 0.0 );
	for ( int n = 1; n <= pose.fold_tree().num_cutpoint(); ++n ) {
		int const cutpoint( pose.fold_tree().cutpoint( n ) );
		if ( cutpoint > static_cast< int >( max_res ) ) continue;
		Residue const & lower_rsd( pose.residue( cutpoint ) );
		if ( ! lower_rsd.has_variant_type( CUTPOINT_LOWER ) ) continue;

		// This logic restricts cutpoints from occurring at branch points, but that may be ok. ~Labonte
		Residue const & upper_rsd( pose.residue( cutpoint + 1 ) );
	debug_assert( upper_rsd.has_variant_type( CUTPOINT_UPPER ) );
		Size const last_mainchain_atm( lower_rsd.mainchain_atoms().size() );

		DistanceSquared current_dev( 0.0 );

		// How well does the 1st main-chain atom of the downstream (upper) residue of this cutpoint overlap
		// with the corresponding virtual atom stemming from the upstream (lower) residue of this cutpoint?
		current_dev += upper_rsd.atom( upper_rsd.mainchain_atoms()[ 1 ] ).xyz().distance_squared(
				lower_rsd.atom( "OVL1" ).xyz() );

		// How well does the last main-chain atom of the upstream (upper) residue of this cutpoint overlap
		// with the corresponding virtual atom stemming from the downstream (lower) residue of this cutpoint?
		current_dev += lower_rsd.atom( lower_rsd.mainchain_atoms()[ last_mainchain_atm ] ).xyz().distance_squared(
				upper_rsd.atom( "OVU1" ).xyz() );

		// If this is a cutpoint on a carbohydrate chain, we don't care at all about a third virtual atom, because the
		// residues connect at a single bond, not a partial-double as is the case with a peptide bond.
		if ( upper_rsd.is_carbohydrate() ) {
			current_dev *= 1.5;  // Normalize the score with that of non-carbohydrate chain breaks.
		} else {
			// How well does the 2nd main-chain atom of the downstream (upper) residue of this cutpoint overlap
			// with the corresponding virtual atom stemming from the upstream (lower) residue of this cutpoint?
			current_dev += upper_rsd.atom( upper_rsd.mainchain_atoms()[ 2 ] ).xyz().distance_squared(
					lower_rsd.atom( "OVL2" ).xyz() );
		}

		total_dev += current_dev;
	}
debug_assert( std::abs( totals[ chainbreak ] ) < 1e-3 );
	totals[ chainbreak ] = total_dev;
}


// Called during gradient-based minimization inside dfunc.
/// @note F1 and F2 are not zeroed -- contributions from this atom are just summed in.
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
	using core::chemical::CUTPOINT_LOWER;
	using core::chemical::CUTPOINT_UPPER;

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		using namespace core::conformation::symmetry;

		Size max_res = pose.n_residue();
		SymmetricConformation const & symm_conf(
				dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
		max_res = symm_info->num_independent_residues() - 1;
		if( id.rsd() > max_res ) return;
	}

	Vector const & xyz_moving( pose.xyz( id ) );  // position of the moving atom
	Vector xyz_fixed;  // where the moving atom should superimpose if the chain is no longer broken
	bool match( false );

	// There are only 6 cases we care about; check for a match with one of them.

	// If the moving atom is on the upstream (lower) side of the cutpoint,...
	if ( pose.fold_tree().is_cutpoint( id.rsd() ) &&
			pose.residue( id.rsd() ).has_variant_type( CUTPOINT_LOWER ) ) {
		// Get the two residues across the cutpoint.
		Residue const & lower_rsd( pose.residue( id.rsd() ) );
		Residue const & upper_rsd( pose.residue( id.rsd() + 1 ) );

		core::uint const last_mainchain_atm( lower_rsd.mainchain_atoms().size() );
		// Case 1: The moving atom is the last main-chain atom on the upstream residue across the cutpoint and
		// should superimpose with virtual atom OVU1 from the downstream residue.
		if ( id.atomno() == lower_rsd.mainchain_atoms()[ last_mainchain_atm ] ) {
			xyz_fixed = upper_rsd.atom( "OVU1" ).xyz();
			match = true;

		// Case 2: The moving atom is the virtual atom OVL1 on the upstream residue across the cutpoint and
		// should superimpose with the 1st main-chain atom from the downstream residue.
		} else if ( lower_rsd.atom_name( id.atomno() ) == "OVL1" ) {
			xyz_fixed = upper_rsd.atom( upper_rsd.mainchain_atoms()[ 1 ] ).xyz();
			match = true;

		// Case 3: The moving atom is the virtual atom OVL2 on the upstream residue across the cutpoint and
		// should superimpose with the 2nd main-chain atom from the downstream residue.
		} else if ( lower_rsd.atom_name( id.atomno() ) == "OVL2" ) {
			xyz_fixed = upper_rsd.atom( upper_rsd.mainchain_atoms()[ 2 ] ).xyz();
			match = true;
		}
	}

	// If the moving atom is on the downstream (upper) side of the cutpoint,...
	// (It is possible that a residue can simultaneously be a CUTPOINT_LOWER and _UPPER variant; hence, no else if
	// here.  However, this code assumes that no atom on such a residue can simultaneously require superimposition with
	// two targets.  Such a case would require a residue with only two main-chain atoms with a partial-double bond con-
	// nection or one with only one main-chain atom and a single bond connection.  Neither case is at all likely.
	// ~Labonte)
	if ( ! match && id.rsd() > 1 && pose.fold_tree().is_cutpoint( id.rsd() - 1 ) &&
			 pose.residue( id.rsd() ).has_variant_type( CUTPOINT_UPPER ) ) {
		// Get the two residues across the cutpoint.
		Residue const & lower_rsd( pose.residue( id.rsd() - 1 ) );
		Residue const & upper_rsd( pose.residue( id.rsd() ) );

		// Case 4: The moving atom is the 1st main-chain atom on the downstream residue across the cutpoint and
		// should superimpose with virtual atom OVL1 from the upstream residue.
		if ( id.atomno() == upper_rsd.mainchain_atoms()[ 1 ] ) {
			xyz_fixed = lower_rsd.atom( "OVL1" ).xyz();
			match = true;

		// Case 5: The moving atom is the 2nd main-chain atom on the downstream residue across the cutpoint and
		// should superimpose with virtual atom OVL2 from the upstream residue.
		} else if ( id.atomno() == upper_rsd.mainchain_atoms()[ 2 ] ) {
			xyz_fixed = lower_rsd.atom( "OVL2" ).xyz();
			match = true;

		// Case 6: The moving atom is the virtual atom OVU1 on the downstream residue across the cutpoint and
		// should superimpose with the last main-chain atom from the upstream residue.
		} else if ( upper_rsd.atom_name( id.atomno() ) == "OVU1" ) {
			core::uint const last_mainchain_atm( lower_rsd.mainchain_atoms().size() );
			xyz_fixed = lower_rsd.atom( lower_rsd.mainchain_atoms()[ last_mainchain_atm ] ).xyz();
			match = true;
		}
	}

	if ( match ) {
		F1 += weights[ chainbreak ] * 2 * cross( xyz_moving, xyz_fixed );
		F2 += weights[ chainbreak ] * 2 * ( xyz_moving - xyz_fixed );
	}
}


/// @note ChainbreakEnergy is context-independent and thus indicates that no context graphs need to be maintained by
/// class Energies.
void
ChainbreakEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{}


core::Size
ChainbreakEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace methods
} // namespace scoring
} // namespace core
