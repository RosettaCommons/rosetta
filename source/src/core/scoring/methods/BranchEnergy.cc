// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/methods/BranchEnergy.cc
/// @brief  Method definitions for BranchEnergyCreator and BranchEnergy classes
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)


// Unit headers
#include <core/scoring/methods/branch_energy_util.hh>
#include <core/scoring/methods/BranchEnergy.hh>
#include <core/scoring/methods/BranchEnergyCreator.hh>

// Project headers
#include <core/id/AtomID.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/EnergyMap.hh>

#include <basic/Tracer.hh>

// Utility header
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

static basic::Tracer tr( "core.scoring.methods.Branch", basic::t_info );

/// @details This must return a fresh instance of the BranchEnergy class, never an instance already in use.
methods::EnergyMethodOP
BranchEnergyCreator::create_energy_method( methods::EnergyMethodOptions const & ) const
{
	return methods::EnergyMethodOP( new BranchEnergy );
}

ScoreTypes
BranchEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( branch_conn );
	return sts;
}


BranchEnergy::BranchEnergy() : parent( methods::EnergyMethodCreatorOP( new BranchEnergyCreator ) )
{}

// Called at the end of the energy evaluation.
/// @details In this case (BranchEnergy), all the calculation is done here.
void
BranchEnergy::finalize_total_energy( pose::Pose & pose, ScoreFunction const &, EnergyMap & totals ) const
{
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		using namespace core::conformation::symmetry;
		auto const & symm_conf(
			dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
	}

	using conformation::Residue;
	using namespace core::chemical;
	DistanceSquared total_dev( 0.0 );
	//utility::vector1< int > cutpoints;
	utility::vector1< ResidueAtomOverlaps > branch_connections;
	// This grabs all the PAIRS of residues that need scoring by this term.
	find_relevant_connections( pose, branch_connections );

	for ( auto const & elem : branch_connections ) {

		// Hm, okay, which atoms to compare?
		DistanceSquared current_dev( 0.0 );

		auto const & lower_rsd = pose.residue( elem.res1 );
		auto const & upper_rsd = pose.residue( elem.res2 );

		current_dev += upper_rsd.atom( elem.res1_ovl1_overlaps ).xyz().distance_squared(
			lower_rsd.atom( "OVL1" ).xyz() );

		current_dev += lower_rsd.atom( elem.res2_ovu1_overlaps ).xyz().distance_squared(
			upper_rsd.atom( "OVU1" ).xyz() );

		current_dev += upper_rsd.atom( elem.res1_ovl2_overlaps ).xyz().distance_squared(
			lower_rsd.atom( "OVL2" ).xyz() );

		total_dev += current_dev;
	}
	debug_assert( std::abs( totals[ branch_conn ] ) < 1e-3 );
	totals[ branch_conn ] = total_dev;
}


// Called during gradient-based minimization inside dfunc.
/// @note F1 and F2 are not zeroed -- contributions from this atom are just summed in.
void
BranchEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const &, // sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	// Like chainbreak, we assume EITHER WAY could happen, so we want to get the
	// direction info from a special version of the prior function

	using conformation::Residue;
	using core::chemical::CUTPOINT_LOWER;
	using core::chemical::CUTPOINT_UPPER;

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		using namespace core::conformation::symmetry;

		auto const & symm_conf(
			dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
		Size max_res = symm_info->num_independent_residues() - 1;
		if ( id.rsd() > max_res ) return;
	}

	ResidueAtomOverlaps branch_connection;
	find_relevant_connections_onersd( pose, id.rsd(), branch_connection );

	Vector const & xyz_moving( pose.xyz( id ) );  // position of the moving atom
	Vector xyz_fixed;  // where the moving atom should superimpose if the chain is no longer broken
	bool match( false );


	// We should catch problems below but somehow it doesn't? explicit here.
	if ( branch_connection.res1 == 0 ) return;

	// There are only 6 cases we care about; check for a match with one of them.

	// If the moving atom is on the upstream (lower) side of the cutpoint,...
	if ( branch_connection.res1 == id.rsd() ) {
		// Get the two residues across the cutpoint.
		Residue const & lower_rsd( pose.residue( id.rsd() ) );
		Residue const & upper_rsd( pose.residue( branch_connection.res2 ) );

		// Case 1: The moving atom is the last main-chain atom on the upstream residue across the cutpoint and
		// should superimpose with virtual atom OVU1 from the downstream residue.
		if ( id.atomno() == lower_rsd.atom_index( branch_connection.res2_ovu1_overlaps ) ) {
			xyz_fixed = upper_rsd.atom( "OVU1" ).xyz();
			match = true;

			// Case 2: The moving atom is the virtual atom OVL1 on the upstream residue across the cutpoint and
			// should superimpose with the 1st main-chain atom from the downstream residue.
		} else if ( lower_rsd.atom_name( id.atomno() ) == "OVL1" ) {
			xyz_fixed = upper_rsd.atom( branch_connection.res1_ovl1_overlaps ).xyz();
			match = true;

			// Case 3: The moving atom is the virtual atom OVL2 on the upstream residue across the cutpoint and
			// should superimpose with the 2nd main-chain atom from the downstream residue.
		} else if ( lower_rsd.atom_name( id.atomno() ) == "OVL2" ) {
			xyz_fixed = upper_rsd.atom( branch_connection.res1_ovl2_overlaps ).xyz();
			match = true;
		}
	}

	// If the moving atom is on the downstream (upper) side of the cutpoint,...
	// (It is possible that a residue can simultaneously be a CUTPOINT_LOWER and _UPPER variant; hence, no else if
	// here.  However, this code assumes that no atom on such a residue can simultaneously require superimposition with
	// two targets.  Such a case would require a residue with only two main-chain atoms with a partial-double bond con-
	// nection or one with only one main-chain atom and a single bond connection.  Neither case is at all likely.
	// ~Labonte)
	if ( ! match && branch_connection.res2 == id.rsd() ) {
		// Get the two residues across the cutpoint.
		Residue const & lower_rsd( pose.residue( branch_connection.res1 ) );
		Residue const & upper_rsd( pose.residue( id.rsd() ) );

		// Case 4: The moving atom is the 1st main-chain atom on the downstream residue across the cutpoint and
		// should superimpose with virtual atom OVL1 from the upstream residue.
		if ( id.atomno() == upper_rsd.atom_index( branch_connection.res1_ovl1_overlaps ) ) {
			xyz_fixed = lower_rsd.atom( "OVL1" ).xyz();
			match = true;

			// Case 5: The moving atom is the 2nd main-chain atom on the downstream residue across the cutpoint and
			// should superimpose with virtual atom OVL2 from the upstream residue, IF the residue is not a sugar.
		} else if ( id.atomno() == upper_rsd.atom_index( branch_connection.res1_ovl2_overlaps ) ) {
			xyz_fixed = lower_rsd.atom( "OVL2" ).xyz();
			match = true;

			// Case 6: The moving atom is the virtual atom OVU1 on the downstream residue across the cutpoint and
			// should superimpose with the last main-chain atom from the upstream residue.
		} else if ( upper_rsd.atom_name( id.atomno() ) == "OVU1" ) {
			xyz_fixed = lower_rsd.atom( lower_rsd.atom_index( branch_connection.res2_ovu1_overlaps ) ).xyz();
			match = true;
		}
	}

	if ( match ) {
		F1 += weights[ branch_conn ] * 2 * cross( xyz_moving, xyz_fixed );
		F2 += weights[ branch_conn ] * 2 * ( xyz_moving - xyz_fixed );
	}
}


/// @note BranchEnergy is context-independent and thus indicates that no context graphs need to be maintained by
/// class Energies.
void
BranchEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{}


core::Size
BranchEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace methods
} // namespace scoring
} // namespace core
