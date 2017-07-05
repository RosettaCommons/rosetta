// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/LinearBranchEnergy.cc
/// @brief  Atom pair energy functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Christopher Miles (cmiles@uw.edu)
/// @author Rhiju Das (a few updates)

#include <core/scoring/methods/branch_energy_util.hh>
#include <core/scoring/methods/LinearBranchEnergy.hh>
#include <core/scoring/methods/LinearBranchEnergyCreator.hh>

#include <basic/options/keys/OptionKeys.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>

#include <basic/Tracer.hh>

#include <algorithm>
#include <iostream>

#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>


namespace core {
namespace scoring {
namespace methods {

using core::Size;
static THREAD_LOCAL basic::Tracer tr( "core.scoring.methods.LinearBranch", basic::t_info );

LinearBranchEnergy::LinearBranchEnergy()
: parent(methods::EnergyMethodCreatorOP( new LinearBranchEnergyCreator )) {
	initialize( std::numeric_limits<core::Size>::max() );
}

LinearBranchEnergy::LinearBranchEnergy(const LinearBranchEnergy& o)
: parent(o) {
	// assignment to non-cache-related variables
	allowable_sequence_sep_ = o.allowable_sequence_sep_;
}

LinearBranchEnergy& LinearBranchEnergy::operator=(const LinearBranchEnergy& o) {

	if ( this == &o ) return *this;

	// base class assignments
	parent::operator=(o);

	// assignment to non-cache-related variables
	allowable_sequence_sep_ = o.allowable_sequence_sep_;

	return *this;
}

LinearBranchEnergy::~LinearBranchEnergy() {}

void LinearBranchEnergy::initialize(Size allowable_sequence_sep) {
	allowable_sequence_sep_ = allowable_sequence_sep;
}

/// called at the end of energy evaluation
/// In this case (LinearBranchEnergy), all the calculation is done here
void LinearBranchEnergy::finalize_total_energy( pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals ) const
{
	using core::Size;
	using conformation::Residue; using utility::vector1;

	Real total_dev = 0.0;

	// AMW: can't use this function for branches
	// Identify all cutpoint variants defined by the caller

	utility::vector1< ResidueAtomOverlaps > branch_connections;
	find_relevant_connections( pose, branch_connections );

	if ( branch_connections.size() == 0 ) return;

	for ( auto const & elem : branch_connections ) {
		Residue const & lower_rsd = pose.residue( elem.res1 );
		Residue const & upper_rsd = pose.residue( elem.res2 );

		Real score_dev = 0;
		score_dev += lower_rsd.atom( elem.res2_ovu1_overlaps ).xyz().distance( upper_rsd.atom( "OVU1" ).xyz() );
		score_dev += upper_rsd.atom( elem.res1_ovl1_overlaps ).xyz().distance( lower_rsd.atom( "OVL1" ).xyz() );
		if ( upper_rsd.is_carbohydrate() ) {
			// If this is a cutpoint on a carbohydrate chain, we don't care at all about a third virtual atom, because the
			// residues connect at a single bond, not a partial-double as is the case with a peptide bond.
			score_dev *= 1.5;  // Normalize the score with that of non-carbohydrate chain breaks.
		} else {
			score_dev +=
				upper_rsd.atom( elem.res1_ovl2_overlaps ).xyz().distance( lower_rsd.atom( "OVL2" ).xyz() );
		}
		// Scoring
		total_dev += score_dev;
	}

	debug_assert( std::abs( totals[ linear_branch_conn ] ) < 1e-3 );
	totals[ linear_branch_conn ] = total_dev / 3.0;  // average over 3 distances
}


/// called during gradient-based minimization inside dfunc
/**
F1 and F2 are not zeroed -- contributions from this atom are
just summed in
**/
void LinearBranchEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const &, // sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using core::Size;
	using conformation::Residue;

	// This method is called on a per-residue basis, which means that we could
	// be looking at a residue with no cutpoint variants, a lower cutpoint
	// variant, or an upper cutpoint variant. This method must address both
	// possibilities.
	ResidueAtomOverlaps branch_connection;
	find_relevant_connections_onersd( pose, id.rsd(), branch_connection );

	// CASE 1: left-hand side of branch_conn (CUTPOINT_LOWER)
	Size const residue = id.rsd();
	if ( branch_connection.res1 == residue ) {
		Residue const & lower_rsd( pose.residue( branch_connection.res1 ) );
		Residue const & upper_rsd( pose.residue( branch_connection.res2 ) );
		Vector const & xyz_moving( pose.xyz( id ) );

		bool match( true );
		Vector xyz_fixed;
		if ( id.atomno() == lower_rsd.atom_index( branch_connection.res2_ovu1_overlaps ) ) {
			xyz_fixed = upper_rsd.atom( "OVU1" ).xyz();
		} else if ( lower_rsd.atom_name( id.atomno() ) == "OVL1" ) {
			xyz_fixed = upper_rsd.atom( branch_connection.res1_ovl1_overlaps ).xyz();
		} else if ( lower_rsd.atom_name( id.atomno() ) == "OVL2" ) {
			xyz_fixed = upper_rsd.atom( branch_connection.res1_ovl2_overlaps ).xyz();
		} else {
			match = false;
		}

		if ( match ) {
			Vector const f2 ( xyz_moving - xyz_fixed );
			Real const dist ( f2.length() );
			if ( dist >= 1.0e-8 ) { // avoid getting too close to singularity...
				Real const invdist( 1.0 / dist );
				F1 += weights[ linear_branch_conn ] * invdist * cross( xyz_moving, xyz_fixed ) / 3;
				F2 += weights[ linear_branch_conn ] * invdist * ( xyz_moving - xyz_fixed ) / 3;
			}
		}
	}

	// CASE 2: right-hand side of branch_conn (CUTPOINT_UPPER)
	if ( branch_connection.res2 == residue ) {
		Residue const & lower_rsd( pose.residue( branch_connection.res1 ) );
		Residue const & upper_rsd( pose.residue( branch_connection.res2 ) );
		Vector const & xyz_moving( pose.xyz( id ) );

		bool match( true );
		Vector xyz_fixed;
		if ( id.atomno() == upper_rsd.atom_index( branch_connection.res1_ovl1_overlaps ) ) {
			xyz_fixed = lower_rsd.atom( "OVL1" ).xyz();
		} else if ( id.atomno() == upper_rsd.atom_index( branch_connection.res1_ovl2_overlaps ) && ! upper_rsd.is_carbohydrate() ) {
			xyz_fixed = lower_rsd.atom( "OVL2" ).xyz();
		} else if ( upper_rsd.atom_name( id.atomno() ) == "OVU1" ) {
			xyz_fixed = lower_rsd.atom( branch_connection.res2_ovu1_overlaps ).xyz();
		} else {
			match = false;
		}

		if ( match ) {
			Vector const f2 ( xyz_moving - xyz_fixed );
			Real const dist ( f2.length() );
			if ( dist >= 1.0e-8 ) { // avoid getting too close to singularity...
				Real const invdist( 1.0 / dist );
				F1 += weights[ linear_branch_conn ] * invdist * cross( xyz_moving, xyz_fixed ) / 3;
				F2 += weights[ linear_branch_conn ] * invdist * ( xyz_moving - xyz_fixed ) / 3;
			}
		}
	}
}

/// @brief LinearBranch Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
LinearBranchEnergy::indicate_required_context_graphs(utility::vector1<bool>&) const {} /*context_graphs_required*/

core::Size LinearBranchEnergy::version() const {
	return 2;
}

} // namespace methods
} // namespace scoring
} // namespace core
