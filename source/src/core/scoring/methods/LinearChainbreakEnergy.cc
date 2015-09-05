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
/// @author Christopher Miles (cmiles@uw.edu)


#include <core/scoring/methods/chainbreak_util.hh>
#include <core/scoring/methods/LinearChainbreakEnergy.hh>
#include <core/scoring/methods/LinearChainbreakEnergyCreator.hh>

#include <basic/options/keys/OptionKeys.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
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
static thread_local basic::Tracer tr( "core.scoring.LinearChainbreak", basic::t_info );

LinearChainbreakEnergy::LinearChainbreakEnergy()
: parent(methods::EnergyMethodCreatorOP( new LinearChainbreakEnergyCreator )) {
	initialize( std::numeric_limits<core::Size>::max() );
}

LinearChainbreakEnergy::LinearChainbreakEnergy(Size allowable_sequence_sep)
: parent(methods::EnergyMethodCreatorOP( new LinearChainbreakEnergyCreator )) {
	initialize(allowable_sequence_sep);
}

LinearChainbreakEnergy::LinearChainbreakEnergy(const LinearChainbreakEnergy& o)
: parent(o) {
	// assignment to non-cache-related variables
	allowable_sequence_sep_ = o.allowable_sequence_sep_;

	// assign cache-related variables to their default values. this will trigger
	// a cache miss on their first use, but greatly simplifies the design and
	// improves reasoning about code paths in the event of copying/cloning
	shortest_paths_.reset();
	previous_hash_value_ = 0;
}

LinearChainbreakEnergy& LinearChainbreakEnergy::operator=(const LinearChainbreakEnergy& o) {
	if ( this != &o ) {
		// base class assignments
		parent::operator=(o);

		// assignment to non-cache-related variables
		allowable_sequence_sep_ = o.allowable_sequence_sep_;

		// assign cache-related variables to their default values. this will trigger
		// a cache miss on their first use, but greatly simplifies the design and
		// improves reasoning about code paths in the event of copying/cloning
		shortest_paths_.reset();
		previous_hash_value_ = 0;
	}
	return *this;
}

LinearChainbreakEnergy::~LinearChainbreakEnergy() {}

void LinearChainbreakEnergy::initialize(Size allowable_sequence_sep) {
	allowable_sequence_sep_ = allowable_sequence_sep;
	previous_hash_value_ = 0;
	shortest_paths_.reset();
}

core::Real LinearChainbreakEnergy::do_score_dev( core::conformation::Residue const & lower_rsd,
		core::conformation::Residue const & upper_rsd,
		core::Size const nbb ) const
{
	Distance score_dev( 0.0 );

	// virtual N and CA on lower_rsd (OVL1 and OVL2) and a virtual C on upper_rsd "OVU1"
	score_dev += lower_rsd.atom( lower_rsd.mainchain_atoms()[ nbb ] ).xyz().distance( upper_rsd.atom( "OVU1" ).xyz() );
	score_dev += upper_rsd.atom( upper_rsd.mainchain_atoms()[ 1 ] ).xyz().distance( lower_rsd.atom( "OVL1" ).xyz() );
	if ( upper_rsd.is_carbohydrate() ) {
		// If this is a cutpoint on a carbohydrate chain, we don't care at all about a third virtual atom, because the
		// residues connect at a single bond, not a partial-double as is the case with a peptide bond.
		score_dev *= 1.5;  // Normalize the score with that of non-carbohydrate chain breaks.
	} else {
		score_dev +=
				upper_rsd.atom( upper_rsd.mainchain_atoms()[ 2 ] ).xyz().distance( lower_rsd.atom( "OVL2" ).xyz() );
	}

	return score_dev;
}

core::Real LinearChainbreakEnergy::do_score_ovp( core::conformation::Residue const & lower_rsd,
		core::conformation::Residue const & upper_rsd,
		core::Size const nbb,
		core::Size const cutpoint,
		core::pose::Pose const & pose) const
{
	using core::id::AtomID;
	using core::kinematics::Stub;

	debug_assert( ! ( lower_rsd.is_carbohydrate() || upper_rsd.is_carbohydrate() ) );

	// now compute the overlap score: this is done by comparing the stub   (lower side )  C, | N, CA (upper side)
	// the stub for the lower_rsd is already in the atom-tree at atom "OVL2 == CA*"
	Stub lower_stub( pose.conformation().atom_tree().atom(
			AtomID( pose.residue( cutpoint ).atom_index( "OVL2" ), cutpoint ) ).get_stub() );

	// the upper stub... let's just compute it for now...
	// could be gained from AtomID( NamedAtomID( "OVU1", cutpoint+1 ).get_stub()
	// and then the correct reversal...
	Stub upper_stub( upper_rsd.atom( upper_rsd.mainchain_atoms()[ 2 ]).xyz(),  // CA
			upper_rsd.atom( upper_rsd.mainchain_atoms()[ 1 ] ).xyz(),          // N
			upper_rsd.atom( "OVU1" ).xyz() );                                  // virtual C

	//for double-checking... ( debug )
	Stub manual_lower_stub( lower_rsd.atom( "OVL2" ).xyz(),                // virtual CA
			lower_rsd.atom( "OVL1" ).xyz(),                                // virtual N
			lower_rsd.atom( lower_rsd.mainchain_atoms()[ nbb ] ).xyz() );  // C

	if ( distance( lower_stub, manual_lower_stub ) > 0.01 ) {
		tr.Warning << "WARNING: mismatch between manual computed and atom-tree stub: "
				<< lower_stub << " " << manual_lower_stub << std::endl;
	}

	return
		manual_lower_stub.M.col_x().distance(upper_stub.M.col_x()) +
		manual_lower_stub.M.col_y().distance(upper_stub.M.col_y()) +
		manual_lower_stub.M.col_z().distance(upper_stub.M.col_z());
}

/// called at the end of energy evaluation
/// In this case (LinearChainbreakEnergy), all the calculation is done here
void LinearChainbreakEnergy::finalize_total_energy( pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals ) const
{
	using core::Size;
	using conformation::Residue;
	using core::kinematics::FoldTree;
	using core::kinematics::ShortestPathInFoldTree;
	using utility::vector1;

	Real total_dev = 0.0;
	Real total_ovp = 0.0;

	// Cached ShortestPathInFoldTree instance is invalid and must be recomputed
	const FoldTree& tree = pose.fold_tree();
	size_t hash_value = tree.hash_value();
	if ( ! shortest_paths_.get() || previous_hash_value_ != hash_value ) {
		shortest_paths_.reset( new ShortestPathInFoldTree( tree ) );
		previous_hash_value_ = hash_value;
	}

	// Identify all cutpoint variants defined by the caller
	vector1< int > cutpoints;
	find_cutpoint_variants( pose, tree, &cutpoints );

	for ( Size i = 1; i <= cutpoints.size(); ++i ) {
		int const cutpoint = cutpoints[ i ];
		Residue const & lower_rsd = pose.residue( cutpoint );
		Residue const & upper_rsd = pose.residue( cutpoint + 1 );
		Size const nbb = lower_rsd.mainchain_atoms().size();

		if ( ! lower_rsd.has_variant_type( core::chemical::CUTPOINT_LOWER ) ||
				! upper_rsd.has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
			continue;
		}

		// Determine whether the separation between <lowed_rsd> and <upper_rsd>,
		// as computed by ShortestPathInFoldTree, exceeds the current allowable
		// sequence separation
		Size const separation = shortest_paths_->dist( cutpoint, cutpoint + 1 );
		if ( separation > allowable_sequence_sep_ ) {
			tr.Trace << "Chainbreak skipped-- "
				<< separation << " > " << allowable_sequence_sep_
				<< std::endl;
			continue;
		}

		// Scoring
		total_dev += do_score_dev( lower_rsd, upper_rsd, nbb );
		if ( ! lower_rsd.is_carbohydrate() ) {
			// An overlap_chainbreak score cannot currently be calculated for sugars.
			total_ovp += do_score_ovp( lower_rsd, upper_rsd, nbb, cutpoint, pose );
		}
	}

	debug_assert( std::abs( totals[ linear_chainbreak ] ) < 1e-3 );
	totals[ linear_chainbreak ] = total_dev / 3.0;  // average over 3 distances
	totals[ overlap_chainbreak ] = total_ovp;
}


/// called during gradient-based minimization inside dfunc
/**
F1 and F2 are not zeroed -- contributions from this atom are
just summed in
**/
void LinearChainbreakEnergy::eval_atom_derivative(
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

	// CASE 1: left-hand side of chainbreak (CUTPOINT_LOWER)
	Size const residue = id.rsd();
	if ( is_lower_cutpoint( residue, pose ) ) {
		Residue const & lower_rsd( pose.residue( residue ) );
		Residue const & upper_rsd( pose.residue( residue + 1 ) );
		Vector const & xyz_moving( pose.xyz( id ) );

		bool match( true );
		Vector xyz_fixed;
		Size const nbb( lower_rsd.mainchain_atoms().size() );
		if ( id.atomno() == lower_rsd.mainchain_atoms()[ nbb ] ) {
			xyz_fixed = upper_rsd.atom( "OVU1" ).xyz();
		} else if ( lower_rsd.atom_name( id.atomno() ) == "OVL1" ) {
			xyz_fixed = upper_rsd.atom( upper_rsd.mainchain_atoms()[ 1 ] ).xyz();
		} else if ( lower_rsd.atom_name( id.atomno() ) == "OVL2" ) {
			xyz_fixed = upper_rsd.atom( upper_rsd.mainchain_atoms()[ 2 ] ).xyz();
		} else {
			match = false;
		}

		if ( match ) {
			Vector const f2 ( xyz_moving - xyz_fixed );
			Real const dist ( f2.length() );
			if ( dist >= 0.01 ) { // avoid getting too close to singularity...
				Real const invdist( 1.0 / dist );
				F1 += weights[ linear_chainbreak ] * invdist * cross( xyz_moving, xyz_fixed ) / 3;
				F2 += weights[ linear_chainbreak ] * invdist * ( xyz_moving - xyz_fixed ) / 3;
			}
		}
	}

	// CASE 2: right-hand side of chainbreak (CUTPOINT_UPPER)
	if ( is_upper_cutpoint( residue,pose ) ) {
		Residue const & lower_rsd( pose.residue( residue - 1 ) );
		Residue const & upper_rsd( pose.residue( residue ) );
		Vector const & xyz_moving( pose.xyz( id ) );

		bool match( true );
		Vector xyz_fixed;
		Size const nbb( lower_rsd.mainchain_atoms().size() );
		if ( id.atomno() == upper_rsd.mainchain_atoms()[ 1 ] ) {
			xyz_fixed = lower_rsd.atom( "OVL1" ).xyz();
		} else if ( id.atomno() == upper_rsd.mainchain_atoms()[ 2 ] && ! upper_rsd.is_carbohydrate() ) {
			xyz_fixed = lower_rsd.atom( "OVL2" ).xyz();
		} else if ( upper_rsd.atom_name( id.atomno() ) == "OVU1" ) {
			xyz_fixed = lower_rsd.atom( lower_rsd.mainchain_atoms()[ nbb ] ).xyz();
		} else {
			match = false;
		}

		if ( match ) {
			Vector const f2 ( xyz_moving - xyz_fixed );
			Real const dist ( f2.length() );
			if ( dist >= 0.01 ) { // avoid getting too close to singularity...
				Real const invdist( 1.0 / dist );
				F1 += weights[ linear_chainbreak ] * invdist * cross( xyz_moving, xyz_fixed ) / 3;
				F2 += weights[ linear_chainbreak ] * invdist * ( xyz_moving - xyz_fixed ) / 3;
			}
		}
	}
}

/// @brief LinearChainbreak Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
LinearChainbreakEnergy::indicate_required_context_graphs(utility::vector1<bool>&) const {} /*context_graphs_required*/

core::Size LinearChainbreakEnergy::version() const {
	return 2;
}

} // namespace methods
} // namespace scoring
} // namespace core
