// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/atomic_depth/util.hh
/// @brief  Util functions for atomic_depth calculations.
/// @author Brian Coventry


#include <core/scoring/atomic_depth/util.hh>

#include <core/pose/util.tmpl.hh>


namespace core {
namespace scoring {
namespace atomic_depth {


/// @brief Calculate depth of all atoms from edge of Sasa surface.
core::id::AtomID_Map< core::Real >
atomic_depth(
	pose::Pose const & pose,
	Real probe_radius, /*= 1.4*/
	bool poly_leu_depth /*= false*/
) {
	AtomicDepthOP depth = nullptr;
	return atomic_depth( pose, depth, probe_radius, poly_leu_depth );
}

/// @brief Calculate depth of all atoms from edge of Sasa surface.
/// @detail If atomic_depth is a nullptr, will create one and return it.
/// @detail  Note that probe_radius is ignored if atomic_depth is not null
core::id::AtomID_Map< core::Real >
atomic_depth(
	pose::Pose const & pose,
	AtomicDepthOP & depth,
	Real probe_radius, /*= 1.4*/
	bool poly_leu_depth /*= false*/
) {
	core::id::AtomID_Map< bool > depth_atoms;
	pose::initialize_atomid_map( depth_atoms, pose.conformation(), true );

	return atomic_depth( pose, depth_atoms, depth, probe_radius, poly_leu_depth );
}

/// @brief Calculate depth of certain atoms from edge of Sasa surface.
core::id::AtomID_Map< core::Real >
atomic_depth(
	pose::Pose const & pose,
	core::id::AtomID_Map< bool > depth_atoms,
	Real probe_radius, /*= 1.4*/
	bool poly_leu_depth /*= false*/
) {
	AtomicDepthOP depth = nullptr;
	return atomic_depth( pose, depth_atoms, depth, probe_radius, poly_leu_depth );
}

/// @brief Calculate depth of certain atoms from edge of Sasa surface.
/// @detail If atomic_depth is a nullptr, will create one and return it.
/// @detail  Note that probe_radius is ignored if atomic_depth is not null
core::id::AtomID_Map< core::Real >
atomic_depth(
	pose::Pose const & pose,
	core::id::AtomID_Map< bool > depth_atoms,
	AtomicDepthOP & depth,
	Real probe_radius, /*= 1.4*/
	bool poly_leu_depth /*= false*/
) {
	core::id::AtomID_Map< core::Real > final_depths;
	pose::initialize_atomid_map( final_depths, pose.conformation(), -1.0 );

	if ( pose.size() == 0 ) return final_depths;

	if ( ! depth ) {
		depth = AtomicDepthOP( new AtomicDepth( pose, probe_radius, poly_leu_depth ) );
	}
	utility::vector1< conformation::Atom > atoms;

	// This for loop must remain identical to the other one
	for ( core::Size seqpos = 1; seqpos <= final_depths.size(); seqpos++ ) {
		for ( core::Size atno = 1; atno < final_depths.n_atom( seqpos ); atno++ ) {
			if ( depth_atoms( seqpos, atno ) ) {
				atoms.push_back( pose.residue(seqpos).atom(atno) );
			}
		}
	}

	chemical::AtomTypeSet const & type_set = pose.residue(1).type().atom_type_set();
	utility::vector1< Real > depths = depth->calcdepth( atoms, type_set );

	// This for loop must remain identical to the other one
	Size current = 1;
	for ( core::Size seqpos = 1; seqpos <= final_depths.size(); seqpos++ ) {
		for ( core::Size atno = 1; atno < final_depths.n_atom( seqpos ); atno++ ) {
			if ( depth_atoms( seqpos, atno ) ) {
				final_depths( seqpos, atno ) = depths.at( current++ );
			}
		}
	}
	runtime_assert( current-1 == depths.size() );

	return final_depths;
}

/// @brief Find all atoms deeper than threshold from the edge of the sasa surface.
core::id::AtomID_Map< bool >
atoms_deeper_than(
	pose::Pose const & pose,
	Real threshold,
	bool invert /*= false*/,
	Real probe_radius /*= 1.4*/,
	bool poly_leu_depth /*= false*/
) {
	core::id::AtomID_Map< core::Real > depths = atomic_depth( pose, probe_radius, poly_leu_depth );

	core::id::AtomID_Map< bool > is_deep;
	pose::initialize_atomid_map( is_deep, pose.conformation(), false );

	for ( core::Size seqpos = 1; seqpos <= depths.size(); seqpos++ ) {
		for ( core::Size atno = 1; atno < depths.n_atom( seqpos ); atno++ ) {
			bool this_is_deep = depths( seqpos, atno ) >= threshold;

			is_deep( seqpos, atno ) = this_is_deep ^ invert;
		}
	}

	return is_deep;
}




} // namespace atomic_depth
} // namespace scoring
} // namespace core
