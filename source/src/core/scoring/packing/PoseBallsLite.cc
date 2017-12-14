// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packing/PoseBallsLite.cc
/// @brief
/// @author

#include <core/scoring/packing/PoseBallsLite.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>

#include <map>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

namespace core {
namespace scoring {
namespace packing {

/// @details Auto-generated virtual destructor
PoseBallsLite::~PoseBallsLite() = default;

using namespace ObjexxFCL;

static basic::Tracer TR( "core.scoring.packing.PoseBallsLite" );

inline core::Real sqr ( core::Real x ) {
	return x*x;
}

// smoothed neighbor is between 9 and 11 A
inline core::Real sigmoidish_neighbor( core::Real sqdist ) {
	if ( sqdist > 121.0 ) {
		return 0.0;
	} else if ( sqdist < 81.0 ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
		return sqr(1.0 - sqr( (dist - 9.0) / (11.0 - 9.0) ) );
	}
}

/// @details Remove spaces from given string.
inline std::string strip_whitespace( std::string const & name )
{
	std::string trimmed_name( name );
	left_justify( trimmed_name ); trim( trimmed_name ); // simpler way to dothis?
	return trimmed_name;
}


PoseBallsLite::PoseBallsLite(
	core::pose::Pose const & pose,
	core::Size Hmode,
	bool /*ignore_water*/
) {
	using namespace numeric;

	// std::cerr << "PoseBallsLite.cc:67 (" << pose.size() << ")" << std::endl;

	// initialize index and vars
	core::Size index = 0;
	Size num_unrec = 0;
	if ( pose.pdb_info() ) num_unrec = pose.pdb_info()->get_num_unrecognized_atoms();
	balls_      .reserve( pose.size()*5 + num_unrec );
	index_to_id_.reserve( pose.size()*5 + num_unrec );
	core::pose::initialize_atomid_map( id_to_index_, pose );
	id_to_index_.resize( pose.size() + num_unrec );

	// std::cerr << "PoseBallsLite.cc:85 (" << ")" << std::endl;

	// add atoms in pose
	core::Size skippedH = 0;
	for ( core::Size ir = 1; ir <= pose.size(); ++ir ) {
		for ( core::Size ia = 1; ia <= pose.residue(ir).natoms(); ++ia ) {
			// std::cerr << "PoseBallsLite.cc:91 (" << ")" << std::endl;
			if ( pose.residue(ir).is_virtual(ia) ) continue;
			if ( Hmode == 0 ) { // no H's
				if ( pose.residue(ir).atom_type(ia).is_hydrogen() ) {
					skippedH++;
					continue;
				}
			} else if ( Hmode == 1 ) { // polar H's only
				if ( pose.residue(ir).atom_type(ia).is_hydrogen() &&
						!pose.residue(ir).atom_type(ia).is_polar_hydrogen() ) {
					skippedH++;
					continue;
				}
			} else if ( Hmode == 2 ) {
				; // all Hs
			}
			atom_num_.push_back(ia);
			res_num_.push_back(ir);

			// std::cerr << "PoseBallsLite.cc:107 (" << ")" << std::endl;
			core::id::AtomID aid(ia,ir);
			id_to_index_[ aid ] = ++index;
			index_to_id_.push_back( aid );
			balls_.push_back( Ball( pose.xyz(aid), pose.residue(ir).atom_type(ia).lj_radius() ) );
		}
	}
	// std::cerr << "PoseBallsLite.cc:134 (" << ")" << std::endl;

	nballs_ = index;
}


PoseBallsLite::PoseBallsLite(
	core::pose::Pose const & pose,
	core::id::AtomID_Mask const & whichatoms
) {
	using namespace numeric;

	// initialize index and vars
	core::Size index = 0;
	Size num_unrec = 0;
	if ( pose.pdb_info() ) num_unrec = pose.pdb_info()->get_num_unrecognized_atoms();
	balls_      .reserve( pose.size()*5 + num_unrec );
	index_to_id_.reserve( pose.size()*5 + num_unrec );
	core::pose::initialize_atomid_map( id_to_index_, pose );
	id_to_index_.resize( pose.size() + num_unrec );

	// add atoms in pose
	for ( core::Size ir = 1; ir <= pose.size(); ++ir ) {
		for ( core::Size ia = 1; ia <= pose.residue(ir).natoms(); ++ia ) {
			core::id::AtomID aid(ia,ir);
			if ( ! whichatoms[aid] ) continue;
			atom_num_.push_back(ia);
			res_num_.push_back(ir);
			id_to_index_[ aid ] = ++index;
			index_to_id_.push_back( aid );
			balls_.push_back( Ball( pose.xyz(aid), pose.residue(ir).atom_type(ia).lj_radius() ) );
		}
	}

	nballs_ = index;
}


} // namespace packing
} // namespace scoring
} // namespace core
