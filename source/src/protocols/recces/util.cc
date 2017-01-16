// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/recces/util.cc
/// @brief
/// @details
///
/// @author Andy Watkins


#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/io/ozstream.hh>

#include <basic/Tracer.hh>

//// C++ headers
#include <string>
#include <map>
#include <cmath>
#include <utility/vector1.hh>

static basic::Tracer TR( "protocols.recces.util" );

using namespace core;

namespace protocols {
namespace recces {

// score types to be recorded
// This is so bad, let's fix it -- rhiju, 2016
utility::vector1<core::scoring::ScoreType> const & get_scoretypes() {
	using namespace core::scoring;
	static utility::vector1< ScoreType > scoretypes;
	if ( !scoretypes.empty() ) return scoretypes;

	// List of score types to be cached
	scoretypes.push_back( fa_atr );
	scoretypes.push_back( fa_rep );
	scoretypes.push_back( fa_intra_rep );
	scoretypes.push_back( fa_stack );
	scoretypes.push_back( rna_torsion );
	scoretypes.push_back( hbond_sc );
	scoretypes.push_back( lk_nonpolar );
	scoretypes.push_back( geom_sol_fast );
	scoretypes.push_back( stack_elec );
	scoretypes.push_back( fa_elec_rna_phos_phos );
	return scoretypes;
}

////////////////////////////////////////////////////////////////////////////////
core::Size data_dim( utility::vector1< core::scoring::ScoreType > const & score_types ) {
	return ( score_types.size() + 2 );
}

////////////////////////////////////////////////////////////////////////////////
core::Size data_dim() {
	using namespace core::scoring;
	utility::vector1<ScoreType> const & score_types( get_scoretypes() );
	return data_dim( score_types );
}


//////////////////////////////////////////////////////////////////////////////
void update_scores(
	utility::vector1<float> & scores,
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scorefxn,
	utility::vector1< core::scoring::ScoreType > const & score_types
) {
	using namespace core::scoring;
	scores.clear();
	scores.push_back( ( *scorefxn )( pose ) );
	for ( core::Size i = 1; i<= score_types.size(); ++i ) {
		scores.push_back( scorefxn->score_by_scoretype( pose, score_types[i], false /*weighted*/ ) );
	}
}

//////////////////////////////////////////////////////////////////////////////
void update_scores(
	utility::vector1<float> & scores,
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scorefxn )
{
	// duh, should be able to read scoretypes from scorefxn -- this
	// is here to match legacy code and will be replaced.
	utility::vector1< core::scoring::ScoreType > const & scoretypes( get_scoretypes() );
	update_scores( scores, pose, scorefxn, scoretypes );
}

//////////////////////////////////////////////////////////////////////////////
void fill_data(
	utility::vector1<float> & data,
	core::Size const count,
	utility::vector1<float> const & scores
) {
	data.push_back( count );
	data.append( scores );
}

//////////////////////////////////////////////////////////////////
// @brief used to compute moments of inertia, phase space volume
void
print_base_centroid_atoms_for_rb_entropy( core::conformation::Residue const & rsd, std::string filename_xyz  )
{
	utility::io::ozstream out_xyz;
	out_xyz.open( filename_xyz );
	for ( Size i = rsd.first_sidechain_atom() + 1; i <= rsd.nheavyatoms(); ++i ) { //rsd.first_sidechain_atom()+1 to not include the O2prime oxygen.
		//  TR << rsd.atom_name( i ) << std::endl;
		if ( rsd.is_virtual( i ) ) continue;
		out_xyz << rsd.xyz(i).x() << " " << rsd.xyz(i).y() << " " << rsd.xyz(i).z() << std::endl;
	}
	out_xyz.close();
	TR << "Outputted xyz values into " << filename_xyz << std::endl;
}


//////////////////////////////////////////////////////////////////////////////
utility::vector1<core::Real>
get_torsions(
	utility::vector1<core::id::TorsionID> const & torsion_ids,
	core::pose::Pose const & pose
) {
	utility::vector1<core::Real> curr_torsions;
	for ( auto const & torsion_id: torsion_ids )  curr_torsions.push_back( pose.torsion( torsion_id )  );
	return curr_torsions;
}

} //recces
} //protocols
