// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/rna/RNA_Motif.cc
/// @brief
/// @detailed Most of the good stuff is in RNA_Motif.hh [for speed]. Some util functions used by rna_motif app in here.
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/rna/RNA_Motif.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <utility/string_util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.scoring.rna.RNA_Motif" );

using namespace core::pose;

namespace core {
namespace scoring {
namespace rna {

/////////////////////////////////////////////////////////////////
/// @details used in RiboDraw (which takes output from rna_motif app)
void
output_rna_motifs(
	std::ostream & out,
	core::pose::Pose const & pose,
	RNA_Motifs const & rna_motifs,
	bool const output_WC_stacked_pair/* = false */ )
{
	using utility::vector1;
	for ( auto const & motif : rna_motifs ) {

		if ( !output_WC_stacked_pair ) {
			if ( rna_motif_bonus.count( motif.type() ) == 0 ) continue;
			if ( rna_motif_bonus.find( motif.type() )->second >= 0.0  ) continue;
		}

		vector1< int > res;
		vector1< char > chain;
		vector1< std::string > segid;
		for ( auto const & m : motif ) {
			res.push_back( pose.pdb_info()->number( m ) );
			chain.push_back( pose.pdb_info()->chain( m ) );
			segid.push_back( pose.pdb_info()->segmentID( m ) );
		}
		out << to_string( motif.type() ) << " " << make_tag_with_dashes( res, chain, segid ) << std::endl;
	}
}

// @brief super-simple helper function for PyMOL commands.
// @details for speed, could also define PyMOL boolean property and then color things at end based on property.
void
output_motifs_to_pymol( std::ostream & out,
	pose::Pose const & pose,
	RNA_Motifs const & rna_motifs ) {

	std::string tag( tag_from_pose( pose ) );
	tag = utility::replace_in( tag, ".pdb", "" );
	for ( auto const & motif : rna_motifs ) {
		for ( auto const & res : motif ) {
			out << "color " << motif_color.find( motif.type() )->second << ", " << tag << " and chain " << pose.pdb_info()->chain( res ) << " and resi " << pose.pdb_info()->number( res ) << std::endl;

		}
	}
}

} //rna
} //scoring
} //core
