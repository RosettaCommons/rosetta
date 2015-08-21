// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/OutputWriter.cc
/// @brief  Implementation for the class to write output matches that have (presumably) passed the output filters.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), Florian Richter (floric@u.washington.edu) porting to mini

// Unit headers
#include <protocols/match/output/OutputWriter.hh>

#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
#include <protocols/match/MatcherTask.hh>
#include <protocols/match/Hit.hh>

//utility headers
#include <utility/string_util.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {

OutputWriter::OutputWriter()
: cst_io_(/* NULL */)
{}

OutputWriter::~OutputWriter() {}

protocols::toolbox::match_enzdes_util::EnzConstraintIOCOP
OutputWriter::cst_io() const
{
	return cst_io_;
}


void
OutputWriter::initialize_from_matcher_task(
	MatcherTaskCOP mtask
){
	cst_io_ = mtask->enz_input_data();
}


void
OutputWriter::determine_redundant_upstream_matchres(
	match_dspos1 const & m,
	std::map< core::Size, core::Size > & redundant_upstream_res
) const
{
	using namespace core;
	using namespace toolbox::match_enzdes_util;

	redundant_upstream_res.clear();
	std::map< Size, Size > seen_scafpos_to_geom_cst;
	std::map< Size, Size > bb_scafpos_to_geom_cst;

	for ( core::Size i = 1; i <= m.upstream_hits.size(); ++i ) {
		core::Size this_scafpos = m.upstream_hits[i].scaffold_build_id();
		MatchConstraintFileInfoCOP cur_mcfi =
			cst_io_->mcfi_list( i )->mcfi( m.upstream_hits[i].external_geom_id() );
		bool bb_interaction( cur_mcfi->is_backbone( cur_mcfi->upstream_res() ) );

		std::map< Size, Size >::iterator seen_it = seen_scafpos_to_geom_cst.find( this_scafpos );
		std::map< Size, Size >::iterator bb_seen_it = bb_scafpos_to_geom_cst.find( this_scafpos );

		//first case: this residue position has already been seen,
		//so either this or the previously seen constraint are redundant
		if ( seen_it != seen_scafpos_to_geom_cst.end() ) {
			if ( bb_interaction ) {
				redundant_upstream_res[ i ] = seen_it->second;
				continue;
			} else {
				if ( bb_seen_it != bb_scafpos_to_geom_cst.end() ) {
					redundant_upstream_res[ bb_seen_it->second ] = i;
				} else {
					utility_exit_with_message("Error in matcher: two sidechains (for geomcsts "+utility::to_string( i )+" and " +utility::to_string( seen_it->second )+") placed at the same scaffold location ("+utility::to_string( this_scafpos)+") considered a match. Hint: if one of these geomcsts is supposed to be a backbone only interaction, this needs to be specified in the cstfile through the \"TEMPLATE::   ATOM_MAP: <residue num> is_backbone \" tag.");
				}
			}
		} else {
			//otherwise this position has not been encountered before
			seen_scafpos_to_geom_cst[ this_scafpos] = i;
			if ( bb_interaction ) {
				bb_scafpos_to_geom_cst[ this_scafpos ] = i;
			}
		}
	} //loop over all matchres
}

}
}
}
