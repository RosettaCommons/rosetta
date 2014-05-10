// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueDistanceFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#include <protocols/simple_filters/ResidueDistanceFilter.hh>
#include <protocols/simple_filters/ResidueDistanceFilterCreator.hh>

#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <basic/Tracer.hh>


namespace protocols {
namespace simple_filters {

static basic::Tracer residue_distance_filter_tracer( "protocols.simple_filters.ResidueDistanceFilter" );

protocols::filters::FilterOP
ResidueDistanceFilterCreator::create_filter() const { return new ResidueDistanceFilter; }

std::string
ResidueDistanceFilterCreator::keyname() const { return "ResidueDistance"; }

ResidueDistanceFilter::~ResidueDistanceFilter(){}

void
ResidueDistanceFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
	res1_ = core::pose::get_resnum( tag, pose, "res1_" );
	res2_ = core::pose::get_resnum( tag, pose, "res2_" );
	distance_threshold_ = tag->getOption<core::Real>( "distance", 8.0 );

	residue_distance_filter_tracer<<"ResidueDistanceFilter with distance threshold of "<<distance_threshold_<<" between residues "<<res1_<<" and "<<res2_<<std::endl;
}

bool
ResidueDistanceFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const distance( compute( pose ) );

	residue_distance_filter_tracer<<"Distance between residues "<<pose.residue( res1_ ).name3()<<res1_<<" and "<<pose.residue( res2_ ).name3()<<res2_<<" is "<<distance<<std::endl;
	return( distance<=distance_threshold_ );
}

void
ResidueDistanceFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const distance( compute( pose ) );

	out<<"Distance between residues "<<pose.residue( res1_ ).name3()<<res1_<<" and "<<pose.residue( res2_ ).name3()<<res2_<<" is "<<distance<<'\n';
}

core::Real
ResidueDistanceFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const distance( compute( pose ) );

	return( distance );
}
core::Real
ResidueDistanceFilter::compute( core::pose::Pose const & pose ) const {
	core::conformation::Residue const res_res1( pose.conformation().residue( res1_ ) );
	core::conformation::Residue const res_res2( pose.conformation().residue( res2_ ) );
	core::Real const distance( res_res1.xyz( res_res1.nbr_atom() ).distance( res_res2.xyz( res_res2.nbr_atom() ) ) );
	return( distance );
}

}
}
