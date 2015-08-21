// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/AveragePathLengthFilter.cc
/// @brief Filter on the average covalent path length between residues, including disulfides
/// @author Gabriel Rocklin (grocklin@gmail.com)

//Unit Headers
#include <protocols/simple_filters/AveragePathLengthFilter.hh>
#include <protocols/simple_filters/AveragePathLengthFilterCreator.hh>
#include <utility/tag/Tag.hh>
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.hh>
//Project Headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>
//Utility Headers
#include <utility/string_util.hh>
#include <protocols/toolbox/NetworkAlgorithms.hh>

namespace protocols {
namespace simple_filters {

using namespace core;
using namespace core::scoring;

static thread_local basic::Tracer TR( "protocols.simple_filters.AveragePathLengthFilter" );

protocols::filters::FilterOP
AveragePathLengthFilterCreator::create_filter() const { return protocols::filters::FilterOP( new AveragePathLengthFilter ); }

std::string
AveragePathLengthFilterCreator::keyname() const { return "AveragePathLength"; }


//default ctor
AveragePathLengthFilter::AveragePathLengthFilter() :
	protocols::filters::Filter( "AveragePathLength" ),
	path_tightness_(1.0),
	max_path_length_(10000.0)
{}

AveragePathLengthFilter::AveragePathLengthFilter(core::Real path_tightness, core::Real max_path_length) :
	protocols::filters::Filter( "AveragePathLength" ),
	path_tightness_(path_tightness),
	max_path_length_(max_path_length)
{}

AveragePathLengthFilter::AveragePathLengthFilter(
	AveragePathLengthFilter const & src
) :
	protocols::filters::Filter( "AveragePathLength" ),
	path_tightness_(src.path_tightness_),
	max_path_length_(src.max_path_length_)

{}

AveragePathLengthFilter::~AveragePathLengthFilter() {}

filters::FilterOP
AveragePathLengthFilter::clone() const {
	return filters::FilterOP( new AveragePathLengthFilter( *this ) );
}

filters::FilterOP
AveragePathLengthFilter::fresh_instance() const {
	return filters::FilterOP( new AveragePathLengthFilter() );
}

core::Real
AveragePathLengthFilter::max_path_length() const {
	return max_path_length_;
}

void
AveragePathLengthFilter::max_path_length(
	core::Real value
) {
	max_path_length_ = value;
}

core::Real
AveragePathLengthFilter::path_tightness() const {
	return path_tightness_;
}

void
AveragePathLengthFilter::path_tightness(
	core::Real value
) {
	path_tightness_ = value;
}


void
AveragePathLengthFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &
) {

	if ( tag->hasOption("path_tightness") ) {
		path_tightness(tag->getOption< core::Real >("path_tightness",1));
		TR << "setting path_tightness" << std::endl;
	}


	if ( tag->hasOption("max_path_length") ) {
		max_path_length(tag->getOption< core::Real >("max_path_length",99999));
		TR << "setting max_path_length" << std::endl;
	}

}


bool
AveragePathLengthFilter::apply(
	core::pose::Pose const & pose
) const {

	//calculate the path length
	core::Real path_length = compute(pose);

	//must be below the maximum path length
	if ( path_length > max_path_length_ ) {
		TR << "Failed average path length filter (current: " << path_length << " | user inputted maximum " << max_path_length_ << ")" << std::endl;
		return false;
	}

	//must also be below the path tightness.
	//calculate number of disulfides
	int n_disulfides = 0;
	for ( core::Size i=1; i != pose.total_residue(); ++i ) {
		for ( core::Size j=i + 2; j < pose.total_residue() + 1; ++j ) {
			if ( pose.residue(i).is_bonded(pose.residue(j)) ) {
				n_disulfides++;
			}
		}
	}

	core::Real expected_path_length = ((0.1429 * pose.total_residue()) + 0.8635);
	if ( path_length > (expected_path_length + path_tightness_) ) {
		TR << "Failed average path length filter (current: " << path_length << " | expected base ";
		TR << expected_path_length << " + tightness " << path_tightness_ << ")" << std::endl;
		return false;
	}

	TR << "Average path length filter success (current: " << path_length << " | expected base ";
	TR << expected_path_length << " + tightness " << path_tightness_ << ")" << std::endl;

	return true;
}

void
AveragePathLengthFilter::report(
	std::ostream & out,
	core::pose::Pose const & pose
) const {
	out << "Average path length: " << compute( pose ) <<std::endl;
}

core::Real
AveragePathLengthFilter::report_sm(
	core::pose::Pose const & pose
) const {
	return compute( pose );
}

core::Real
AveragePathLengthFilter::compute(
	core::pose::Pose const & pose
) const {

	//generate nodes
	protocols::toolbox::CovalentResidueNetwork network;
	network.create_from_pose(pose);

	return network.average_shortest_path_length();
}


} // namespace
} // namespace
