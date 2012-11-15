// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueCountFilter.cc
/// @brief Filter on the total number of residues in the structure
/// @author Matthew O'Meara (mattjomeara@gmail.com)


//Unit Headers
#include <protocols/simple_filters/ResidueCountFilter.hh>
#include <protocols/simple_filters/ResidueCountFilterCreator.hh>
#include <utility/tag/Tag.hh>
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.hh>
//Project Headers
#include <basic/Tracer.hh>
namespace protocols{
namespace simple_filters {

using namespace core;
using namespace core::scoring;

static basic::Tracer TR( "protocols.simple_filters.ResidueCountFilter" );

protocols::filters::FilterOP
ResidueCountFilterCreator::create_filter() const { return new ResidueCountFilter; }

std::string
ResidueCountFilterCreator::keyname() const { return "ResidueCount"; }

//default ctor
ResidueCountFilter::ResidueCountFilter() :
	protocols::filters::Filter( "ResidueCount" ),
	max_residue_count_(0),
	enable_max_residue_count_(false),
	min_residue_count_(0),
	enable_min_residue_count_(false)
{}

ResidueCountFilter::ResidueCountFilter(
	ResidueCountFilter const & src
) :
	protocols::filters::Filter( "ResidueCount" ),
	max_residue_count_(src.max_residue_count_),
	enable_max_residue_count_(src.enable_max_residue_count_),
	min_residue_count_(src.min_residue_count_),
	enable_min_residue_count_(src.enable_min_residue_count_)
{}

ResidueCountFilter::~ResidueCountFilter() {}

filters::FilterOP
ResidueCountFilter::clone() const {
	return new ResidueCountFilter( *this );
}

filters::FilterOP
ResidueCountFilter::fresh_instance() const {
	return new ResidueCountFilter();
}


void
ResidueCountFilter::parse_my_tag(
	utility::tag::TagPtr const tag,
	moves::DataMap &,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &
) {
	if(tag->hasOption("max_residue_count")){
		enable_max_residue_count();
		max_residue_count(tag->getOption< core::Size >("max_residue_count"));
	}

	if(tag->hasOption("min_residue_count")){
		enable_min_residue_count();
		min_residue_count(tag->getOption< core::Size >("min_residue_count"));
	}
}

bool
ResidueCountFilter::apply(
	core::pose::Pose const & pose
) const {
	if(enable_max_residue_count() && compute(pose) > max_residue_count()){
		return false;
	}

	if(enable_min_residue_count() && compute(pose) < min_residue_count()){
		return false;
	}

	return true;
}

void
ResidueCountFilter::report(
	std::ostream & out,
	core::pose::Pose const & pose
) const {
	out << "Residue Count: " << compute( pose ) <<std::endl;
}

core::Real
ResidueCountFilter::report_sm(
	core::pose::Pose const & pose
) const {
	return compute( pose );
}

core::Real
ResidueCountFilter::compute(
	core::pose::Pose const & pose
) const {
	return pose.total_residue();
}

core::Size
ResidueCountFilter::max_residue_count() const {
	return max_residue_count_;
}

void
ResidueCountFilter::max_residue_count(
	core::Size value
) {
	max_residue_count_ = value;
}

bool
ResidueCountFilter::enable_max_residue_count() const {
	return enable_max_residue_count_;
}

void
ResidueCountFilter::enable_max_residue_count(
	bool value
) {
	enable_max_residue_count_ = value;
}

core::Size
ResidueCountFilter::min_residue_count() const {
	return min_residue_count_;
}

void
ResidueCountFilter::min_residue_count(
	core::Size value
) {
	min_residue_count_ = value;
}

bool
ResidueCountFilter::enable_min_residue_count() const {
	return enable_min_residue_count_;
}

void
ResidueCountFilter::enable_min_residue_count(
	bool value
) {
	enable_min_residue_count_ = value;
}



} // namespace
} // namespace
