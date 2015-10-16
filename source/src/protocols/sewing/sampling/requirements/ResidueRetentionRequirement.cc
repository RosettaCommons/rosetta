// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ResidueRetentionRequirement.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/sewing/sampling/requirements/ResidueRetentionRequirement.hh>
#include <protocols/sewing/sampling/requirements/ResidueRetentionRequirementCreator.hh>

//Package headers
#include <protocols/sewing/conformation/Assembly.hh>

//Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace sewing  {
namespace sampling {
namespace requirements {

static basic::Tracer TR("protocols.sewing.sampling.requirements.ResidueRetentionRequirement");

//////Creator methods/////
GlobalRequirementOP
ResidueRetentionRequirementCreator::create_requirement() const{
	return GlobalRequirementOP( new ResidueRetentionRequirement() );
}

std::string
ResidueRetentionRequirementCreator::type_name() const {
	return "ResidueRetentionRequirement";
}
//////End Creator methods/////

ResidueRetentionRequirement::ResidueRetentionRequirement():
	model_id_(0)
{}

ResidueRetentionRequirement::ResidueRetentionRequirement(
	int model_id
):
	model_id_(model_id)
{}

void
ResidueRetentionRequirement::model_id(
	int model_id
){
	model_id_ = model_id;
}

void
ResidueRetentionRequirement::required_resnums(
	std::set<core::Size> required_resnums
){
	required_resnums_ = required_resnums;
}

void
ResidueRetentionRequirement::add_resnum(
	core::Size resnum
){
	required_resnums_.insert(resnum);
}


bool
ResidueRetentionRequirement::satisfies(
	AssemblyCOP assembly
) const {
	return !violates(assembly);
}

bool
ResidueRetentionRequirement::violates(
	AssemblyCOP assembly
) const {
	std::set<core::Size>::const_iterator it = required_resnums_.begin();
	std::set<core::Size>::const_iterator it_end = required_resnums_.end();
	for ( ; it != it_end; ++it ) {

		bool found = false;
		utility::vector1<SewSegment> const & segments = assembly->segments();
		for ( core::Size i = 1; i <= segments.size(); ++i ) {
			if ( segments[i].model_id_ == model_id_ ) {
				for ( core::Size j = 1; j <= segments[i].residues_.size(); ++j ) {
					if ( segments[i].residues_[j].resnum_ == *it ) {
						found = true;
						break;
					}
				}
			}
		}

		if ( !found ) {
			if ( TR.Debug.visible() ) {
				TR.Debug << "Failed to find resnum " << *it << " for model " << model_id_ << std::endl;
			}
			return true;
		}
	}
	return false;
}

void
ResidueRetentionRequirement::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){

	if ( !tag->hasOption("model_id") ) {
		utility_exit_with_message("ResidueRetenetionRequirement tag requires a 'model_id' option");
	}
	model_id_ = tag->getOption<int>("model_id");

	if ( !tag->hasOption("required_resnums") ) {
		utility_exit_with_message("ResidueRetenetionRequirement tag requires a whitespace delimited 'resnums' option");
	}
	utility::vector1<std::string> required_resnums_strings = utility::split(tag->getOption<std::string>("resquired_resnums"));
	for ( core::Size i=1; i<=required_resnums_strings.size(); ++i ) {
		std::stringstream ss ( std::stringstream::in | std::stringstream::out );
		core::Size resnum;
		ss << required_resnums_strings[i];
		ss >> resnum;
		required_resnums_.insert(resnum);
	}
}


void
ResidueRetentionRequirement::show(
	std::ostream & out
) const {
	out << "/////// ResidueRetentionRequirement - Assembly must have a segment with ID " << model_id_
		<< " and residues ";
	std::set<core::Size>::const_iterator it = required_resnums_.begin();
	std::set<core::Size>::const_iterator it_end = required_resnums_.end();
	for ( ; it != it_end; ++it ) {
		out << " " << *it;
	}
	out << std::endl;
}


} //requirements namespace
} //sampling namespace
} //sewing namespace
} //protocols namespace
