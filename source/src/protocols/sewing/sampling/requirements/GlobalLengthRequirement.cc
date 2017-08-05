// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file GlobalLengthRequirement.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/sewing/sampling/requirements/GlobalLengthRequirement.hh>
#include <protocols/sewing/sampling/requirements/GlobalLengthRequirementCreator.hh>

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

static THREAD_LOCAL basic::Tracer TR("protocols.sewing.sampling.requirements.GlobalLengthRequirement");

//////Creator methods/////
GlobalRequirementOP
GlobalLengthRequirementCreator::create_requirement() const{
	return GlobalRequirementOP( new GlobalLengthRequirement() );
}

std::string
GlobalLengthRequirementCreator::type_name() const {
	return "GlobalLengthRequirement";
}
//////End Creator methods/////

GlobalLengthRequirement::GlobalLengthRequirement():
	min_length_(0),
	max_length_(-1)
{}

GlobalLengthRequirement::GlobalLengthRequirement(
	std::set<std::string> valid_dssp_codes,
	core::Size min_length,
	core::Size max_length
):
	dssp_codes_(valid_dssp_codes),
	min_length_(min_length),
	max_length_(max_length)
{}

bool
GlobalLengthRequirement::satisfies(
	AssemblyCOP assembly
) const {
	return !violates(assembly);
}

bool
GlobalLengthRequirement::violates(
	AssemblyCOP assembly
) const {
	utility::vector1<SewSegment> const & segments = assembly->segments();
	for ( core::Size i = 1; i <= segments.size(); ++i ) {
		if ( dssp_codes_.size() == 0 || dssp_codes_.find(std::string(1, segments[i].dssp_)) != dssp_codes_.end() ) {
			if ( segments[i].residues_.size() < min_length_ || segments[i].residues_.size() > max_length_ ) {
				if ( TR.Debug.visible() ) {
					TR.Debug << "Segment " << i << ", with dssp " << segments[i].dssp_ << " and "
						<< segments[i].residues_.size() << " residues fails GlobalLengthRequirement" << std::endl;
				}
				return true;
			}
		}
	}
	return false;
}

void
GlobalLengthRequirement::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){
	if ( tag->hasOption("dssp") ) {
		dssp_codes_ = utility::split_to_set(tag->getOption<std::string>("dssp"));
	}

	if ( tag->hasOption("min_length") ) {
		min_length_ = tag->getOption<core::Size>("min_length");
	}

	if ( tag->hasOption("max_length") ) {
		max_length_ = tag->getOption<core::Size>("max_length");
	}
}


void
GlobalLengthRequirement::show(
	std::ostream & out
) const {
	out << "/////// GlobalLengthRequirement - ";
	if ( dssp_codes_.size() > 0 ) {
		out << "Segments with DSSP codes ";
		std::set<std::string>::const_iterator it = dssp_codes_.begin();
		std::set<std::string>::const_iterator it_end = dssp_codes_.end();
		for ( ; it != it_end; ++it ) {
			out << *it << ",";
		}
	} else {
		out << "All segments";
	}
	out << " must have between " << min_length_ << " and " << max_length_ << " residues." << std::endl;
}


} //requirements namespace
} //sampling namespace
} //sewing namespace
} //protocols namespace
