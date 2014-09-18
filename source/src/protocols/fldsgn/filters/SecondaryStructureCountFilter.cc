// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/SecondaryStructureCountFilter.cc
/// @brief filter structures by number of secondary structures
/// @detailed
/// @author Lei Shi ( shilei@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/SecondaryStructureCountFilter.hh>
#include <protocols/fldsgn/filters/SecondaryStructureCountFilterCreator.hh>

// Project Headers
#include <protocols/jd2/parser/BluePrint.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <protocols/fldsgn/topology/HSSTriplet.hh> // REQUIRED FOR WINDOWS

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <core/scoring/dssp/Dssp.hh>


//// C++ headers
static thread_local basic::Tracer tr( "protocols.fldsgn.filters.SecondaryStructureCountFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
SecondaryStructureCountFilter::SecondaryStructureCountFilter( ): Filter( "SecondaryStructureCount" )
{}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool SecondaryStructureCountFilter::apply( core::pose::Pose const & pose ) const
{
	compute( pose );

	bool helix_filter=false;
	bool sheet_filter=false;
	bool loop_filter=false;
	bool helix_sheet_filter=false;
	
	if ( filter_helix_ ) {
			if  ( num_helix_pose_ >= num_helix_ )
					helix_filter=true;			
	} else {
			helix_filter=true;
	}

  if ( filter_sheet_ ) {
      if  ( num_sheet_pose_ >= num_sheet_ )
          sheet_filter=true;
  } else {
      sheet_filter=true;
  }

  if ( filter_helix_sheet_ ) {
      if  ( num_helix_pose_ + num_sheet_pose_ >= num_helix_sheet_ )
          helix_sheet_filter=true;
  } else {
      helix_sheet_filter=true;
  }

  if ( filter_loop_ ) {
      if  ( num_loop_pose_ >= num_loop_ )
          loop_filter=true;
  } else {
      loop_filter=true;
  }

	if ( helix_filter && sheet_filter && helix_sheet_filter && loop_filter )
		  return true;
	else
			return false;
} // apply_filter

/// @brief parse xml
void
SecondaryStructureCountFilter::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	using protocols::jd2::parser::BluePrint;
	num_helix_ = tag->getOption<core::Size>( "num_helix", 0 );
	num_sheet_ = tag->getOption<core::Size>( "num_sheet", 0 );
	num_loop_ = tag->getOption<core::Size>( "num_loop", 0 );
	min_helix_length_ = tag->getOption<core::Size>( "min_helix_length", 4 );
	min_sheet_length_ = tag->getOption<core::Size>( "min_sheet_length", 3 );
	min_loop_length_ = tag->getOption<core::Size>( "min_loop_length", 0 );
	max_helix_length_ = tag->getOption<core::Size>( "max_helix_length", 9999 );
	max_sheet_length_ = tag->getOption<core::Size>( "max_sheet_length", 9999 );
	max_loop_length_ = tag->getOption<core::Size>( "max_loop_length", 9999 );
	num_helix_sheet_ = tag->getOption<core::Size>( "num_helix_sheet", 0);
	filter_helix_ = tag->getOption<bool>( "filter_helix", 0 );
	filter_sheet_ = tag->getOption<bool>( "filter_sheet", 0 );
	filter_loop_ = tag->getOption<bool>( "filter_loop", 0 );
	filter_helix_sheet_ = tag->getOption<bool>( "filter_helix_sheet", 1 );

	if (filter_helix_) {
		tr << "filter on "<< num_helix_ << " helics with length: "<<min_helix_length_<<"-"<< max_helix_length_ << std::endl;
		runtime_assert( num_helix_ > 0 );
	}

  if (filter_sheet_) {
    tr << "filter on "<< num_sheet_ << " sheet with length: "<<min_sheet_length_<<"-"<< max_sheet_length_ << std::endl;
    runtime_assert( num_sheet_ > 0 );
  }

  if (filter_loop_) {
    tr << "filter on "<< num_loop_ << " loop with length: "<<min_loop_length_<<"-"<< max_loop_length_ << std::endl;
    runtime_assert( num_loop_ > 0 );
  }

  if (filter_helix_sheet_) {
		tr << "filter on Sum of"<< num_helix_sheet_ << " helics with length: "<<min_helix_length_<<"-"<< max_helix_length_ << " AND sheet with length: "<<min_sheet_length_<<"-"<< max_sheet_length_ << std::endl;
		runtime_assert( num_helix_sheet_ > 0 );
	}

}

core::Size SecondaryStructureCountFilter::compute( core::pose::Pose const & pose ) const {
	core::pose::Pose pose_copy=pose;
	core::scoring::dssp::Dssp dssp( pose_copy );
	std::string dssp_ss=dssp.get_dssp_secstruct();
	tr << dssp_ss << std::endl;

	num_helix_pose_=0;
	num_sheet_pose_=0;
	num_loop_pose_=0;

	core::Size tmp_count=0;
	std::string::const_iterator iter=dssp_ss.begin();
	while (iter< dssp_ss.end()) {
		if (*iter=='H') {
				tmp_count=0;
				while (*iter=='H' && iter< dssp_ss.end() ) {
						++tmp_count;
						++iter;
				}

        if (tmp_count >= min_helix_length_ && tmp_count <= max_helix_length_) {
						num_helix_pose_+=1;
				}
				
		} else if ( *iter=='E') {
        tmp_count=0;
        while (*iter=='E' && iter< dssp_ss.end() ) {
            ++tmp_count;
            ++iter;
        }   
        
        if (tmp_count >= min_sheet_length_ && tmp_count <= max_sheet_length_) {
            num_sheet_pose_+=1;
        } 

		} else {
        tmp_count=0;
        while (*iter=='L' && iter< dssp_ss.end() ) {
            ++tmp_count;
            ++iter;
        }

        if (tmp_count >= min_loop_length_ && tmp_count <= max_loop_length_) {
            num_loop_pose_+=1;
        }
				++iter;
		}
	}

	tr << " Pose has " << num_helix_pose_ << " helix, " << num_sheet_pose_  << " sheet, " << num_helix_pose_ << " loop, according to dssp_reduced definition" <<  std::endl;
	return 0;
}

core::Real SecondaryStructureCountFilter::report_sm( core::pose::Pose const & pose ) const {
	compute( pose );
	return 0;
}

void SecondaryStructureCountFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	compute( pose );
	out << " Pose has " << num_helix_pose_ << " helix " << num_sheet_pose_  << " sheet " << num_helix_pose_ << " loop according to dssp_reduced definition" <<  std::endl;
}

protocols::filters::FilterOP
SecondaryStructureCountFilterCreator::create_filter() const { return new SecondaryStructureCountFilter; }

std::string
SecondaryStructureCountFilterCreator::keyname() const { return "SecondaryStructureCount"; }


} // filters
} // fldsgn
} // protocols
