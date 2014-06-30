// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/SSamountFilter.hh
/// @brief Calculate SS fraction, doesn't assume that poses are statick unmutable objects
/// @author Daniel Silva (dadriano@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_SSamountFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_SSamountFilter_hh

//Include Rosetta Core Stuff
#include <core/pose/Pose.hh>
#include <core/types.hh>

//Include Rosetta protocols 
#include <protocols/filters/Filter.hh>

//Include ObjexxFCL
#include <ObjexxFCL/FArray.all.hh>

//Include Rosetta utilities
#include <utility/vector1.hh>

//Include Rosetta XML tag reader
#include <utility/tag/Tag.fwd.hh>


namespace protocols {
namespace protein_interface_design {
namespace filters {

class SSamountFilter : public protocols::filters::Filter
{
public:
	SSamountFilter();
	
	SSamountFilter(
					core::Real const upper_threshold,
					core::Real const lower_threshold,
					core::Size target_chain,
					bool b_target_chain,
					bool b_discard_lonely_SS);
	bool apply( core::pose::Pose const & pose ) const;
	protocols::filters::FilterOP clone() const;
	protocols::filters::FilterOP fresh_instance() const{
		return new SSamountFilter();
	}
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~SSamountFilter();
	void parse_my_tag( utility::tag::TagCOP tag, 
							basic::datacache::DataMap & data_map, 
							protocols::filters::Filters_map const &, 
							protocols::moves::Movers_map const &, 
							core::pose::Pose const & reference_pose );
private:
	core::Real upper_threshold_;
	core::Real lower_threshold_;
	core::Size target_chain_;
	bool b_target_chain_;
	bool b_discard_lonely_SS_;
	core::pose::PoseOP reference_pose_;
};

} // filters
} // protein_interface_design
} // protocols

#endif //INCLUDED_protocols_protein_interface_design_filters_SSamountFilter_HH_
