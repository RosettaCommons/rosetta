// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/SequenceDistanceFilter.hh

#ifndef INCLUDED_protocols_simple_filters_SequenceDistanceFilter_hh
#define INCLUDED_protocols_simple_filters_SequenceDistanceFilter_hh

//unit headers
#include <protocols/simple_filters/SequenceDistanceFilter.fwd.hh>

// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace simple_filters {

/// @brief test whether a pose contains a comment that evaluates to a predefined value. This is useful in controlling execution flow in RosettaScripts.
class SequenceDistance : public filters::Filter
{
public:
	SequenceDistance();
	~SequenceDistance() override;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new SequenceDistance( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new SequenceDistance() );
	}

	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & pose) override;
	core::Size compute( core::pose::Pose const & pose ) const;

	std::string sequence_comment_id() const { return sequence_comment_id_; }
	void sequence_comment_id( std::string const & s ) { sequence_comment_id_ = s; }

	std::string target_seq() const { return target_seq_; }
	void target_seq( std::string const & s ) { target_seq_ = s; }

	std::string pose_seq() const { return pose_seq_; }
	void pose_seq( std::string const & s ) { pose_seq_ = s; }


	core::Size threshold() const { return threshold_; }
	void threshold( core::Size const & t ) { threshold_ = t; }

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string sequence_comment_id_; //dflt ""; define the comment name
	std::string target_seq_; // dflt ""; the sequence you want to compare to.
	std::string pose_seq_; // dflt ""; the sequence you want to compare to.
	core::Size threshold_; //dflt ""; define the comment value
};
}
}

#endif
