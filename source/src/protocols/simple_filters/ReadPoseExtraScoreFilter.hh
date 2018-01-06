// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ReadPoseExtraScoreFilter.hh
/// @brief definition of filter class ReadPoseExtraScoreFilter.
/// @author Jack Maguire, jack@med.unc.edu

#ifndef INCLUDED_protocols_simple_filters_ReadPoseExtraScoreFilter_hh
#define INCLUDED_protocols_simple_filters_ReadPoseExtraScoreFilter_hh

//unit headers
#include <protocols/simple_filters/ReadPoseExtraScoreFilter.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>


namespace protocols {
namespace simple_filters {

class ReadPoseExtraScoreFilter : public filters::Filter
{
public:
	/// @brief Constructor
	ReadPoseExtraScoreFilter();

	/// @brief Copy constructor
	ReadPoseExtraScoreFilter( ReadPoseExtraScoreFilter const & src );

	/// @brief Constructor with parameters
	ReadPoseExtraScoreFilter( std::string term, core::Real threshold );

	~ReadPoseExtraScoreFilter() override;

	/// @brief Returns false if the score term is greater than the threshold
	bool apply( core::pose::Pose const & pose ) const override;

	filters::FilterOP clone() const override {
		return filters::FilterOP( new ReadPoseExtraScoreFilter( *this ) );
	}

	filters::FilterOP fresh_instance() const override {
		return filters::FilterOP( new ReadPoseExtraScoreFilter() );
	}

	/// @brief Sets the name of the term being searched for in the pose
	inline void set_term_name( std::string term_name ) {
		term_name_ = std::move( term_name );
	}

	/// @brief apply() returns false if the score is greater than this threshold
	inline void set_threshold( core::Real const & threshold ) {
		threshold_ = threshold;
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;

	/// @brief polymorphic way to call compute()
	core::Real report_sm( core::pose::Pose const & pose ) const override;

	/// @brief attempt to extract the score term from the pose. This method does no calculation.
	core::Real compute( core::pose::Pose const &pose ) const;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	std::string
	name() const override {
		//return "ReadPoseExtraScoreFilter";
		return class_name();
	}

	static inline
	std::string
	class_name() {
		return "ReadPoseExtraScoreFilter";
	}

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string term_name_;
	core::Real threshold_;
};

} // simple_filters
} // protocols

#endif
