// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/LoopAnalyzerFilter.hh
/// @brief LoopAnalyzerFilter examines loop structures and packages extra scores into a Job object; you can filter on its "LAM score".  This was originally LoopAnalyzerMover (which treated its Pose as const anyway), but has been converted to a filter.
/// @author Steven Lewis (smlewi@gmail.com)

#ifndef INCLUDED_protocols_loops_filters_LoopAnalyzerFilter_hh
#define INCLUDED_protocols_loops_filters_LoopAnalyzerFilter_hh

// Unit headers
#include <protocols/analysis/LoopAnalyzerFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

#include <protocols/loops/Loops.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace loops {
namespace filters {

///@brief LoopAnalyzerFilter examines loop structures and packages extra scores into a Job object; maybe you can filter on its "LAM score" (not yet implemented)
class LoopAnalyzerFilter : public protocols::filters::Filter {

public:
	LoopAnalyzerFilter();

	/// @brief loops is the loops object to say what to analyze.  tracer controls whether output is dumped to a tracer or a Job object (for output with the Pose).
	LoopAnalyzerFilter( protocols::loops::Loops const & loops, bool const tracer = false );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~LoopAnalyzerFilter();

	LoopAnalyzerFilter( LoopAnalyzerFilter const & rhs );

	/// @brief returns true if the structure passes the filter, false otherwise
	bool
	apply( core::pose::Pose const & pose ) const override;

	/// @brief required for reporting score values
	core::Real
	report_sm( core::pose::Pose const & pose ) const override;

	/// @brief allows printing data to a stream
	void
	report( std::ostream & os, core::pose::Pose const & pose ) const override;

public:
	std::string
	get_name() const;

	/// @brief parse XML tag (to use this Filter in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	clone() const override;

public: ///////////////////getters, setters/////////////
	/// @brief set loops object, because public setters/getters are a rule
	void set_loops( protocols::loops::LoopsCOP loops );

	/// @brief get loops object, because public setters/getters are a rule
	protocols::loops::LoopsCOP const & get_loops( void ) const ;

	/// @brief set tracer bool, because public setters/getters are a rule
	inline void set_use_tracer( bool tracer ) { tracer_ = tracer; }

	/// @brief get tracer bool, because public setters/getters are a rule
	inline bool get_use_tracer( void ) const { return tracer_; }

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	/// @brief used to store a copy of the input loops
	protocols::loops::LoopsCOP loops_;

	/// @brief output to tracer or PDB/silent file
	bool tracer_;

};

} //protocols
} //loops
} //filters

#endif //INCLUDED_protocols_loops_filters_LoopAnalyzerFilter_hh
