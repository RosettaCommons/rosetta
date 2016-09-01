// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/moves/FilterReportAsPoseExtraScoresMover.hh
/// @brief This Mover runs a Filter and dumps the report_sm value to Pose's extra scores (for later JD reporting)
/// @author Steven Lewis (smlewi@gmail.com)

#ifndef INCLUDED_protocols_moves_FilterReportAsPoseExtraScoresMover_hh
#define INCLUDED_protocols_moves_FilterReportAsPoseExtraScoresMover_hh

// Unit headers
#include <protocols/moves/FilterReportAsPoseExtraScoresMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

#include <string>

namespace protocols {
namespace moves {

///@brief This Mover runs a Filter and dumps the report_sm value to Pose's extra scores (for later JD reporting)
class FilterReportAsPoseExtraScoresMover : public protocols::moves::Mover {

public:

	FilterReportAsPoseExtraScoresMover();

	/// @details with args: takes filter and the string to report by
	FilterReportAsPoseExtraScoresMover( protocols::filters::FilterOP filter, std::string report_as );

	// copy constructor (not needed unless you need deep copies)
	//NOT defined, which means filter_ will do shallow copies; change that if you need it.
	//FilterReportAsPoseExtraScoresMover( FilterReportAsPoseExtraScoresMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~FilterReportAsPoseExtraScoresMover();

	static std::string
	class_name();

public:
	// mover virtual API
	virtual void
	apply( core::pose::Pose & pose );

	virtual void
	show( std::ostream & output = std::cout ) const;

	virtual std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//not defined, which means filter_ will do shallow copies, change that if you need it.
	//FilterReportAsPoseExtraScoresMover & operator=( FilterReportAsPoseExtraScoresMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;

	/// @brief set filter we report from
	void set_filter( protocols::filters::FilterOP filter );
	void set_report_as( std::string const & report_as );

	protocols::filters::FilterOP get_filter() const { return filter_; }
	std::string const & get_report_as() const { return report_as_; }


private:

	protocols::filters::FilterOP filter_;
	std::string report_as_;

};

std::ostream &
operator<<( std::ostream & os, FilterReportAsPoseExtraScoresMover const & mover );

} //protocols
} //moves

#endif //protocols/moves_FilterReportAsPoseExtraScoresMover_hh
