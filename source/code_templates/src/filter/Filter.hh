// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.hh
/// @brief --brief--
/// @author --name-- (--email--)

#ifndef INCLUDED_--path_underscore--_--class--_hh
#define INCLUDED_--path_underscore--_--class--_hh

// Unit headers
#include <--path--/--class--.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

--namespace--

///@brief --brief--
class --class-- : public protocols::filters::Filter {

public:
	--class--();

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~--class--();

	/// @brief returns true if the structure passes the filter, false otherwise
	virtual bool
	apply( core::pose::Pose const & pose ) const;

	/// @brief required for reporting score values
	virtual core::Real
	report_sm( core::pose::Pose const & pose ) const;

	/// @brief allows printing data to a stream
	virtual void
	report( std::ostream & os, core::pose::Pose const & pose ) const;

public:
	virtual std::string
	get_name() const;

	/// @brief parse XML tag (to use this Filter in Rosetta Scripts)
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::filters::FilterOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::filters::FilterOP
	clone() const;

private:

};

--end_namespace--

#endif //INCLUDED_--path_underscore--_--class--_hh
