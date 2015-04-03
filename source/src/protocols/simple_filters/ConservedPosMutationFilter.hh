// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/simple_filters/ConservedPosMutationFilter.hh
/// @brief header file for ConservedPosMutationFitler class.
/// @details
/// @author Florian Richter (floric@u.washington.edu), may 2011


#ifndef INCLUDED_protocols_simple_filters_ConservedPosMutationFilter_hh
#define INCLUDED_protocols_simple_filters_ConservedPosMutationFilter_hh

// Unit Headers
#include <protocols/simple_filters/ConservedPosMutationFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

// package headers
#include <protocols/toolbox/task_operations/SeqprofConsensusOperation.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers

// Utility headers

//// C++ headers

namespace protocols {
namespace simple_filters {


class ConservedPosMutationFilter : public filters::Filter {

public:
	typedef filters::Filter parent;
	typedef filters::Filter Filter;
	typedef filters::FilterOP FilterOP;

public:

	ConservedPosMutationFilter();

	~ConservedPosMutationFilter();

	FilterOP clone() const;

	FilterOP fresh_instance() const;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &  );


	/// @brief Returns true if the given pose passes the filter, false otherwise.
	/// In this case, a pose passes if it less than max_allowed_conserved_pos_mutations_
	/// mutations at conserved position. the decision whether a given position
	/// counts as conserved is made by the RestrictConservedLowDdgOperation
	/// task operation to prevent duplication of code
	virtual
	bool apply( core::pose::Pose const & pose ) const;

	void set_max_allowed_conserved_pos_mutations( core::Size value ){
		max_allowed_conserved_pos_mutations_ = value;
	}

	virtual std::string name() const {
		return "ConservedPosMutationFilter";
	}

private:

	toolbox::task_operations::RestrictConservedLowDdgOperationOP conserved_pos_taskop_;
	core::Size max_allowed_conserved_pos_mutations_;

};

} // filters
} // protocols

#endif
