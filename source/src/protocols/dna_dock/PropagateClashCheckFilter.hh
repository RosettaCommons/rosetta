// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/dna_dock/PropagateClashCheckFilter.hh
/// @brief  header file for PropagateClashCheckFilter class
/// @author Carl walkey (cwalkey@u.washington.edu)


#ifndef INCLUDED_protocols_dna_dock_PropagateClashCheckFilter_hh
#define INCLUDED_protocols_dna_dock_PropagateClashCheckFilter_hh

// Unit Headers
#include <protocols/dna_dock/PropagateClashCheckFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

//// C++ headers

namespace protocols {
namespace dna_dock {

class PropagateClashCheckFilter : public protocols::filters::Filter {
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	// @brief default constructor
	PropagateClashCheckFilter();

	// @brief constructor with arguments
	PropagateClashCheckFilter( core::Real const fa_rep_thresh,
		std::string default_bridge_type,
		core::Real omega,
		core::Real rise,
		core::Size num_repeats,
		core::Size prop_dir );

	// @brief copy constructor
	PropagateClashCheckFilter( PropagateClashCheckFilter const & rval );

	~PropagateClashCheckFilter() override;


public:// virtual constructor


	// @brief make clone
	protocols::filters::FilterOP clone() const override;

	// @brief make fresh instance
	protocols::filters::FilterOP fresh_instance() const override;


public:// accessor

	// @brief get name of this filter
	std::string name() const override { return "PropagateClashCheck"; }

public:// setters
	void set_fa_rep_thresh( core::Real fa_rep_thresh );
	void set_default_bridge_type( std::string default_bridge_type );
	void set_omega( core::Real omega );
	void set_rise( core::Real rise );
	void set_num_repeats( core::Size num_repeats );

public:// getters
	core::Real get_fa_rep_thresh() const;
	std::string get_default_bridge_type() const;
	core::Real get_omega() const;
	core::Real get_rise() const;
	core::Size get_num_repeats() const;

public:// parser

	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		Movers_map const &,
		Pose const & ) override;


public:// virtual main operation

	// @brief returns true if the given pose passes the filter, false otherwise.
	bool apply( core::pose::Pose const & pose ) const override;

private:

	core::Real fa_rep_thresh_;
	std::string default_bridge_type_;
	core::Real omega_;
	core::Real rise_;
	core::Size num_repeats_;
	core::Size prop_dir_;

};

} // dna_dock
} // protocols

#endif
