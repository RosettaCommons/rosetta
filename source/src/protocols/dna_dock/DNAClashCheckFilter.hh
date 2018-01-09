// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/dna_dock/DNAClashCheckFilter.hh
/// @brief  header file for DNAClashCheckFilter class
/// @author Carl walkey (cwalkey@u.washington.edu)


#ifndef INCLUDED_protocols_dna_dock_DNAClashCheckFilter_hh
#define INCLUDED_protocols_dna_dock_DNAClashCheckFilter_hh

// Unit Headers
#include <protocols/dna_dock/DNAClashCheckFilter.fwd.hh>

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

class DNAClashCheckFilter : public protocols::filters::Filter {
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
	DNAClashCheckFilter();

	// @brief constructor with arguments
	DNAClashCheckFilter( core::Real const fa_rep_thresh, std::string dna_a_path, std::string dna_b_path );

	// @brief copy constructor
	DNAClashCheckFilter( DNAClashCheckFilter const & rval );

	~DNAClashCheckFilter() override;


public:// virtual constructor


	// @brief make clone
	protocols::filters::FilterOP clone() const override;

	// @brief make fresh instance
	protocols::filters::FilterOP fresh_instance() const override;


public:// accessor

	// @brief get name of this filter
	std::string name() const override { return "DNAClashCheck"; }

public:// setters
	void set_fa_rep_thresh( core::Real fa_rep_thresh );
	void set_dna_a_path( std::string dna_a_path );
	void set_dna_b_path( std::string dna_b_path );

public:// getters
	core::Real get_fa_rep_thresh() const;
	std::string get_dna_a_path() const;
	std::string get_dna_b_path() const;

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
	std::string dna_a_path_;
	std::string dna_b_path_;

};

} // dna_dock
} // protocols

#endif
