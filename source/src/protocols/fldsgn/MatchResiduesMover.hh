// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fldsgn/MatchResiduesMover.hh
/// @brief  header file for MatchResiduesMover class
//  @brief  This mover returns the RMSD between a subset of residues of the movered pose against a list of residues in the reference pose.
/// @author Javier Castellanos ( javiercv@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_MatchResiduesMover_hh
#define INCLUDED_protocols_fldsgn_MatchResiduesMover_hh

// Unit Headers
#include <protocols/fldsgn/MatchResidues.hh>
#include <protocols/fldsgn/MatchResiduesMover.fwd.hh>

// Package Headers
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// Boost headers
#include <boost/tuple/tuple.hpp>


namespace protocols {
namespace fldsgn {

class MatchResiduesMover : public protocols::fldsgn::MatchResidues, public protocols::moves::Mover {
public:

	typedef protocols::moves::Mover Mover;
	typedef protocols::moves::MoverOP MoverOP;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;

	typedef utility::tag::TagCOP TagCOP;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::filters::Filters_map Filters_map;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor

	// @brief default constructor
	MatchResiduesMover();


	~MatchResiduesMover() override;

	// @brief make clone
	MoverOP clone() const override;

	// @brief make fresh instance
	MoverOP fresh_instance() const override;

	// @brief get name of this mover
	// XRW TEMP  std::string get_name() const override { return "MatchResiduesMover"; }


	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		Movers_map const & movers,
		Pose const & pose) override;

	void  apply( core::pose::Pose & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	bool superimpose_;


};

} // fldsgn
} // protocols

#endif
