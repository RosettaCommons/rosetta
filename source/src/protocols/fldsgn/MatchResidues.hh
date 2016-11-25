// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fldsgn/MatchResidues.hh
/// @brief  header file for MatchResidues class
//  @brief  returns the RMSD between a subset of residues of the movered pose against a list of residues in the reference pose.
/// @author Javier Castellanos ( javiercv@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_MatchResidues_hh
#define INCLUDED_protocols_fldsgn_MatchResidues_hh

// Unit Headers

// Package Headers
#include <protocols/moves/Mover.fwd.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// Boost headers
#include <boost/tuple/tuple.hpp>


namespace protocols {
namespace fldsgn {

class MatchResidues {
public:

	typedef protocols::moves::Mover Mover;
	typedef protocols::moves::MoverOP MoverOP;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;

	typedef utility::tag::TagCOP TagCOP;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;
	typedef protocols::filters::Filters_map Filters_map;
	typedef utility::vector1< Size > VecSize;
	typedef utility::vector1< VecSize > VecVecSize;


public:// constructor/destructor

	// @brief default constructor
	MatchResidues();

	virtual ~MatchResidues();


public:// parser

	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & );

	static void provide_attributes_and_subelements( utility::tag::AttributeList & attlist, utility::tag::XMLSchemaSimpleSubelementList & ssl );


public:
	core::Real compute( core::pose::Pose const & pose, VecSize & best_fit ) const;
	core::Real compute_comb( core::pose::Pose const & pose, VecSize const & comb ) const;
	core::Real superimpose_comb( core::pose::Pose & pose, VecSize const & comb ) const;
	core::Real threshold() const { return threshold_; }
	void threshold(core::Real value) {  threshold_ = value; }

private:

	std::map< std::string, boost::tuple<Size, Size> >  map_ss_segments( std::string const & ss) const;

	void cart_product( VecVecSize& rvvi, VecSize&  rvi, VecVecSize::const_iterator me, VecVecSize::const_iterator end ) const;
	VecVecSize cart_product( VecVecSize const & input) const;

private:
	core::pose::Pose reference_pose_;
	VecSize reference_residues_indexes_;
	VecVecSize mod_segment_prod_;

	core::Real threshold_;

};

} // fldsgn
} // protocols

#endif
