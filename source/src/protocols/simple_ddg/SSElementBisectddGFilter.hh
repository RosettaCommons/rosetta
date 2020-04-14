// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_ddg/SSElementBisectddGFilter.hh
/// @brief  Disconnects a protein by secondary structure elements and calculates the DDG between the elements
/// @author TJ Brunette (tjbrunette@gmail.com)


#ifndef INCLUDED_protocols_simple_ddg_SSElementBisectddGFilter_hh
#define INCLUDED_protocols_simple_ddg_SSElementBisectddGFilter_hh

// Unit Headers
#include <protocols/simple_ddg/SSElementBisectddGFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
// Utility headers

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>
#include <set>

//// C++ headers

namespace protocols {
namespace simple_ddg {

class SSElementBisectddGFilter : public protocols::filters::Filter{
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef std::string String;

	typedef utility::tag::TagCOP TagCOP;
	typedef basic::datacache::DataMap DataMap;


public:// constructor/destructor

	// @brief constructor w variables
	SSElementBisectddGFilter(core::scoring::ScoreFunctionOP scorefxn,Real threshold,bool report_avg,core::Size ignore_terminal_SS,bool only_n_term,bool only_c_term,bool skip_ss_element,bool report_sasa_instead,bool convert_charged_res_to_ala,protocols::moves::MoverOP relax_mover);

	// @brief default constructor
	SSElementBisectddGFilter();

	// @brief copy constructor
	SSElementBisectddGFilter( SSElementBisectddGFilter const & rval );


public:


	// @brief make clone
	filters::FilterOP clone() const override { return utility::pointer::make_shared< SSElementBisectddGFilter >(*this);}
	// @brief make fresh instance
	filters::FilterOP fresh_instance() const override { return utility::pointer::make_shared< SSElementBisectddGFilter >();}


public:// accessor


	// @brief get name of this filter


public:// virtual main operation


	Real report_sm(const Pose & pose ) const override;
	void report( std::ostream & out,const Pose & pose ) const override;
	protocols::loops::Loops get_ss_elements(const Pose & pose) const;
	Real get_ddg_bisect_score(core::Size element,protocols::loops::Loops ssElements, const Pose & pose) const;
	Real compute( const Pose & orig_pose ) const;
	bool apply(const Pose & pose ) const override;


public:// parser

	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	core::scoring::ScoreFunctionOP scorefxn_; //dflt NULL
	core::Size ignore_terminal_SS_;
	bool only_n_term_;
	bool only_c_term_;
	Real threshold_;
	bool report_avg_;
	bool skip_ss_element_;
	bool report_sasa_instead_;
	bool convert_charged_res_to_ala_;
	protocols::moves::MoverOP relax_mover_; //dflt NULL; in the unbound state, prior to taking the energy, should we do any relaxation

};

} // filters
} // protocols

#endif
