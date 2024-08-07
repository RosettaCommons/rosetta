// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/Worst9merFilter.hh
/// @brief  Find the rmsd to the worst 9mer in the chain
/// @author TJ Brunette (tjbrunette@gmail.com)


#ifndef INCLUDED_protocols_simple_filters_Worst9merFilter_hh
#define INCLUDED_protocols_simple_filters_Worst9merFilter_hh

// Unit Headers
#include <protocols/simple_filters/Worst9merFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <protocols/indexed_structure_store/SSHashedFragmentStore.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>


//// C++ headers

namespace protocols {
namespace simple_filters {

class Worst9merFilter : public protocols::filters::Filter{
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


	// @brief default constructor
	Worst9merFilter();


	~Worst9merFilter() override;


public:// virtual constructor

	// @brief make clone
	filters::FilterOP clone() const override { return filters::FilterOP(new Worst9merFilter(*this));}
	// @brief make fresh instance
	filters::FilterOP fresh_instance() const override { return filters::FilterOP(new Worst9merFilter());}

	// @brief get name of this filter
	Real report_sm(const Pose & pose ) const override;
	void report( std::ostream & out,const Pose & pose ) const override;
	Real compute( const Pose & pose ) const;
	bool apply(const Pose & pose ) const override;
	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &
	) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	void
	set_report_mean_median( bool report ) { report_mean_median_ = report; }

	void
	set_residue_selector( core::select::residue_selector::ResidueSelector const & selector );

private:

	void
	write_mean_median( Pose const & pose, core::Real mean, core::Real median ) const;

private:

	Real filtered_value_;
	Real rmsd_lookup_thresh_;
	protocols::indexed_structure_store::SSHashedFragmentStoreOP SSHashedFragmentStoreOP_;
	bool ignore_terminal_res_;
	bool only_helices_;
	bool report_mean_median_;
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;
	std::string fragment_store_path_;
	std::string fragment_store_format_;
	std::string fragment_store_compression_;

};

} // filters
} // protocols

#endif

