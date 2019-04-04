// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/SecretionPredictionFilter.hh
/// @brief  Filter for greasy helices  https://www.nature.com/articles/nature06387
/// @author Designed by John Wang(jyjwang@uw.edu) converted to a filter by TJ Brunette(tjbrunette@gmail.com)


#ifndef INCLUDED_protocols_simple_filters_SecretionPredictionFilter_hh
#define INCLUDED_protocols_simple_filters_SecretionPredictionFilter_hh

// Unit Headers
#include <protocols/simple_filters/SecretionPredictionFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
// Utility headers

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>

//// C++ headers

namespace protocols {
namespace simple_filters {

class SecretionPredictionFilter : public protocols::filters::Filter{
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef std::string String;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;

	struct mutt
	{
		Size resnum; //residue number
		std::string aa; //amino acid it's mutated to
		Real ddG_ins; //change in dG_ins that the mutation causes
		Real dscore; //change in score that the mutation causes
		bool operator< (const mutt & rhs) const
		{
			return ddG_ins < rhs.ddG_ins; //to compare mutants by their ddG_ins
		}
	} ;

	struct tm_region //each region of interest within a protein
	{
		Size index; //where the region begins
		std::string sequence; //the sequence of that region
		core::Real dG_ins; //the dG_ins of that region
		bool operator< (const tm_region & a) const
		{
			return dG_ins < a.dG_ins; //to compare regions by their dG_ins
		}
		core::Real score;
		std::vector<SecretionPredictionFilter::mutt> possible_mutants; //each region of interest can have multiple possible mutants
	} ;


public:// constructor/destructor

	// @brief default constructor
	SecretionPredictionFilter();

	// @brief copy constructor
	SecretionPredictionFilter( SecretionPredictionFilter const & rval );

	virtual ~SecretionPredictionFilter();


public:// virtual constructor


	// @brief make clone
	filters::FilterOP clone() const override { return filters::FilterOP(utility::pointer::make_shared<SecretionPredictionFilter>(*this));}
	// @brief make fresh instance
	filters::FilterOP fresh_instance() const override { return filters::FilterOP(utility::pointer::make_shared<SecretionPredictionFilter>());}


public:// accessor


	// @brief get name of this filter


public:// virtual main operation
	double get_coeff(const std::string & str1, const std::string & str2) const;
	core::Real dG_ins_for_window( const std::string & window ) const;
	std::string get_window_at_index( const Pose & pose , const core::Size & wlen, const core::Size & start ) const;
	utility::vector1<SecretionPredictionFilter::tm_region> find_tm_regions( std::ostream & out, const Pose & pose ) const;


	Real report_sm(const Pose & pose ) const override;
	void report( std::ostream & out,const Pose & pose ) const override;
	Real compute( const Pose & pose ) const;
	bool apply(const Pose & pose ) const override;


public:// parser

	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		filters::Filters_map const &,
		Movers_map const &,
		Pose const & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	Real threshold_;
	Size window_size_;
	Real dG_ins_threshold_;
	std::string keep_;
};

} // filters
} // protocols

#endif
