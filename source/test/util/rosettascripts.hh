// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/util/rosettascripts.hh
/// @brief  helper functions for unit testing RosettaScripts functionality.
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Note: I don't claim this is a complete toolkit for unittesting RosettaScripts components.
// Feel free to add anything here that makes it simpler for you to add Unit tests for RosettaScripts functionality

// // Sample parse_my_tag() test
//
// DataMap data;
// Filters_map filters;
// Movers_map movers;
//
// prime_Movers( movers ); // Adds "null" to movers map - optional
// prime_Filters( filters ); // Adds true_filter, false_filter - optional
// prime_Data( data ); // Adds score12, commandline scorefunctions - optional
//
// filters["stubby"] = new StubFilter( true, -3.14 );
//
// ScoreFunctionOP rep_only = ScoreFunctionFactory::create_score_function( "fa_rep_only" );
// data.add( "scorefxns", "rep", rep_only );
//
// MyMover testmover;
// TagPtr tag = tagptr_from_string("<MyTest name=test filter=stubby mover=null reps=4>\n"
//		"</MyTest>"); // Remember that C++ has implicit string literal concatenation, but note that the \n is required for the tag parser
// testmover.parse_my_tag( tag, data, filters, movers, pose );
//
//
// TS_ASSERT_EQUALS( testmover.reps(), 4 )
// //... etc.

#ifndef INCLUDED_util_rosettascripts_HH
#define INCLUDED_util_rosettascripts_HH

#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/moves/NullMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/types.hh>

#include <utility/tag/Tag.hh>

#include <string>
#include <sstream>

//Bring out items to ease typing
using protocols::filters::Filter;
using protocols::filters::FilterOP;
using protocols::moves::Mover;
using protocols::moves::MoverOP;

using utility::tag::TagPtr;
using protocols::moves::DataMap;
using protocols::moves::Movers_map; // A std::map of string to MoverOP
using protocols::filters::Filters_map; // A std::map string to FilterOP

using protocols::filters::TrueFilter;
using protocols::filters::FalseFilter;
using protocols::moves::NullMover;

///@brief Generate a tagptr from a string
/// For parse_my_tag tests, only do the relevant tag, not the full <ROSETTASCRIPTS> ... </ROSETTASCRIPTS> wrapped tag.
inline TagPtr tagptr_from_string(std::string input) {
	std::stringstream instream( input );
	return utility::tag::Tag::create( instream );
}

///@brief setup filters map with some of the the RosettaScript defaults
inline void prime_Filters( Filters_map & filters ) {
	filters["true_filter"] = new protocols::filters::TrueFilter;
	filters["false_filter"] = new protocols::filters::FalseFilter;
}

///@brief setup movers map with some of the the RosettaScript defaults
inline void prime_Movers( Movers_map & movers ) {
	movers["null"] = new protocols::moves::NullMover;
}

///@brief setup data map with *some* of the the RosettaScript defaults

inline void prime_Data( DataMap & data ) {
	core::scoring::ScoreFunctionOP commandline_sfxn = core::scoring::getScoreFunction();
	core::scoring::ScoreFunctionOP talaris2013 = core::scoring::getScoreFunction();

	data.add( "scorefxns", "commandline", commandline_sfxn );
	data.add( "scorefxns", "talaris2013", talaris2013 );
}

///@brief A simple filter for helping to test nested classes
/// will apply() with the given truth value,
/// report_sm() with the given value,
/// and report() with the truth,

class StubFilter : public Filter {
public:
	StubFilter( bool truth = true, core::Real value = 0, std::string tag = "" ) :
		truth_(truth),
		value_(value),
		tag_(tag)
	{}
	FilterOP clone() const { return new StubFilter( *this ); }
	FilterOP fresh_instance() const {	return new StubFilter(); }
	void set( bool truth, core::Real value) { truth_ = truth; value_ = value; }
	bool apply( core::pose::Pose const & ) const { return truth_; }
	core::Real report_sm( core::pose::Pose const & ) const { return value_;}
	void report( std::ostream & ostream, core::pose::Pose const & ) const {
		ostream << "StubFilter " << tag_ << ": " << truth_ << " " << value_ << std::endl;
	}
public: // Yes, public - deliberately so people can easily change them, if they want to
	bool truth_;
	core::Real value_;
	std::string tag_;
};

typedef utility::pointer::owning_ptr< StubFilter >  StubFilterOP;
typedef utility::pointer::owning_ptr< StubFilter const >  StubFilterCOP;

///@brief A simple filter for helping to test nested classes
/// will apply() with the given truth value,
/// When called, report_sm() will cycle through the given list of values

class StubMultiFilter : public Filter {
public:
	StubMultiFilter( bool truth = true,std::string tag = "" ) :
		truth_(truth),
		pos_(1),
		tag_(tag)
	{}
	FilterOP clone() const { return new StubMultiFilter( *this ); }
	FilterOP fresh_instance() const {	return new StubMultiFilter(); }
	void set( utility::vector1<core::Real> const & values, bool truth=true, std::string tag="") { values_ = values, truth_ = truth; tag_ = tag; }
	void push_back( core::Real value ) { values_.push_back( value); }
	void set_pos( core::Size pos = 1) { pos_ = pos; }
	bool apply( core::pose::Pose const & ) const { return truth_; }
	core::Real report_sm( core::pose::Pose const & ) const {
		if ( pos_ == 0 || pos_ > values_.size() ) { pos_ = 1; }
		++pos_; // Overflow taken care of next time around
		return values_[ pos_ - 1 ];
	}
	void report( std::ostream & ostream, core::pose::Pose const & ) const {
		ostream << "StubMultiFilter " << tag_ << ": " << truth_ << " ";
		for( core::Size ii(1); ii <= values_.size(); ++ii ) {
			ostream << values_[ ii ] << " ";
		}
		ostream << std::endl;
	}
public: // Yes, public - deliberately so people can easily change them, if they want to
	bool truth_;
	mutable core::Size pos_;
	std::string tag_;
	utility::vector1<core::Real> values_;
};

typedef utility::pointer::owning_ptr< StubMultiFilter >  StubMultiFilterOP;
typedef utility::pointer::owning_ptr< StubMultiFilter const >  StubMultiFilterCOP;

#endif
