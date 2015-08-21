// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/filters/ParallelBetaPairingPreferenceFilter.hh
/// @brief header file for ParallelBetaPairingPreferenceFilter class.
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_filters_ParallelBetaPairingPreferenceFilter_hh
#define INCLUDED_protocols_fldsgn_filters_ParallelBetaPairingPreferenceFilter_hh

// Unit Headers
#include <protocols/fldsgn/filters/ParallelBetaPairingPreferenceFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


//// C++ headers

namespace protocols {
namespace fldsgn {
namespace filters {

class ParallelBetaPairingPreferenceFilter : public protocols::filters::Filter {
public: // typedef

	typedef protocols::filters::Filter Super;

	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::chemical::AA AA;
	typedef protocols::filters::Filter Filter;
	typedef std::string String;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::pose::Pose Pose;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	// @brief default constructor
	ParallelBetaPairingPreferenceFilter();

	// @brief constructor with arguments
	// Undefinded comminting out to fix PyRosetta build  ParallelBetaPairingPreferenceFilter( String const & ss );

	// @brief copy constructor
	ParallelBetaPairingPreferenceFilter( ParallelBetaPairingPreferenceFilter const & rval );

	// @brief destructor
	virtual ~ParallelBetaPairingPreferenceFilter(){}


public:// virtual constructor


	// @brief make clone
	virtual FilterOP clone() const { return FilterOP( new ParallelBetaPairingPreferenceFilter( *this ) ); }

	// @brief make fresh instance
	virtual FilterOP fresh_instance() const { return FilterOP( new ParallelBetaPairingPreferenceFilter() ); }


public:// set filter value


	void filter_value( Real const value );


public:// main calculator


	/// @brief compute number of contacts
	Real compute( Pose const & pose ) const;


public:// helper functions


	/// @brief
	Real score_pairmatrix( AA aa1, AA aa2 ) const;


public:// virtual main operations


	/// @brief used to report score
	virtual Real report_sm( Pose const & pose ) const;

	/// @brief used to report score
	virtual void report( std::ostream & out, Pose const & pose ) const;

	// @brief returns true if the given pose passes the filter, false otherwise.
	virtual bool apply( Pose const & pose ) const;


public:// parser


	virtual void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & );


private:


	/// @brief
	Real filter_value_;

	/// @brief
	utility::vector1< utility::vector1< Real > >  score_pairmatrix_;

	/// @brief
	bool verbose_;


};

} // filters
} // fldsgn
} // protocols

#endif
