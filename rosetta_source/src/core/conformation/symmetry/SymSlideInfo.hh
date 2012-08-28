// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief  SymSlideInfo data container
/// @file   core/conformation/symmetry/SymSlideInfo.hh
/// @author Ingemar Andre


#ifndef INCLUDED_core_conformation_symmetry_SymSlideInfo_hh
#define INCLUDED_core_conformation_symmetry_SymSlideInfo_hh

// Utility headers
#include <core/conformation/symmetry/SymSlideInfo.fwd.hh>
#include <core/types.hh>
// C++ headers
#include <string>
#include <vector>

namespace core {
namespace conformation {
namespace symmetry {

	enum SlideType {
		SEQUENTIAL = 1,
		ORDERED_SEQUENTIAL,
		RANDOM
	};

	enum SlideCriteriaType {
		CEN_DOCK_SCORE = 1,
		FA_REP_SCORE,
		CONTACTS,
		TOTAL_NUM_CRITERIA
	};

class SymSlideInfo {

	public:

	/// @brief constructor
	SymSlideInfo();

	/* SymSlideInfo(
		SlideType slide_type,
		SlideCriteriaType score_criteria,
		std::string SlideCriteriaVal = "AUTOMATIC",
		std::vector<core::Size> slide_order = std::vector<core::Size>()
	); */

	/// @brief copy constructor
	SymSlideInfo( SymSlideInfo const & src );

	SymSlideInfo &
  operator=( SymSlideInfo const & src );

	~SymSlideInfo();

	// setter functions
	void set_slide_type( SlideType slide_type );
	void set_SlideCriteriaType( SlideCriteriaType score_criteria );
	void set_SlideCriteriaVal( std::string SlideCriteriaVal );
	void set_slide_order( std::vector<core::Size> slide_order );

	// get functions
	SlideType get_slide_type() const;
	SlideCriteriaType get_SlideCriteriaType() const;
	std::string get_SlideCriteriaVal() const;
	std::vector<core::Size> get_slide_order() const;

	friend
	bool
	operator==(SymSlideInfo const & a, SymSlideInfo const & b);

	friend
	bool
	operator!=(SymSlideInfo const & a, SymSlideInfo const & b);


	private:
		SlideType slide_type_;
		SlideCriteriaType score_criteria_;
		std::string SlideCriteriaVal_;
		std::vector<core::Size> slide_order_;

};

} // symmetry
} // conformation
} // core
#endif
