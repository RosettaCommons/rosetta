// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // symmetry
} // conformation
} // core
#endif
