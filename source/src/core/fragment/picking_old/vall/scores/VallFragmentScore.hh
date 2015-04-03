// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/vall/scores/VallFragmentScore.hh
/// @brief  the base Vall FragmentScore struct
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_scores_VallFragmentScore_hh
#define INCLUDED_core_fragment_picking_old_vall_scores_VallFragmentScore_hh

// unit headers
#include <core/fragment/picking_old/vall/scores/VallFragmentScore.fwd.hh>

// type headers
#include <core/types.hh>

// package headers
#include <core/fragment/picking_old/vall/VallSection.hh>

// C++ headers
#include <ostream>
#include <sstream>

#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace scores {


/// @brief  the base Vall FragmentScore struct
struct VallFragmentScore {


public: // typedefs


	typedef core::Real Real;


public: // concept typedefs


	/// @brief typedef for Bookmark concept
	typedef VallSection::PageConstIterator PageConstIterator;


public: // concept translation typedefs


	typedef PageConstIterator VallResidueConstIterator;


public: // friends


	/// @brief friend declaration for stream operator << using virtual friend idiom
	friend std::ostream & operator <<( std::ostream & o, VallFragmentScore const & s );


public: // construct/destruct


	/// @brief default constructor
	inline
	VallFragmentScore() :
		score( 0.0 )
	{}


	/// @brief copy constructor
	inline
	VallFragmentScore( VallFragmentScore const & rval ) :
		extent_begin( rval.extent_begin ),
		extent_end( rval.extent_end ),
		score( rval.score )
	{}


	/// @brief default destructor
	inline
	virtual
	~VallFragmentScore() {}


public: // assignment


	/// @brief copy assignment
	inline
	VallFragmentScore & operator =( VallFragmentScore const & rval ) {
		if ( this != &rval ) {
			extent_begin = rval.extent_begin;
			extent_end = rval.extent_end;
			score = rval.score;
		}
		return *this;
	}


public: // comparator


	/// @brief '<' comparison
	inline
	bool operator <( VallFragmentScore const & rval ) const {
		return score < rval.score;
	}


public: // status


	/// @brief return string describing contents or status
	/// @remarks overriding this will automatically give useable
	///  stream output due to friend operator << in base class
	virtual
	std::string to_string() const {
		std::ostringstream s;
		s << "score = " << score;
		return s.str();
	}


public: // convenience


	/// @brief compute distance (effectively the length of the extent)
	///  from begin -> end
	inline
	Size distance() const {
		return std::distance( extent_begin, extent_end );
	}


public: // data


	/// @brief points to the beginning of the fragment
	VallResidueConstIterator extent_begin;


	/// @brief points just beyond the end of the fragment
	VallResidueConstIterator extent_end;


	/// @brief the cumulative score
	Real score;


};


/// @brief stream output <<
/// @remarks uses virtual friend idiom, calls to_string()
inline
std::ostream & operator <<( std::ostream & o, VallFragmentScore const & s ) {
	o << s.to_string();
	return o;
}


} // scores
} // vall
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_vall_scores_VallFragmentScore_HH */
