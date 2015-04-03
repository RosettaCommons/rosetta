// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/hbtrie/HBCPData.hh
/// @brief  This class enforces the sc/bb hbond exclusion rule for use in the trie 
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_hbonds_hbtrie_HBCPData_hh
#define INCLUDED_core_scoring_hbonds_hbtrie_HBCPData_hh

// Unit Headers
#include <core/scoring/hbonds/hbtrie/HBCPData.fwd.hh>

// Project Headers

// STL Headers
#include <iosfwd>

namespace core {
namespace scoring {
namespace hbonds {
namespace hbtrie {

class HBCPData
{
public:
	HBCPData();

	bool avoid_sc_hbonds() const { return avoid_sc_hbonds_; }
	bool is_sc() const { return is_sc_; }

	void avoid_sc_hbonds( bool setting ) { avoid_sc_hbonds_ = setting; }
	void is_sc( bool setting ) { is_sc_ = setting; }

	inline
	bool operator < ( HBCPData const & other ) const
	{
		return ( (avoid_sc_hbonds_ < other.avoid_sc_hbonds_) ||
			((avoid_sc_hbonds_ == other.avoid_sc_hbonds_) && (is_sc_ < other.is_sc_ )) );
	}

	inline
	bool operator == ( HBCPData const & other ) const
	{
		return avoid_sc_hbonds_ == other.avoid_sc_hbonds_ && is_sc_ == other.is_sc_;
	}

	void print( std::ostream & os ) const;

	//void
	//set_count_pair_data_to_use( Size connection_id ) const;


private:
	bool avoid_sc_hbonds_;
	bool is_sc_;


};

std::ostream & operator << ( std::ostream & os, HBCPData const & cpdat );

} // namespace hbtrie
} // namespace hbonds
} // namespace scoring
} // namespace core

#endif
