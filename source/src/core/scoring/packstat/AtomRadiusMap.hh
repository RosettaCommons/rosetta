// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packstat/AtomRadiusMap.hh
///
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_packstat_AtomRadiusMap_hh
#define INCLUDED_core_scoring_packstat_AtomRadiusMap_hh

#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/AtomRadiusMap.fwd.hh>


#include <map>
#include <string>

namespace core {
namespace scoring {
namespace packstat {


/// @brief
class AtomRadiusMap
{

	// friend std::istream & operator>> ( std::istream & in , AtomRadiusMap       & pdb  );
	// friend std::ostream & operator<< ( std::ostream & out, AtomRadiusMap const & pdb  );

public: // Creation

	AtomRadiusMap( int include_water = -1 ); // -1 (default) means read from global options
	~AtomRadiusMap() {}

	PackstatReal get_radius( std::string type, std::string res ) const;

private:

	std::map< std::pair<std::string const,std::string const>, PackstatReal > type_map_;

}; // AtomRadiusMap


} // namespace packstat
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_packstat_AtomRadiusMap_HH
