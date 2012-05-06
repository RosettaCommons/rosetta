// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_core_scoring_symmetry_NBListCache_hh
#define INCLUDED_core_scoring_symmetry_NBListCache_hh

//  Unit headers
#include <core/scoring/symmetry/NBListCache.fwd.hh>

// Package Headers
#include <basic/datacache/CacheableData.hh>
//#include <basic/DataCache.hh>

#include <core/scoring/NeighborList.hh>

namespace core {
namespace scoring {
namespace symmetry {

class NBListCache : public basic::datacache::CacheableData {

public:

	NBListCache()
	: CacheableData()
	{};

	NBListCache( NeighborListOP nblist )
	: CacheableData()
	{
		nblist_ = nblist;
	};

	~NBListCache(){};

	basic::datacache::CacheableDataOP
	clone() const
	{
		return new NBListCache( *this );
	}

	NeighborListOP
	get_nblist()
	{
		return nblist_;
	}

	void
    set_nblist(
		NeighborListOP nblist
	)
    {
        nblist_ = nblist;
    }

private:

	NeighborListOP nblist_;

};

} // namespace symmetry
} // namespace scoring
} // namespace core
#endif
