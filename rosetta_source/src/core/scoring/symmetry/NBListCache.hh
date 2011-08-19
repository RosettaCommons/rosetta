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
