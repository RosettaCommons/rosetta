// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   CanonicalFragmentHistory.hh
/// @brief
/// @author TJ Brunette


#ifndef INCLUDED_protocols_simple_moves_VallLookbackData_hh
#define INCLUDED_protocols_simple_moves_VallLookbackData_hh

#include <core/scoring/methods/vall_lookback/VallLookbackData.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/CacheableData.hh>
#include <core/types.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace methods {

///////////////////////////////////////////////////////////////////////////////
// ncs residue mapping
//   - stored in the pose and used by other movers (fragment insertion for example)
class VallLookbackData:  public basic::datacache::CacheableData {
public:
	VallLookbackData( core::pose::Pose &pose );

	basic::datacache::CacheableDataOP  clone() const{
		return basic::datacache::CacheableDataOP(new VallLookbackData(*this));
	}

	void set_rmsd( core::Size resid, core::Real rmsd);
	core::Real get_rmsd( core::Size resid ) const;
	void set_res_changed(core::Size resid, bool changed);
	bool get_res_changed(core::Size resid) const;
	core::Size size() { return rmsd_history_.size(); }

private:
	utility::vector1< core::Real > rmsd_history_;  //residue rmsd
	utility::vector1< bool > res_changed_history_; //residue changed
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	VallLookbackData();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


}//end methods
}//end scoring
}//end core
#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_methods_vall_lookback_VallLookbackData )
#endif // SERIALIZATION


#endif

