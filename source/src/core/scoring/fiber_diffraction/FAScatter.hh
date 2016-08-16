// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Cache for scattering factors in all atom mode
/// @author Wojciech Potrzebowski and Ingemar Andre

#ifndef INCLUDED_core_scoring_fiber_diffraction_FAScatter_hh
#define INCLUDED_core_scoring_fiber_diffraction_FAScatter_hh

#include <core/scoring/fiber_diffraction/FAScatter.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace fiber_diffraction {

/// silly vector1 wrapper class so we can derive from PoseCachedData

class FAScatter : public basic::datacache::CacheableData {
public:

	FAScatter( utility::vector0< utility::vector1< utility::vector1< core::Real > > > & form_factors_in ):
		form_factors_( form_factors_in )
	{}


	basic::datacache::CacheableDataOP clone() const {
		return basic::datacache::CacheableDataOP( new FAScatter( *this ) );
	}

	utility::vector0< utility::vector1< utility::vector1< core::Real > > > getValues() {
		return form_factors_;
	}

private:
	utility::vector0< utility::vector1< utility::vector1< core::Real > > > form_factors_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	FAScatter();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

FAScatter &
retrieve_fa_scatter_from_pose( pose::Pose & pose );

} // namespace fiber_diffraction
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_fiber_diffraction_FAScatter )
#endif // SERIALIZATION


#endif
