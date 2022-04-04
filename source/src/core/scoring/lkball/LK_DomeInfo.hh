// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/LK_BallInfo.hh
/// @brief  Orientation dependent variant of the LK Solvation using
/// @author Phil Bradley

#ifndef INCLUDED_core_scoring_methods_LK_DomeInfo_HH
#define INCLUDED_core_scoring_methods_LK_DomeInfo_HH


// Unit headers
#include <core/scoring/lkball/LK_DomeEnergy.fwd.hh>
#include <core/scoring/lkball/LK_BallInfo.hh>

// // Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/RestypeDestructionEvent.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <basic/datacache/CacheableData.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/fixedsizearray1.hh>
#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace lkball {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef utility::fixedsizearray1< Real, MAX_N_WATERS_PER_ATOM > WaterOcclusions;
typedef utility::fixedsizearray1< Real, MAX_N_WATERS_PER_ATOM > WaterSolValues;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details  Holds the locations of ideal waters attached to the atoms of a Residue
class LKD_ResidueInfo; // fwd
typedef utility::pointer::shared_ptr< LKD_ResidueInfo > LKD_ResidueInfoOP;

class LKD_ResidueInfo : public LKB_ResidueInfo {
public:
	~LKD_ResidueInfo() override;

public:

	LKD_ResidueInfo( conformation::Residue const & rsd, LK_DomeEnergy const * lk_dome );

	LKD_ResidueInfo( LKD_ResidueInfo const & src );

	LKD_ResidueInfo &
	operator=(LKD_ResidueInfo const & src);

	LKD_ResidueInfo();

	basic::datacache::CacheableDataOP
	clone() const override;

	utility::vector1< WaterOcclusions > &
	water_occlusions();

	utility::vector1< WaterOcclusions > const &
	water_occlusions() const;

	utility::vector1< Real > const &
	water_sol_values() const;

	void clear_occlusions();



private:
	utility::vector1< WaterOcclusions > water_occlusions_;
	utility::vector1< Real > water_sol_values_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

typedef utility::pointer::shared_ptr< LKD_ResidueInfo > LKD_ResidueInfoOP;
typedef utility::pointer::shared_ptr< const LKD_ResidueInfo > LKD_ResidueInfoCOP;



}
}
}
#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_lkball_LK_DomeInfo )
#endif // SERIALIZATION


#endif
