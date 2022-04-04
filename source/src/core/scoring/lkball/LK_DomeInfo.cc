// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/LK_BallInfo.cc
/// @brief  Orientation dependent variant of the LK Solvation using
/// @author Phil Bradley

// Unit headers
#include <core/scoring/lkball/LK_DomeInfo.hh>
#include <core/scoring/lkball/LK_DomeEnergy.hh>

// // Package headers
//#include <core/pack/rotamer_set/WaterPackingInfo.hh>

// // Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/util.hh>


#include <basic/Tracer.hh>

#include <sstream>

// #include <utility/vector1.functions.hh> // HACK

//#ifdef WIN32
//#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
//#endif

#ifdef    SERIALIZATION

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/fixedsizearray1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace lkball {

//bool const sidechain_only_hack( false ); // make this configurable

/// @details Auto-generated virtual destructor
LKD_ResidueInfo::~LKD_ResidueInfo() = default;


static basic::Tracer TR("core.scoring.methods.LK_DomeInfo" );


/// LAZY
using namespace chemical;
using namespace pose;
using namespace conformation;


LKD_ResidueInfo::LKD_ResidueInfo(
	conformation::Residue const & rsd,
	LK_DomeEnergy const * lk_dome
)
: LKB_ResidueInfo( rsd ),
	water_occlusions_( rsd.nheavyatoms() ),
	water_sol_values_( rsd.nheavyatoms() )
{

	if ( ! has_waters() ) return;

	// debug_assert( waters().size() == rsd.nheavyatoms() );

	for ( Size iatom = 1; iatom <= rsd.nheavyatoms(); iatom++ ) {
		water_sol_values_[iatom] = lk_dome->get_sol_value( rsd.atom_type( iatom ), n_attached_waters()[iatom] );
	}
}

LKD_ResidueInfo::LKD_ResidueInfo( LKD_ResidueInfo const & src )
: LKB_ResidueInfo( src ),
	water_occlusions_( src.water_occlusions_ ),
	water_sol_values_( src.water_sol_values_ )
{
}

LKD_ResidueInfo &
LKD_ResidueInfo::operator=(LKD_ResidueInfo const & src) {
	LKB_ResidueInfo::operator=( src );

	water_occlusions_ = src.water_occlusions_;
	water_sol_values_ = src.water_sol_values_;

	return *this;
}

LKD_ResidueInfo::LKD_ResidueInfo() = default;



basic::datacache::CacheableDataOP
LKD_ResidueInfo::clone() const {
	return utility::pointer::make_shared<LKD_ResidueInfo>( *this );
}


utility::vector1< WaterOcclusions > &
LKD_ResidueInfo::water_occlusions() {
	return water_occlusions_;
}

utility::vector1< WaterOcclusions > const &
LKD_ResidueInfo::water_occlusions() const {
	return water_occlusions_;
}

utility::vector1< Real > const &
LKD_ResidueInfo::water_sol_values() const {
	return water_sol_values_;
}



void
LKD_ResidueInfo::clear_occlusions() {
	for ( Size i = 1; i <= water_occlusions_.size(); i++ ) {
		for ( Size j = 1; j <= n_attached_waters()[i]; j++ ) {
			water_occlusions_[i][j] = 0;
		}
	}
}

}
}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::lkball::LKD_ResidueInfo::save( Archive & arc ) const {
	arc( cereal::base_class< LKB_ResidueInfo >( this ) );
	arc( CEREAL_NVP( water_occlusions_ ) );
	arc( CEREAL_NVP( water_sol_values_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::lkball::LKD_ResidueInfo::load( Archive & arc ) {
	arc( cereal::base_class< LKB_ResidueInfo >( this ) );
	arc( water_occlusions_ );
	arc( water_sol_values_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::lkball::LKD_ResidueInfo );
CEREAL_REGISTER_TYPE( core::scoring::lkball::LKD_ResidueInfo )



CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_lkball_LK_DomeInfo )
#endif // SERIALIZATION
