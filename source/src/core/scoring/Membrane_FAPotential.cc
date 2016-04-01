// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/Membrane_FAPotential.cc
///
/// @brief  Membrane FA Potential - Class for Fullatom Membrane Scoring Methods
/// @details Compute High resolution energy terms and high resolution embedding corrections
///    for penalties. Also contains pass-through methods for accessing and updating
///    mp framework supported data in a membrane conformation.
///    Last Modified: 3/11/14
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @author  Patrick Barth (original)

// Unit headers
#include <core/scoring/Membrane_FAPotential.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/EnvPairPotential.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh>

// Utility headers
#include <core/types.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.Membrane_FAEmbed" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

Membrane_FAPotential::Membrane_FAPotential() :
	membrane_potential_( ScoringManager::get_instance()->get_MembranePotential() )
{
}


/// @brief Copy Constructor for Membrane Fullatom Embedding
Membrane_FAEmbed::Membrane_FAEmbed( Membrane_FAEmbed const & src ) :
	CacheableData()
{
	fa_center_ = src.fa_center_;
	fa_proj_ = src.fa_proj_;
	fa_proj_coord_ = src.fa_proj_coord_;
	fa_proj_deriv_ = src.fa_proj_deriv_;
	fa_depth_ = src.fa_depth_;
	calculated_ = src.calculated_;
	Fa_Membed_update_ = src.Fa_Membed_update_;
}

/// @brief Initialize Mmebrane Fullatom Embedding from Options
void Membrane_FAEmbed::initialize(pose::Pose const & pose)
{
	thickness_ = basic::options::option[ basic::options::OptionKeys::membrane::thickness ]();
	steepness_ = 10.0;
	Fa_Membed_update_ = basic::options::option[ basic::options::OptionKeys::membrane::Fa_Membed_update ]();
	allocate_appropriate_memory( pose );
}

/// @brief Setup Data Members for Appropriate Sizing
void Membrane_FAEmbed::allocate_appropriate_memory(pose::Pose const & pose) const
{
	fa_center_ = 0.0;
	Size const nres( pose.total_residue() );
	fa_proj_.resize( nres );
	fa_proj_coord_.resize( nres );
	fa_proj_deriv_.resize( nres );
	fa_depth_.resize( nres );

	// Tryptophan is the largest...
	static Size const MAX_AMINOACID_SIZE = 15; ///Todo Fix this Hack. Membrane Design failed because these vectors are presized assuming residue atom count remains constant
	for ( Size i = 1; i <= nres; ++i ) {
		Size const max_size = std::max( MAX_AMINOACID_SIZE, pose.residue( i ).nheavyatoms() ); // in case of ligands or NCAAs
		fa_proj_[i].resize( max_size );
		fa_proj_coord_[i].resize( max_size );
		fa_proj_deriv_[i].resize( max_size );
		fa_depth_[i].resize( max_size);
		for ( Size j = 1; j <= max_size; ++j ) {
			fa_proj_[i][j] = 0.0;
			fa_proj_coord_[i][j].assign(0.0,0.0,0.0);
			fa_proj_deriv_[i][j] = 0.0;
			fa_depth_[i][j] = 0.0;
		}
	}
}

/// @brief Finalize Setup of MP Potential
void Membrane_FAPotential::finalize( pose::Pose & pose ) const
{
	CenListInfo & cenlist( nonconst_cenlist_from_pose( pose ));
	cenlist.calculated() = true; //pba

	core::scoring::MembraneEmbed & membrane_embed( core::scoring::nonconst_MembraneEmbed_from_pose( pose ));
	membrane_embed.calculated() = false; //pba 112209
	Membrane_FAEmbed & membrane_faembed( nonconst_Membrane_FAEmbed_from_pose( pose ));
	membrane_faembed.calculated() = false;
}

/// @brief Compute FullAtom TM projection
void Membrane_FAPotential::compute_fa_projection(pose::Pose & pose) const
{
	core::scoring::MembraneEmbed & membrane_embed(core::scoring::nonconst_MembraneEmbed_from_pose( pose ));
	Membrane_FAEmbed & membrane_faembed(nonconst_Membrane_FAEmbed_from_pose( pose ));

	membrane_faembed.initialize(pose);

	if ( membrane_faembed.Fa_Membed_update() ) membrane_embed.calculated() = false;
	membrane_potential_.compute_membrane_embedding(pose);

	Vector const normal(MembraneEmbed_from_pose( pose ).normal());
	Vector const center(MembraneEmbed_from_pose( pose ).center());
	Real const penalty(MembraneEmbed_from_pose( pose ).penalty());
	Real const thickness(Membrane_FAEmbed_from_pose( pose ).thickness());
	Real const steepness(Membrane_FAEmbed_from_pose( pose ).steepness());

	fa_projection(pose,normal,center,thickness,steepness,penalty);
}

/// @brief Helper function called by compute_fa_projection
void Membrane_FAPotential::fa_projection(
	pose::Pose & pose,
	Vector const & normal,
	Vector const & center,
	Real const & thickness,
	Real const & steepness,
	Real const & penalty
) const {
	// mjo commenting out 'topology' because it is unused and causes a warning
	Membrane_FAEmbed & membrane_faembed(nonconst_Membrane_FAEmbed_from_pose( pose ));

	Size nres=pose.total_residue();
	Real internal_product(0), z(0), zn(0), znm1(0);

	membrane_faembed.fa_center() = std::abs(dot(center, normal));
	membrane_faembed.fa_penalty() = penalty;

	//pbadebug
	if ( !membrane_faembed.calculated() ) {
		// USE TRACER OUTPUT ! ! !
		//     std::cout << "CENTER " << center.x() <<  " " << center.y() << " " << center.z() << "\n";
		//     std::cout << "NORMAL " << normal.x() <<  " " << normal.y() << " " << normal.z() << "\n";
		//     std::cout << "thickness steepness " << thickness <<  " " << steepness << "\n";
		//     std::cout << "fa_center " << membrane_faembed.fa_center() << "\n";
		// std::cout << "fa_penalty " << membrane_faembed.fa_penalty() << "\n";
	}

	for ( Size i = 1; i <= nres; ++i ) {
		for ( Size j = 1, j_end = pose.residue( i ).nheavyatoms(); j <= j_end; ++j ) {
			Vector const xyz( pose.residue( i ).xyz( j ) );
			membrane_faembed.fa_depth( i, j ) = dot(xyz-center, normal);

			internal_product = std::abs(membrane_faembed.fa_depth( i, j ));
			z = internal_product;
			z /= thickness;
			zn = std::pow( z, steepness );
			membrane_faembed.fa_proj( i, j ) = zn/(1 + zn);

			//deriv
			znm1 = std::pow( z, (steepness-1) );
			membrane_faembed.fa_proj_deriv( i, j ) = steepness * znm1 * std::pow((1 + zn),-2);
			membrane_faembed.fa_proj_deriv( i, j ) /= thickness;

			Vector proj_i = center + membrane_faembed.fa_depth(i,j)*normal;
			Vector i_ip = proj_i - xyz;
			Vector proj_center = center - i_ip;
			membrane_faembed.fa_proj_coord(i,j) = proj_center;
		}
	}
	membrane_faembed.calculated() = true;
}


/// @brief Return a Const Reference to the Embedding Object from the Pose Cache
Membrane_FAEmbed const & Membrane_FAEmbed_from_pose( pose::Pose const & pose )
{
	return *( utility::pointer::static_pointer_cast< core::scoring::Membrane_FAEmbed const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED ) ));
}

/// @brief Return a Non Const Reference to the Embedding Object from the Pose Cache
Membrane_FAEmbed & nonconst_Membrane_FAEmbed_from_pose( pose::Pose & pose )
{
	if ( pose.data().has( core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED ) ) {
		return *( utility::pointer::static_pointer_cast< core::scoring::Membrane_FAEmbed > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED ) ));
	}

	Membrane_FAEmbedOP membrane_faembed( new Membrane_FAEmbed );
	pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED, membrane_faembed );
	return *membrane_faembed;
}

} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::Membrane_FAEmbed::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( fa_proj_ ) ); // utility::vector1<utility::vector1<Real> >
	arc( CEREAL_NVP( fa_depth_ ) ); // utility::vector1<utility::vector1<Real> >
	arc( CEREAL_NVP( fa_proj_coord_ ) ); // utility::vector1<utility::vector1<Vector> >
	arc( CEREAL_NVP( fa_proj_deriv_ ) ); // utility::vector1<utility::vector1<Real> >
	arc( CEREAL_NVP( calculated_ ) ); // _Bool
	arc( CEREAL_NVP( fa_center_ ) ); // Real
	arc( CEREAL_NVP( fa_penalty_ ) ); // Real
	arc( CEREAL_NVP( thickness_ ) ); // Real
	arc( CEREAL_NVP( steepness_ ) ); // Real
	arc( CEREAL_NVP( Fa_Membed_update_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::Membrane_FAEmbed::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( fa_proj_ ); // utility::vector1<utility::vector1<Real> >
	arc( fa_depth_ ); // utility::vector1<utility::vector1<Real> >
	arc( fa_proj_coord_ ); // utility::vector1<utility::vector1<Vector> >
	arc( fa_proj_deriv_ ); // utility::vector1<utility::vector1<Real> >
	arc( calculated_ ); // _Bool
	arc( fa_center_ ); // Real
	arc( fa_penalty_ ); // Real
	arc( thickness_ ); // Real
	arc( steepness_ ); // Real
	arc( Fa_Membed_update_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::Membrane_FAEmbed );
CEREAL_REGISTER_TYPE( core::scoring::Membrane_FAEmbed )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_Membrane_FAPotential )
#endif // SERIALIZATION
