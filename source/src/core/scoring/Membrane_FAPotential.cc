// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   core/scoring/methods/EnvPairPotential.cc
/// @brief  Membrane FA Potential
/// @author Patrick Barth
///


// Unit headers
#include <core/scoring/Membrane_FAPotential.hh>
#include <core/scoring/Membrane_FAPotential.fwd.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/EnvPairPotential.hh>

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


#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>



// just for debugging
//#include <ObjexxFCL/formatted.o.hh>

// C++


namespace core {
namespace scoring {
//static numeric::random::RandomGenerator RG(280628);  // <- Magic number, do not change it!


//pba

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

void
Membrane_FAEmbed::initialize(pose::Pose const & pose)
{
  thickness_=basic::options::option[ basic::options::OptionKeys::membrane::thickness ]();
  steepness_=10.0;
  Fa_Membed_update_=basic::options::option[ basic::options::OptionKeys::membrane::Fa_Membed_update ]();
  allocate_appropriate_memory( pose );
}

void
Membrane_FAEmbed::allocate_appropriate_memory(pose::Pose const & pose) const
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

void
Membrane_FAPotential::finalize( pose::Pose & pose ) const
{
  CenListInfo & cenlist( nonconst_cenlist_from_pose( pose ));
  cenlist.calculated() = true; //pba
  core::scoring::MembraneEmbed & membrane_embed( core::scoring::nonconst_MembraneEmbed_from_pose( pose ));
  membrane_embed.calculated() = false; //pba 112209
  Membrane_FAEmbed & membrane_faembed( nonconst_Membrane_FAEmbed_from_pose( pose ));
  membrane_faembed.calculated() = false; //pba
}

void
Membrane_FAPotential::compute_fa_projection(pose::Pose & pose) const
{

  core::scoring::MembraneEmbed & membrane_embed(core::scoring::nonconst_MembraneEmbed_from_pose( pose ));
  Membrane_FAEmbed & membrane_faembed(nonconst_Membrane_FAEmbed_from_pose( pose ));

  membrane_faembed.initialize(pose);

  if (membrane_faembed.Fa_Membed_update()) membrane_embed.calculated() = false;
  membrane_potential_.compute_membrane_embedding(pose);

  Vector const normal(MembraneEmbed_from_pose( pose ).normal());
  Vector const center(MembraneEmbed_from_pose( pose ).center());
  Real const penalty(MembraneEmbed_from_pose( pose ).penalty());
  Real const thickness(Membrane_FAEmbed_from_pose( pose ).thickness());
  Real const steepness(Membrane_FAEmbed_from_pose( pose ).steepness());

  fa_projection(pose,normal,center,thickness,steepness,penalty);
}

void
Membrane_FAPotential::fa_projection(
  pose::Pose & pose,
  Vector const & normal,
  Vector const & center,
  Real const & thickness,
  Real const & steepness,
  Real const & penalty
) const
{
	// mjo commenting out 'topology' because it is unused and causes a warning
	//core::scoring::MembraneTopology const & topology( core::scoring::MembraneTopology_from_pose(pose) );
  Membrane_FAEmbed & membrane_faembed(nonconst_Membrane_FAEmbed_from_pose( pose ));

  Size nres=pose.total_residue();
  Real internal_product(0), z(0), zn(0), znm1(0);

  membrane_faembed.fa_center() = std::abs(dot(center, normal));
  membrane_faembed.fa_penalty() = penalty;


  //pbadebug
  if(!membrane_faembed.calculated()) {
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


Membrane_FAEmbed const &
Membrane_FAEmbed_from_pose( pose::Pose const & pose )
{
  // ////using core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED;
  return *( static_cast< Membrane_FAEmbed const * >( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED )() ));
}

/// @details Either returns a non-const reference to the cenlist object already stored
/// in the pose, or creates a new cenist object, places it in the pose, and returns
/// a non-const reference to it.
Membrane_FAEmbed &
nonconst_Membrane_FAEmbed_from_pose( pose::Pose & pose )
{
  // ////using core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED;

  if ( pose.data().has( core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED ) ) {
    return *( static_cast< Membrane_FAEmbed * >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED )() ));
  }
  // else
  Membrane_FAEmbedOP membrane_faembed = new Membrane_FAEmbed;
  pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED, membrane_faembed );
  return *membrane_faembed;
}

}
}
