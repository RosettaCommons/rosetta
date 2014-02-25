// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_CentroidInfo.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_rna_RNA_CentroidInfo_hh
#define INCLUDED_core_scoring_rna_RNA_CentroidInfo_hh

#include <core/types.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.hh>

// Utility headers

// Numceric Headers
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


// C++

namespace core {
namespace scoring {
namespace rna {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.
//// Rhiju move this to its own namespace!
class RNA_CentroidInfo: public basic::datacache::CacheableData  {

public:

RNA_CentroidInfo(): calculated_( false ) {};

  RNA_CentroidInfo( RNA_CentroidInfo const & src );

  basic::datacache::CacheableDataOP
  clone() const
  {
    return new RNA_CentroidInfo( *this );
  }

	void
	update( pose::Pose const & pose );

  Size
  size() const {
    return base_centroids_.size();
  }

  bool
  calculated() const
  {
    return calculated_;
  }

  bool &
  calculated()
  {
    return calculated_;
  }

  void
  set_calculated( bool const & setting )
  {
    calculated_ = setting;
  }

  utility::vector1< Vector > const &
	base_centroids() const
	{
		return base_centroids_;
	}

  utility::vector1< kinematics::Stub > const &
	base_stubs() const
	{
		return base_stubs_;
	}

 	Vector
	get_base_centroid( conformation::Residue const & rsd ) const;

	kinematics::Stub
	get_base_coordinate_system( conformation::Residue const & rsd, Vector const & centroid ) const;

	kinematics::Stub
	get_base_coordinate_system( conformation::Residue const & rsd ) const;

private:

  void
  initialize_base_centroids_and_stubs( pose::Pose const & pose );

  utility::vector1< Vector > base_centroids_;
  utility::vector1< kinematics::Stub > base_stubs_;
  bool calculated_;

};


}
}
}
#endif
