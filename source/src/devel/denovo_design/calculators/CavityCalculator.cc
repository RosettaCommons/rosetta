// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/denovo_design/calculators/CavityCalculator.cc
/// @brief calculate centrality for all (or a subset of) residues
/// @details
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit Headers
#include <devel/denovo_design/calculators/CavityCalculator.hh>

// Project Headers
#include <basic/MetricValue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/simple_moves/AddCavitiesMover.hh>


// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <boost/lexical_cast.hpp>

//// C++ headers
#include <ObjexxFCL/format.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>


static thread_local basic::Tracer TR( "devel.denovo_design.calculators.CavityCalculator" );

namespace devel {
namespace denovo_design {
namespace calculators {

/// @brief default constructor
CavityCalculator::CavityCalculator() :
  core::pose::metrics::StructureDependentCalculator(),
  total_volume_( -1.0 )
{
}

/// @brief destructor
CavityCalculator::~CavityCalculator(){}

/// @brief
void
CavityCalculator::lookup( std::string const & key,
			  basic::MetricValueBase * valptr ) const
{
  if ( key == "volume" ) {
    basic::check_cast( valptr, &total_volume_, "Total Cavity Volume" );
    (static_cast< basic::MetricValue< core::Real > *>(valptr))->set( total_volume_ );
  } else if ( key == "cavities" ) {
		basic::check_cast( valptr, &clusters_, "Cavities" );
		(static_cast< basic::MetricValue< utility::vector1< core::scoring::packstat::CavityBallCluster > > *>(valptr))->set( clusters_ );
	} else {
    TR << "CavityCalculator cannot compute the requested metric " << key << std::endl;
    utility_exit();
  }
} //lookup


// @brief
std::string
CavityCalculator::print( std::string const & key ) const
{
  std::string result;
  if ( key == "volume" ) {
    result += boost::lexical_cast<std::string>( total_volume_ ) + " ";
	} else {
    basic::Error() << "CavityCalculator cannot compute metric " << key << std::endl;
  }
  return result;
} // print

/// @brief recompute cavity volumes
void
CavityCalculator::recompute( core::pose::Pose const & pose )
{
  // convert to centroid and compare to saved pose; if it is the same, don't do anything.
  total_volume_ = 0;

	// setup sasa options -- same as in AddCavitiesMover
	core::scoring::packstat::SasaOptions opts;
	opts.prune_max_iters = 5;
	opts.prune_max_delta = 0;
	opts.num_cav_ball_layers = 2;
	opts.frac_cav_ball_required_exposed = 0.0;
	opts.area_cav_ball_required_exposed = 0.0;
	opts.surrounding_sasa_smoothing_window = 1;

  protocols::simple_moves::AddCavitiesMover add_cavities;
  // parameters come from AddCavitiesMover
  // 10.0 is hard coded, 150 is default, 3.0 is hard coded
  core::scoring::packstat::CavBalls cbs( add_cavities.get_cavities( pose, 10.0, 150, 3.0 ) );
  // hack required so CavBalls object is not modified
  core::scoring::packstat::CavBalls tmp_cbs( cbs );
  compute_cav_ball_volumes( tmp_cbs, opts );

	clusters_.clear();
	clusters_ = compute_cav_ball_clusters( tmp_cbs, opts );
  for( Size i = 1; i <= clusters_.size(); i++ ) {
    for( Size j = 1; j <= clusters_[i].cavballs.size(); j++ ) {
      if( clusters_[i].cavballs[j].radius() > 0.6 )
				TR    << clusters_[i].cavballs[j].hetero_atom_line( pose.total_residue()+i, i, 0.0 ) << std::endl;
    }
		TR << "Cluster " << i << ": volume=" << clusters_[i].volume << ", surface_area=" << clusters_[i].surface_area << ", surface_accessibility=" << clusters_[i].surface_accessibility << ", center=" << clusters_[i].center.x() << "," << clusters_[i].center.y() << "," << clusters_[i].center.z() << std::endl;
		total_volume_ += clusters_[i].volume;
  }
}

} // calculators
} // denovo_design
} // devel

