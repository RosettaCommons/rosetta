// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/PoseMetricCalculators/SurfaceCalculator.cc
/// @brief  A class which will keep track of surface energy as a Pose metric
/// @author Ron Jacak

// Unit headers
#include <protocols/toolbox/pose_metric_calculators/SurfaceCalculator.hh>
#include <core/chemical/ResidueType.hh>
#include <basic/MetricValue.hh>
#include <core/pose/Pose.hh>
#include <core/pack/interaction_graph/SurfacePotential.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

#include <ObjexxFCL/format.hh>

#include <utility/assert.hh>
#include <iostream>
#include <sstream>

#include <utility/vector1.hh>


using namespace core;

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

SurfaceCalculator::SurfaceCalculator( bool remove_nonprotein_res )
: total_surface_energy_(0.0), remove_nonprotein_res_(remove_nonprotein_res )
{}

SurfaceCalculator::~SurfaceCalculator(){}

void
SurfaceCalculator::lookup( std::string const & key, basic::MetricValueBase* valptr ) const {

	if ( key == "total_surface" ) {
		basic::check_cast( valptr, &total_surface_energy_, "total_surface expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( total_surface_energy_ );

	} else if ( key == "residue_surface" ) {
		basic::check_cast( valptr, &residue_surface_energy_, "residue_surface expects to return a utility::vector1< Real >" );
		(static_cast<basic::MetricValue<utility::vector1<Real > > *>(valptr) )->set( residue_surface_energy_ );

	} else {
		basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
		utility_exit();
	}

}


std::string
SurfaceCalculator::print( std::string const & key ) const {

	if ( key == "total_surface" ) {
		std::ostringstream ss;
		ss << utility::to_string( total_surface_energy_ ) << " (UNWEIGHTED)";
		return ss.str();
		//return utility::to_string( total_surface_energy_ );

	} else if ( key == "residue_surface" ) {
		std::ostringstream ss;
		ss << "[";
		for ( Size ii=1; ii < residue_surface_energy_.size(); ++ii ) {
			ss << ii << ":";
			if ( residue_surface_energy_[ ii ] != 0.00 ) {
				ss << ObjexxFCL::format::F(8,4, residue_surface_energy_[ ii ]);
			} else {
				ss << "---";
			}
			ss << ",  ";
		}
		ss << "] (UNWEIGHTED)";
		return ss.str();
		//return utility::to_string( residue_surface_energy_ );
	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	return "";

}


void
SurfaceCalculator::recompute( pose::Pose const & this_pose ) {

	// use the SurfacePotential class to recompute the surface energy
	if ( !remove_nonprotein_res_ ) {
		pack::interaction_graph::SurfacePotential::get_instance()->compute_pose_surface_energy( this_pose, total_surface_energy_, residue_surface_energy_ );

	} else {
		bool has_nonprot_res(false);
		for ( core::Size i = 1; i <= this_pose.total_residue(); ++i ) {
			if ( ! this_pose.residue_type(i).is_protein() ) {
				has_nonprot_res = true;
				break;
			}
		}
		if ( has_nonprot_res ) {
			core::pose::PoseOP pureprotpose( new core::pose::Pose( this_pose ) );
			pose_manipulation::remove_non_protein_residues( *pureprotpose );
			pureprotpose->update_residue_neighbors();
			pack::interaction_graph::SurfacePotential::get_instance()->compute_pose_surface_energy( *pureprotpose, total_surface_energy_, residue_surface_energy_ );
		} else {
			pack::interaction_graph::SurfacePotential::get_instance()->compute_pose_surface_energy( this_pose, total_surface_energy_, residue_surface_energy_ );
		}
	}
}


} // PoseMetricCalculators
} // toolbox
} // protocols

