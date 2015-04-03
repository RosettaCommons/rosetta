// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/PoseMetricCalculators/SurfaceCalculator.cc
/// @brief  A class which will keep track of the SASA-based hpatch score of a Pose object
/// @author Ron Jacak

// Unit headers
#include <protocols/toolbox/pose_metric_calculators/HPatchCalculator.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pack/interaction_graph/SurfacePotential.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>


#include <utility/assert.hh>
#include <iostream>
#include <sstream>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


using namespace core;

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

HPatchCalculator::HPatchCalculator( bool remove_nonprotein_res )
	: total_hpatch_score_(0.0), remove_nonprotein_res_(remove_nonprotein_res )
{}

HPatchCalculator::~HPatchCalculator(){}

void
HPatchCalculator::lookup( std::string const & key, basic::MetricValueBase* valptr ) const {

	if ( key == "total_hpatch" ) {
		basic::check_cast( valptr, &total_hpatch_score_, "total_hpatch expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( total_hpatch_score_ );

	} else if ( key == "patch_scores" ) {
		basic::check_cast( valptr, &patch_scores_, "patch_scores expects to return a std::map< Size, Real >" );
		(static_cast<basic::MetricValue<std::map< Size, std::pair< Real, Real > > > *>(valptr) )->set( patch_scores_ );

	} else if ( key == "atoms_in_patches" ) {
		basic::check_cast( valptr, &atoms_in_patches_, "atoms_in_patches expects to return a std::map< Size, utility::vector1< id::AtomID > >" );
		(static_cast<basic::MetricValue< std::map< Size, utility::vector1< id::AtomID > > > *>(valptr) )->set( atoms_in_patches_ );

	} else {
		basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
		utility_exit();
	}

}


std::string
HPatchCalculator::print( std::string const & key ) const {

	if ( key == "total_hpatch" ) {
		std::ostringstream ss;
		ss << utility::to_string( total_hpatch_score_ ) << " (UNWEIGHTED)";
		return ss.str();
		//return utility::to_string( total_hpatch_score_ );

	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	return "";

}


void
HPatchCalculator::recompute( pose::Pose const & this_pose ) {

	// use the SurfacePotential class to recompute the hpatch score
	if( !remove_nonprotein_res_ ){
		pack::interaction_graph::SurfacePotential::get_instance()->compute_pose_hpatch_score( this_pose, total_hpatch_score_, patch_scores_, atoms_in_patches_ );

	} else {
		// iterate over all residues to see if this pose has any non-protein residues
		bool has_nonprot_res(false);
		for ( core::Size ii = 1; ii <= this_pose.total_residue(); ++ii ) {
			if ( ! this_pose.residue_type(ii).is_protein() ) {
				has_nonprot_res = true;
				break;
			}
		}
		if ( has_nonprot_res ) {
			core::pose::PoseOP pureprotpose( new core::pose::Pose( this_pose ) );
			pose_manipulation::remove_non_protein_residues( *pureprotpose );
			pack::interaction_graph::SurfacePotential::get_instance()->compute_pose_hpatch_score( *pureprotpose, total_hpatch_score_, patch_scores_, atoms_in_patches_ );

		} else {
			pack::interaction_graph::SurfacePotential::get_instance()->compute_pose_hpatch_score( this_pose, total_hpatch_score_, patch_scores_, atoms_in_patches_ );
		}
	}
}


} // PoseMetricCalculators
} // toolbox
} // protocols

