// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Define a suitable replacement for lrint() on Windows
#if (defined WIN32)
	#include <boost/math/special_functions/round.hpp>
	int lrint(double x) {
	return boost::math::iround(x); }
#endif

// Unit Headers
#include <protocols/noesy_assign/StructureDependentPeakCalibrator.hh>

// Package Headers
// AUTO-REMOVED #include <protocols/noesy_assign/CrossPeakList.hh>
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>

// Project Headers
#include <core/scoring/constraints/Constraint.hh>

#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>

//// C++ headers
#include <cmath>

#include <utility/vector1.hh>


static basic::Tracer tr("protocols.noesy_assign.crosspeaks");

using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

namespace protocols {
namespace noesy_assign {




void StructureDependentPeakCalibrator::init_calibrator() {
	generate_constraints();
}

void StructureDependentPeakCalibrator::generate_constraints() {
	core::Size npeaks( peaks().size() );
	constraints_.resize( npeaks, NULL );
	core::scoring::constraints::ConstraintOP dummy( NULL );
	core::pose::Pose dummy_pose;
	core::Size ct( 1 );
	core::pose::Pose const& pose( **(structures_.begin()) );
	for ( utility::vector1< CrossPeakOP >::const_iterator it = peaks().begin(); it != peaks().end(); ++it, ++ct ) {
		(*it)->create_fa_and_cen_constraint( constraints_[ ct ], dummy, pose, dummy_pose, 1, true /*only fa cst*/ );
	}
}

void StructureDependentPeakCalibrator::collect_upperbound_statistics( core::Size peak, TypeCumulator const& types ) {
	Size violated( 0 );
	Real inv_n_struct( 1.0 / structures_.size() );
	runtime_assert( peak <= peaks().size() );
	runtime_assert( constraints_.size() == peaks().size() );
	Size pose_ct( 1 );
	if ( constraints_[ peak ] ) {
		for ( PoseVector::const_iterator pose_it = structures_.begin(); pose_it != structures_.end(); ++pose_it, ++pose_ct ) {
			violated += ( constraints_[ peak ]->dist( **pose_it ) - peaks()[ peak ]->distance_bound() ) > dcalibrate_;
			if ( pose_ct == 1 ) tr.Trace << peaks()[ peak ]->peak_id() << " " << peaks()[ peak ]->filename() << " sum_dist " << constraints_[ peak ]->dist( **pose_it ) << std::endl;
		}
		//	tr.Debug << "peak: " << peaks()[ peak ]->peak_id() << " " << peaks()[ peak ]->filename() << " violated: " << violated << " " << 1.0/inv_n_struct << " " << std::endl;
		collect_target_statistics( violated*inv_n_struct, types );
// 		for ( core::Size type = BACKBONE; type < MAX_TYPE; ++type ) {
// 			if ( types.test( type ) ) {
// 				accumulated_count_[ type ] += 1;
// 				accumulated_target_[ type ] += violated * inv_n_struct;
// 			}
// 		}
	}
}

void StructureDependentPeakCalibrator::eliminate_violated_constraints() {
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );

	Size ct( 1 );

	typedef utility::vector1< Real > RealVector;
	RealVector distance_deltas( structures_.size(), 0.0 );


	for ( utility::vector1< CrossPeakOP >::const_iterator it = peaks().begin(); it != peaks().end(); ++it, ++ct ) {
		if ( !constraints_[ ct ] ) continue;
		Size violated( 0 );
		Size pose_ct( 1 );
		for ( PoseVector::const_iterator pose_it = structures_.begin(); pose_it != structures_.end(); ++pose_it, ++pose_ct ) {
			Real delta( constraints_[ ct ]->dist( **pose_it ) - peaks()[ ct ]->distance_bound() );
			distance_deltas[ pose_ct ] = delta;
			violated += delta > params.dcut_;
		}
		if ( !params.use_local_distviol_ )  {
			tr.Debug << "peak: " << (*it)->peak_id() <<" " << (*it)->filename() << " violations: " << violated << std::endl;
			(*it)->set_eliminated_due_to_dist_violations( violated > ( params.nr_conformers_violatable_*structures_.size() ) );
		} else {  //local dist viol

			//first sort to get a 90% distribution length --i.e., ignore 5% on each side and take distance between those
			std::sort( distance_deltas.begin(), distance_deltas.end() );

			//get extension between high and low.
			Size const ind_low_5( 1+lrint( params.local_distviol_range_*distance_deltas.size() ) );   //lower 5% -
			Size const ind_high_5( lrint( 1.0*distance_deltas.size()*(1-params.local_distviol_range_) ) ); //upper 5%
			Real max_extension( distance_deltas[ ind_high_5 ] - distance_deltas[ ind_low_5 ] );

			tr.Debug << "ind_low_5 " << ind_low_5 << " ind_high_5 " << ind_high_5 << " min_delta: "
							 << distance_deltas[ 1 ] << " max_delta "
							 << distance_deltas.back() << std::endl;

			Size viol_count( 0 );
			for ( RealVector::const_iterator delta_it = distance_deltas.begin(); delta_it != distance_deltas.end(); ++delta_it ) {
				viol_count += ( *delta_it > ( max_extension + params.local_distviol_global_buffer_ ) ) ? 1 : 0;
			}

			tr.Debug << "peak: " << (*it)->peak_id() <<" " << (*it)->filename()
							 << " max_extension " << max_extension << " viol_count " << viol_count << std::endl;

			(*it)->set_eliminated_due_to_dist_violations( viol_count > (params.nr_conformers_violatable_*distance_deltas.size() ) );
		} // use_local_distviol
	} // for peaks
}

}
}
