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
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/StructureDependentPeakCalibrator.hh>

// Package Headers
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>

// Project Headers
#include <core/scoring/constraints/Constraint.hh>

#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/numbers.hh>

//// C++ headers
#include <cmath>
#include <iomanip>

static thread_local basic::Tracer tr( "protocols.noesy_assign.calibration" );

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
	runtime_assert( structures_.size() )
		core::pose::Pose const& pose( **(structures_.begin()) );
	for ( utility::vector1< CrossPeakOP >::const_iterator it = peaks().begin(); it != peaks().end(); ++it, ++ct ) {
		(*it)->create_fa_and_cen_constraint( constraints_[ ct ], dummy, pose, dummy_pose, 1, 0.0 /*padding*/, true /*only fa cst*/ );
	}
}

void StructureDependentPeakCalibrator::collect_upperbound_statistics( core::Size peak, TypeCumulator const& types ) {
	Size violated( 0 );
	Real inv_n_struct( 1.0 / structures_.size() );
	runtime_assert( peak <= peaks().size() );
	runtime_assert( constraints_.size() == peaks().size() );
	Size pose_ct( 1 );
	Real stddev( 0.0);
	Real mean( 0.0 );
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	if ( constraints_[ peak ] &&
			!( params.calibration_ignore_eliminated_peaks_  && peaks()[ peak ]->eliminated() ) ) {
		for ( PoseVector::const_iterator pose_it = structures_.begin(); pose_it != structures_.end(); ++pose_it, ++pose_ct ) {
			Real dist( constraints_[ peak ]->dist( **pose_it ) );
			stddev += dist*dist;
			mean += dist;
			violated += ( dist - peaks()[ peak ]->distance_bound() ) > dcalibrate_;
			//   if ( pose_ct == 1 && tr.Trace.visible() ) tr.Trace << peaks()[ peak ]->peak_id() << " " << peaks()[ peak ]->filename() << " sum_dist " << constraints_[ peak ]->dist( **pose_it ) << std::endl;
		}
		mean *= inv_n_struct;
		stddev = stddev*inv_n_struct - mean*mean;
		// tr.Debug << "peak: " << peaks()[ peak ]->peak_id() << " " << peaks()[ peak ]->filename() << " violated: " << violated << " " << 1.0/inv_n_struct << " " << std::endl;
		if ( stddev < params.calibration_convergence_ || params.calibration_convergence_ < 0.01 ) {
			collect_target_statistics( violated*inv_n_struct, types );
		}
		//   for ( core::Size type = BACKBONE; type < MAX_TYPE; ++type ) {
		//    if ( types.test( type ) ) {
		//     accumulated_count_[ type ] += 1;
		//     accumulated_target_[ type ] += violated * inv_n_struct;
		//    }
		//   }
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

		if ( params.calibration_max_nudging_ > 1.0 && violated > params.calibration_start_nudging_*structures_.size() ) {
			tr.Trace << "Check peak " << (*it)->peak_id() << " for nudging... "<< std::endl;
			Real const CORRECTION_STEP( 0.1 );
			Real const max_correction( peaks()[ ct ]->distance_bound()*( params.calibration_max_nudging_ - 1) );
			Size const old_violated( violated );
			for ( Real correction = 0.1; correction <= max_correction; correction += CORRECTION_STEP ) {
				violated = 0;
				Size pose_ct( 1 );
				for ( PoseVector::const_iterator pose_it = structures_.begin(); pose_it != structures_.end(); ++pose_it, ++pose_ct ) {
					Real delta( distance_deltas[ pose_ct ] - correction );
					violated += delta > params.dcut_;
				}
				if ( violated <= params.calibration_stop_nudging_*structures_.size() ) {
					peaks()[ ct ]->nudge_distance_bound( correction );
					Size pose_ct( 1 );
					for ( PoseVector::const_iterator pose_it = structures_.begin(); pose_it != structures_.end(); ++pose_it, ++pose_ct ) {
						distance_deltas[ pose_ct ] -= correction;
					}
					tr.Debug << "peak: " << (*it)->peak_id() <<" " << (*it)->filename()
						<< " original violations: " << old_violated
						<< std::setprecision(2) << " new distance: " << peaks()[ ct ]->distance_bound()
						<< " nudged by: " << correction
						<< std::setprecision(2) << " of max " << max_correction
						<< " new violations: " << violated << std::endl;
					break;
				}
				violated=old_violated;
			}
		}

		if ( !params.use_local_distviol_ )  {
			tr.Debug << "peak: " << (*it)->peak_id() <<" " << (*it)->filename() << " violations: " << violated << std::endl;
			(*it)->set_eliminated_due_to_dist_violations( violated > ( params.nr_conformers_violatable_*structures_.size() ) );
			std::ostringstream elim_msg;
			std::sort( distance_deltas.begin(), distance_deltas.end() );
			core::Size median_position( (core::Size)(utility::round( 0.5*distance_deltas.size()+0.5 )) );
			elim_msg << violated << " ("<<distance_deltas.size()<<") violated by >" << distance_deltas[median_position] << "A (" << params.dcut_ << "A) ";
			(*it)->set_elimination_comment( elim_msg.str() );
		} else {  //local dist viol

			//first sort to get a 90% distribution length --i.e., ignore 5% on each side and take distance between those
			std::sort( distance_deltas.begin(), distance_deltas.end() );

			//find smallest interval that fits 99% of the deltas
			//with default setting of 99% this is basically the length difference between shortest and longest distance
			Size const num_element_cluster( (core::Size)(utility::round( 1.0*distance_deltas.size() * params.local_distviol_range_  )) );
			Size const low_quartil_pos( (core::Size)(utility::round( 1.0*distance_deltas.size()*0.25 )) );
			Real const low_quartil_dist( distance_deltas[ low_quartil_pos ]+(*it)->distance_bound() );
			tr.Debug << "peak: " << (*it)->peak_id() << " " << (*it)->filename() << " check " << num_element_cluster << " of a total " << distance_deltas.size() << " distances for max-extension " << std::endl;
			Real max_extension( 1000 );
			if ( low_quartil_dist > params.local_distviol_cutoff_ || low_quartil_dist > (*it)->distance_bound() + params.local_distviol_cutoff_buffer_ ) {
				tr.Debug << "peak: " << (*it)->peak_id() << " " << (*it)->filename() << " dist " << (*it)->distance_bound() << " REMOVED due to large Q1 dist of " << low_quartil_dist << std::endl;
				(*it)->set_eliminated_due_to_dist_violations( true );
				std::ostringstream elim_msg;
				elim_msg << "Q1 dist to high: " << low_quartil_dist;
				(*it)->set_elimination_comment( elim_msg.str() );
			} else {
				for ( Size start_cluster = 1; start_cluster+num_element_cluster-1 <= distance_deltas.size(); start_cluster++ ) {
					Real ext = distance_deltas[ start_cluster+num_element_cluster-1 ] - distance_deltas[ start_cluster ];
					if ( max_extension > ext ) max_extension = ext;
				}

				tr.Debug << num_element_cluster << " distances are in an interval of only " << max_extension << " with a Q1 dist of " << low_quartil_dist << std::endl;
				//get extension between high and low.
				//   Size const ind_low_5( 1+utility::round( params.local_distviol_range_*distance_deltas.size() ) );   //lower 5% -
				//   Size const ind_high_5( utility::round( 1.0*distance_deltas.size()*(1-params.local_distviol_range_) ) ); //upper 5%
				//   Real max_extension( distance_deltas[ ind_high_5 ] - distance_deltas[ ind_low_5 ] );

				//   tr.Debug << "ind_low_5 " << ind_low_5 << " ind_high_5 " << ind_high_5 << " min_delta: "
				//        << distance_deltas[ 1 ] << " max_delta "
				//<< distance_deltas.back() << std::endl;

				Size viol_count( 0 );
				tr.Trace << " dist: " << (*it)->distance_bound() << "| " ;
				core::Real violation_cutoff( max_extension * params.local_distviol_global_factor_ + params.local_distviol_global_buffer_ );
				for ( RealVector::const_iterator delta_it = distance_deltas.begin(); delta_it != distance_deltas.end(); ++delta_it ) {
					tr.Trace << " " << *delta_it;
					viol_count += ( *delta_it > violation_cutoff ) ? 1 : 0;
				}
				tr.Trace << std::endl;

				tr.Debug << "peak: " << (*it)->peak_id() <<" " << (*it)->filename() << " dist: " << (*it)->distance_bound()
					<< " max_extension " << max_extension << " viol_count " << viol_count
					<<  ( viol_count > (params.nr_conformers_violatable_*distance_deltas.size() ) ? " REMOVED " : "" ) << std::endl;

				(*it)->set_eliminated_due_to_dist_violations( viol_count > ( params.nr_conformers_violatable_*distance_deltas.size() ) );
				std::ostringstream elim_msg;
				elim_msg << viol_count << " ("<<distance_deltas.size()<<") violated by >" << distance_deltas[1] << "A (" << violation_cutoff << "A) ";
				(*it)->set_elimination_comment( elim_msg.str() );
			}
			//what is an elimination candidate ?
			(*it)->set_elimination_candidate( violated > ( params.nr_conformers_violatable_*structures_.size() ) );
		} // use_local_distviol
	} // for peaks
}

}
}
