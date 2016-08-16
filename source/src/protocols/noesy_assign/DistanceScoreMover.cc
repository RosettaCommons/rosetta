// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/DistanceScoreMover.hh>

// Package Headers
#include <protocols/noesy_assign/CrossPeakList.hh>
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>
#include <protocols/noesy_assign/PeakAssignment.hh>
#include <core/id/Exceptions.hh>
// Project Headers
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>
#include <core/scoring/constraints/AmbiguousNMRConstraint.hh>
// Utility headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>

#include <protocols/noesy_assign/CrossPeak.hh>
#include <utility/vector1.hh>
#include <cmath>


//// C++ headers


static THREAD_LOCAL basic::Tracer tr( "protocols.noesy_assign.DistanceScoreMover" );

using core::Real;
using namespace core;
using namespace basic;

namespace protocols {
namespace noesy_assign {

DistanceScoreMover::DistanceScoreMover( CrossPeakList& cpl, pose::Pose const& pose, core::Real dcut ) :
	cross_peaks_( cpl ),
	count_decoys_( 0 ),
	nr_assignments_( cpl.count_assignments() ),
	peak_violation_counts_( cross_peaks_.size(), 0 ),
	final_dist_power_( 6.0 ), //default
	dcut_( dcut )
{
	basic::ProfileThis doit( basic::NOESY_ASSIGN_DIST_INIT );
	constraints_.reserve( nr_assignments_ );
	// peak_constraints_.reserve( cross_peaks_.size() );

#ifndef WIN32

	for ( CrossPeakList::iterator it = cross_peaks_.begin(); it != cross_peaks_.end(); ++it ) {
		for ( CrossPeak::iterator ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
			try {
				constraints_.push_back( (*ait)->create_constraint( pose ) );
			} catch ( core::id::EXCN_AtomNotFound& excn ) {
				tr.Error << "while setting up constraints in DistanceScoreMover: " << excn << std::endl;
				constraints_.push_back( NULL );
			}
		}
	}
#endif
}


void DistanceScoreMover::prepare_scoring( bool use_for_calibration /*default false */ ) {
#ifndef WIN32
	basic::ProfileThis doit( basic::NOESY_ASSIGN_DIST_PREP_SCORE );
	//assignment_distances_ = VectorReal( nr_assignments_, 0.0 );
	used_for_calibration_ = use_for_calibration;

	Size ct( 1 );
	for ( CrossPeakList::iterator it = cross_peaks_.begin(); it != cross_peaks_.end(); ++it, ++ct ) {
		if ( !used_for_calibration_ ) {
			for ( CrossPeak::iterator ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
				(*ait)->set_decoy_compatibility( 0.0 );
			}
		}
		peak_violation_counts_[ ct ] = 0;
		total_violation_count_ = 0;
	}
	count_decoys_ = 0;
	total_assigned_distances_ = 0;
#endif
}


// void DistanceScoreMover::find_violators_with_individual_dist_cutoff( PoseVector poses ) {
//  basic::ProfileThis doit( basic::NOESY_ASSIGN_DIST_APPLY );
//  PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );

//   SingleConstraints::const_iterator constraint_it( constraints_.begin() );
//  active_peaks_ = 0;

//  Size ct_peaks( 1 );
//  typedef utility::vector1< Real > RealVector;
//  RealVector distance_deltas( poses.size(), 0.0 );

//   for ( CrossPeakList::iterator it = cross_peaks_.begin(); it != cross_peaks_.end(); ++it, ++ct_peaks ) {
//   //for each CrossPeak do the following:
//   // ct_peaks .. number of peak, ct ... number of assignment / peak
//   // a) iterate over all initial assignments and get effective distance (i.e., QD1 is evaluated with d^-6 averaging )
//   // b1) average over all assignments (d^-6) (--> sum_dist )
//   // b2) average over all selected assignments ( calibration phase --> sum_dist_filt )
//   // c) cache distance in dist_buf to compute Dk in the end..
//    if ( (*it)->n_assigned() == 0 ) continue;

//   Size pose_ct( 1 );
//   SingleConstraints::const_iterator pre_pose_cst_it = constraint_it;
//   for ( PoseVector::const_iterator pose_it = poses.begin(); pose_it != poses.end(); ++pose_it, ++pose_ct ) {
//    Real sum_dist( 0.0 ); //accumulate only distances where Vk > Vmin  --> used for calibration ( eq. (12)  ).
//    constraint_it = pre_pose_cst_it;
//    for ( CrossPeak::iterator ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
//     runtime_assert( constraint_it != constraints_.end() );
//     Real dist,invd6;
//     if ( !( *constraint_it )) {
//      dist=0;
//      invd6=0;
//     } else {
//      {
//       dist = (*constraint_it)->dist( **pose_it );
//      }
//      Real invd = 1.0/dist;
//      Real invd3 = invd*invd*invd;
//      invd6 = invd3*invd3;
//     }
//     sum_dist += ( (*ait)->normalized_peak_volume() > params.min_volume_ )*invd6;
//     ++constraint_it; //used this constraint.... we have created these constraints in constructor in same sequence as used here
//    }
//    total_assigned_distances_ += sum_dist;
//    if ( sum_dist > 0 )  {
//     sum_dist=pow( sum_dist, -1.0/6.0 );
//     ++active_peaks_;
//    }
//    distance_deltas[ pose_ct ] = sum_dist - (*it)->distance_bound();
//   } // for all poses

//   //now figure out distribution of deltas and make call if this means we violated

//   //first sort to get a 90% distribution length --i.e., ignore 5% on each side and take distance between those
//   std::sort( distance_deltas.begin(), distance_deltas.end() );
//   Size const ind_low_5( 1+lrint( params.local_distviol_range_*distance_deltas.size() ) );   //lower 5% -
//   Size const ind_high_5( lrint( 1.0*distance_deltas.size()*(1-params.local_distviol_range_) ) ); //upper 5%
//   Real max_extension( distance_deltas[ ind_high_5 ] - distance_deltas[ ind_low_5 ] );
//   tr.Debug << "ind_low_5 " << ind_low_5 << " ind_high_5 " << ind_high_5 << " min_delta: " << distance_deltas[ 1 ] << " max_delta "
//        << distance_deltas.back() << std::endl;
//   Size viol_count( 0 );
//   for ( RealVector::const_iterator delta_it = distance_deltas.begin(); delta_it != distance_deltas.end(); ++delta_it ) {
//    viol_count += ( *delta_it > ( max_extension + params.local_distviol_global_buffer_ ) ) ? 1 : 0;
//   }
//   tr.Debug << "peak: " << (*it)->peak_id() <<" " << (*it)->filename() << " max_extension " << max_extension << " viol_count " << viol_count << std::endl;
//   (*it)->set_eliminated_due_to_dist_violations( viol_count > (params.nr_conformers_violatable_*distance_deltas.size() ) );
//  }
// }

void DistanceScoreMover::apply( pose::Pose& pose ) {
#ifndef WIN32
	basic::ProfileThis doit( basic::NOESY_ASSIGN_DIST_APPLY );
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );

	count_decoys_++;
	SingleConstraints::const_iterator constraint_it( constraints_.begin() );
	//  PeakConstraints::const_iterator peak_constraint_it( peak_constraints_.begin() );
	active_peaks_ = 0;
	Size ct_peaks( 1 );
	if ( !used_for_calibration_ ) {
		tr.Debug << "DistanceScoreMover is not used in calibration mode " << std::endl;
	}
	for ( CrossPeakList::iterator it = cross_peaks_.begin(); it != cross_peaks_.end(); ++it, ++ct_peaks ) {
		//for each CrossPeak do the following:
		// ct_peaks .. number of peak, ct ... number of assignment / peak
		// a) iterate over all initial assignments and get effective distance (i.e., QD1 is evaluated with d^-6 averaging )
		// b1) average over all assignments (d^-6) (--> sum_dist )
		// b2) average over all selected assignments ( calibration phase --> sum_dist_filt )
		// c) cache distance in dist_buf to compute Dk in the end..

		Real dist_buf[ 2000 ];//some cheap memory buffer for distances
		runtime_assert( (*it)->n_assigned() < 2000 );
		if ( (*it)->n_assigned() == 0 ) continue;
		Size ct_assignments( 1 );
		Real sum_dist( 0.0 );
		Real sum_dist_filt( 0.0 ); //accumulate only distances where Vk > Vmin  --> used for calibration ( eq. (12)  ).
		for ( CrossPeak::iterator ait = (*it)->begin(); ait != (*it)->end(); ++ait, ++ct_assignments ) {
			runtime_assert( constraint_it != constraints_.end() );
			Real dist,invd6;
			if ( !( *constraint_it ) ) {
				dist=0;
				invd6=0;
				//tr.Trace << "PeakAssignment " << (*it)->peak_id() << " " << " assignment: " << ct_assignments << " "
				//     << (*ait)->resonance_id( 1 ) << " " << (*ait)->resonance_id( 2 )
				//     << " has invalid constraint assign 0" << std::endl;
			} else {
				//    { basic::ProfileThis doit( basic::NOESY_ASSIGN_DIST_CST_EVAL );
				dist = (*constraint_it)->dist( pose );
				//    }
				//    tr.Trace << " dist: " << dist << std::endl;
				Real invd = 1.0/dist;
				Real invd3 = invd*invd*invd;
				invd6 = invd3*invd3;
				//    tr.Trace << "invd6: " << invd6 << std::endl;
			}
			sum_dist+=invd6;
			tr.Trace << "sum_dist " << sum_dist << " dist " << dist << std::endl;
			///this resembles eq( 12 ) after filtering by peak volume --> used for calibration...
			//careful: these values are messed up when !used_for_calibration since peak_volume() can return 0 or similar ....
			sum_dist_filt += ( (*ait)->normalized_peak_volume() > params.min_volume_ )*invd6;
			dist_buf[ ct_assignments ] = dist;
			++constraint_it; //used this constraint.... we have created these constraints in prepare_scoring() in same sequence as used here
		}
		/// is upper bound violated ?
		if ( used_for_calibration_ ) {
			sum_dist=sum_dist_filt;
		}
		total_assigned_distances_ += sum_dist;
		//  tr.Trace << "sum_dist " << sum_dist << std::endl;
		if ( sum_dist > 0 )  {
			sum_dist=pow( sum_dist, -1.0/6.0 );
			++active_peaks_;
		}
		//  tr.Trace << "sum_dist " << sum_dist << std::endl;
		bool violated(  dcut_ > 0 && ( sum_dist - (*it)->distance_bound() ) > dcut_ );
		peak_violation_counts_[ ct_peaks ] += violated ? 1 : 0;
		total_violation_count_ += violated ? 1 : 0;

		if ( !used_for_calibration_ ) { ////NOTE THERE IS STILL A PROBLEM DUE TO SKIPPED CONSTRAINTS...
			basic::ProfileThis doit( basic::NOESY_ASSIGN_DIST_SET_COMPABILITY_SCORE );
			//now add to Dk in Assignments  --- formulas (6)+(7) on p.214  eta --> final_dist_power
			//d_ak,bk is distance of individual PeakAssignment (computed in for-loop above --> dist_buf[] ).
			ct_assignments = 1;
			for ( CrossPeak::iterator ait = (*it)->begin(); ait != (*it)->end(); ++ait, ++ct_assignments ) {
				//    if ( dist_buf[ ct_assignments ] == 0 ) {
				//          tr.Trace << "Crosspeak: " << (*it)->peak_id() << " assignment " << ct_assignments << " has zero dist_buf " << std::endl;
				//    }
				if ( dist_buf[ ct_assignments ] ) {
					//     tr.Trace << " dist_buf " << dist_buf[ ct_assignments ] << " over " << sum_dist << " is " << dist_buf[ ct_assignments ]/sum_dist << std::endl;
					// runtime_assert( (*ait)->decoy_compatibility() < 0.01 ); -- ah it is not 0 because we go through the different decoys and sum up!
					(*ait)->set_decoy_compatibility( (*ait)->decoy_compatibility() + pow( dist_buf[ ct_assignments ]/sum_dist, -final_dist_power_) );
				}
				if ( (*ait)->decoy_compatibility() == 0 ) {
					//     tr.Trace <<" oh drat, no compatiblity... " << (*it)->peak_id() << std::endl;
					//     tr.Trace << ct_assignments << " " << dist_buf[ct_assignments] << " " << sum_dist << " " << dist_buf[ ct_assignments ]/sum_dist << " " << pow( dist_buf[ ct_assignments ]/sum_dist, -final_dist_power_) << std::endl;
				}
			} ///... divide by M (aka count_decoys_ ) in finalize_scoring
		}
	} //loop over peaks
#endif
}

void DistanceScoreMover::finalize_scoring() const {
#ifndef WIN32
	Size ct_peaks( 1 );
	for ( CrossPeakList::iterator it = cross_peaks_.begin(); it != cross_peaks_.end(); ++it, ++ct_peaks ) {
		for ( CrossPeak::iterator ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
			(*ait)->set_decoy_compatibility( (*ait)->decoy_compatibility()/count_decoys_ );
		}
	}
#endif
}

// void DistanceScoreMover::eliminate_violated_constraints() const {
//  PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
//  Size ct_peaks( 1 );
//   for ( CrossPeakList::iterator it = cross_peaks_.begin(); it != cross_peaks_.end(); ++it, ++ct_peaks ) {
//   tr.Debug << "peak: " << (*it)->peak_id() <<" " << (*it)->filename() << " violations: " << peak_violation_counts_[ ct_peaks ] << std::endl;
//   (*it)->set_eliminated_due_to_dist_violations( peak_violation_counts_[ ct_peaks ] > (params.nr_conformers_violatable_*count_decoys_) );
//  }
// }

// core::Real DistanceScoreMover::compute_violation_percentage() const {
//  Real percent_violated( 0 );
//  Size ct_peaks( 1 );
//   for ( CrossPeakList::iterator it = cross_peaks_.begin(); it != cross_peaks_.end(); ++it, ++ct_peaks ) {
//   percent_violated+=peak_violation_counts_[ ct_peaks ];
//  }
//  return percent_violated / active_peaks_ / count_decoys_;
// }


}
}
