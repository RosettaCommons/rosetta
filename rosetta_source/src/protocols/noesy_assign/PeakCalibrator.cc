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

// Unit Headers
#include <protocols/noesy_assign/PeakCalibrator.hh>

// Package Headers
#include <protocols/noesy_assign/CrossPeakList.hh>
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>

// Project Headers

// Utility headers
#include <basic/Tracer.hh>

//// C++ headers
#include <cmath>
#include <map>

#include <utility/vector1.hh>


static basic::Tracer tr("protocols.noesy_assign.crosspeaks");

using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

#ifdef _WIN32
#include <float.h>  // REQUIRED FOR WINDOWS
#endif

namespace protocols {
namespace noesy_assign {

PeakCalibratorMap::PeakCalibratorMap( CrossPeakList& list, PeakCalibratorOP calibrator_template ) {
	PeakCalibratorOP fresh_calibrator = calibrator_template->fresh_instance();
	for ( CrossPeakList::iterator it = list.begin(); it != list.end(); ++it ) {
		std::pair< CalibratorMap::iterator, bool > last_insert;
		last_insert = calibrators_.insert( std::make_pair( (*it)->filename(), fresh_calibrator ) );
		if ( last_insert.second ) { //true if new element was inserted
			fresh_calibrator = calibrator_template->fresh_instance();
		}
		(last_insert.first->second)->add_peak( *it );
	}
}

void PeakCalibratorMap::set_new_upper_bounds() {
	for ( CalibratorMap::iterator it=calibrators_.begin(); it!=calibrators_.end(); ++it ) {
		it->second->set_new_upper_bounds();
	}
}

void PeakCalibratorMap::do_calibration() {
	for ( CalibratorMap::iterator it=calibrators_.begin(); it!=calibrators_.end(); ++it ) {
		tr.Info << "Calibrate " << it->first << "..." << std::endl;
		it->second->do_calibration();
	}
}

void PeakCalibratorMap::set_target_and_tolerance( core::Real target, core::Real tolerance ) {
	for ( CalibratorMap::iterator it=calibrators_.begin(); it!=calibrators_.end(); ++it ) {
		it->second->set_target_and_tolerance( target, tolerance );
	}
}

void PeakCalibratorMap::eliminate_violated_constraints() {
	for ( CalibratorMap::iterator it=calibrators_.begin(); it!=calibrators_.end(); ++it ) {
		it->second->eliminate_violated_constraints();
	}
}

void PeakCalibrator::set_target_and_tolerance( core::Real target, core::Real tolerance ) {
	target_ = target;
	tolerance_ = tolerance;
}

void PeakCalibrator::reset_statistics() {
	for ( Size type=BACKBONE; type < MAX_TYPE; ++type ) {
		accumulated_count_[ type ] = 0;
		accumulated_target_[ type ] = 0;
	}
}


bool PeakCalibrator::interpolate_on_statistics() {
	bool finished = true;
	for ( Size type = BACKBONE; type < MAX_TYPE; ++type ) {
		if ( accumulated_count_[ type ] ) {
			//			tr.Debug << " acc. target: " << accumulated_target_[ type ] << " acc. count: " << accumulated_count_[ type ] << std::endl;
			core::Real average_target = accumulated_target_[ type ] / accumulated_count_[ type ];
#ifdef _WIN32
			if ( _isnan(average_target) || !_finite( average_target)) return true;  // REQUIRED FOR WINDOWS
#else
			if ( std::isnan(average_target) || std::isinf( average_target)) return true;
#endif
			if ( target_sign_* ( average_target - target_) < -tolerance_ ) {
				interpolate_too_small( type );
				finished = false;
			} else if ( target_sign_* ( average_target - target_ ) > tolerance_ ) {
				interpolate_too_big( type );
				finished = false;
			}
		}
	}
	return finished;
}

void PeakCalibrator::interpolate_too_small( core::Size type ) {
	calibration_constant_lows_[ type ] = calibration_constant_[ type ];
  calibration_constant_[ type ] = exp( 0.5*( log( calibration_constant_lows_[ type ] ) + log( calibration_constant_highs_[ type ] ) ) );
}

void PeakCalibrator::interpolate_too_big( core::Size type ) {
	calibration_constant_highs_[ type ] = calibration_constant_[ type ];
  calibration_constant_[ type ] = exp( 0.5*( log( calibration_constant_lows_[ type ] )+ log( calibration_constant_highs_[ type ] ) ) );
}

void PeakCalibrator::collect_target_statistics( core::Real target, TypeCumulator const& types ) {
	for ( core::Size type = BACKBONE; type < MAX_TYPE; ++type ) {
		if ( types.test( type ) ) {
			accumulated_count_[ type ] += 1;
			accumulated_target_[ type ] += target;
		}
	}
	//	tr.Debug << "acc. " << accumulated_target_[ BACKBONE ] << std::endl;
}

// void PeakCalibrator::reset_statistics() {
// 	for ( Size type=BACKBONE; type < MAX_TYPE; ++type ) {
// 		accumulated_count_[ type ] = 0;
// 		accumulated_dist_[ type ] = 0;
// 	}
// }

void PeakCalibrator::reset_calibration_constants() {
	for ( Size type=BACKBONE; type < MAX_TYPE; ++type ) {
		calibration_constant_[ type ] = 1e10;
		calibration_constant_lows_[ type ] = 1;
		calibration_constant_highs_[ type ] = 1e20;
	}
}

void PeakCalibrator::do_calibration() {
  bool finished = false;
  Size max_cycles = 30;

	init_calibrator();
	reset_calibration_constants();

  tr.Info << "Calibration .... for " << peaks_.size() << " crosspeaks " << std::endl;

	//	Q_backbone_ = calibration_constant_[ BACKBONE ];
	tr.Info << "value   target   BACKBONE   SIDECHAIN  METHYL "<< std::endl;
  while ( !finished && max_cycles ) {
    --max_cycles;
		reset_statistics();
		set_new_upper_bounds(); //compute statistics about upper bounds
		//		show_statistics( tr.Info );
		tr.Info << accumulated_target_[ BACKBONE ]/accumulated_count_[ BACKBONE ] << " " << target_;
		tr.Info << " " << calibration_constant_[ BACKBONE ] << " " << calibration_constant_[ SIDECHAIN ] << " " << calibration_constant_[ METHYL ] << std::endl;

		finished = interpolate_on_statistics();
		//		Q_backbone_ = calibration_constant_[ BACKBONE ];
	}
}

// void PeakCalibrator::interpolate( PeakCalibrator const& cal1, PeakCalibrator const& cal2 ) {
//   Q_backbone_ = exp( 0.5*( log( cal1.Q_backbone_) + log( cal2.Q_backbone_ ) ));
// 	Q_nonmethyl_beta_ = exp( 0.5*( log(cal1.Q_nonmethyl_beta_) + log( cal2.Q_nonmethyl_beta_) ));
//   Q_methyl_ = 3.0 * Q_backbone_;
//   Q_nonmethyl_sidechain_ = 1.5 * Q_nonmethyl_beta_;
// }

void PeakCalibrator::add_peak( CrossPeakOP peak ) {
	peaks_.push_back( peak );
}


void PeakCalibrator::set_new_upper_bounds() {
	core::Size ct( 1 );
	for ( utility::vector1< CrossPeakOP >::iterator it = peaks_.begin(); it != peaks_.end(); ++it, ++ct ) {
		TypeCumulator types;
		(*it)->calibrate( *this, types );
		core::Real dist( (*it)->distance_bound() );
		if ( dist == 0 ) continue;
		collect_upperbound_statistics( ct, types );
	}
}


CALIBRATION_ATOM_TYPE PeakCalibrator::atom_type( core::id::NamedAtomID const& atom ) {
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	if ( !params.atom_dependent_calibration_ ) return BACKBONE;
	if ( atom.atom() == "HA" || atom.atom() == "H" ) {
		return BACKBONE;
	} else if ( atom.atom().find( "Q" ) != std::string::npos ) {
		return METHYL;
	} else return SIDECHAIN;
}


}
}
