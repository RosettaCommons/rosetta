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
#include <utility/numbers.hh>

static THREAD_LOCAL basic::Tracer tr( "protocols.noesy_assign.crosspeaks" );

using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

#ifdef _WIN32
#include <float.h>  // REQUIRED FOR WINDOWS
#endif

const char* CALIBRATOR_TYPE_NAMES[]={"NONE","BACKBONE","BETA_NON_METHYL","METHYL","SIDECHAIN" };

namespace protocols {
namespace noesy_assign {

/// @details Auto-generated virtual destructor
PeakCalibratorMap::~PeakCalibratorMap() = default;

/// @details Auto-generated virtual destructor
PeakCalibrator::~PeakCalibrator() = default;

PeakCalibratorMap::PeakCalibratorMap( CrossPeakList& list, PeakCalibratorOP calibrator_template ) {
	PeakCalibratorOP fresh_calibrator = calibrator_template->fresh_instance();
	for ( auto it = list.begin(); it != list.end(); ++it ) {
		std::pair< CalibratorMap::iterator, bool > last_insert;
		last_insert = calibrators_.insert( std::make_pair( (*it)->filename(), fresh_calibrator ) );
		if ( last_insert.second ) { //true if new element was inserted
			fresh_calibrator = calibrator_template->fresh_instance();
		}
		(last_insert.first->second)->add_peak( *it );
	}
	for ( auto & calibrator : calibrators_ ) {
		calibrator.second->init_calibrator();
	}
}

void PeakCalibratorMap::set_new_upper_bounds() {
	for ( auto & calibrator : calibrators_ ) {
		calibrator.second->set_new_upper_bounds();
	}
}

void PeakCalibratorMap::do_calibration() {
	for ( auto & calibrator : calibrators_ ) {
		tr.Info << "Calibrate " << calibrator.first << "..." << std::endl;
		calibrator.second->do_calibration();
	}
}

void PeakCalibratorMap::set_target_and_tolerance( core::Real target, core::Real tolerance ) {
	for ( auto & calibrator : calibrators_ ) {
		calibrator.second->set_target_and_tolerance( target, tolerance );
	}
}

void PeakCalibratorMap::eliminate_violated_constraints() {
	for ( auto & calibrator : calibrators_ ) {
		calibrator.second->eliminate_violated_constraints();
	}
}

PeakCalibrator::PeakCalibrator( int target_sign )
: max_type_direct_( BETA_NON_METHYL + 1 ),
	target_sign_( target_sign )
{
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	if ( !params.atom_dependent_calibration_ ) {
		max_type_ = BACKBONE + 1;
		max_type_direct_ = BACKBONE + 1;
	} else {
		max_type_ = MAX_TYPE;
	}
}


void PeakCalibrator::set_target_and_tolerance( core::Real target, core::Real tolerance ) {
	target_ = target;
	tolerance_ = tolerance;
}

void PeakCalibrator::reset_statistics() {
	for ( Size type=BACKBONE; type < max_type_; ++type ) {
		accumulated_count_[ type ] = 0;
		accumulated_target_[ type ] = 0;
		target_values_[ type ].clear();
	}
}


bool PeakCalibrator::interpolate_on_statistics() {
	noesy_assign::PeakAssignmentParameters& params( *noesy_assign::PeakAssignmentParameters::get_nonconst_instance() );
	bool const use_median( params.calibration_use_median_ );
	bool finished = true;
	for ( Size type = BACKBONE; type < max_type_direct_; ++type ) {
		if ( accumulated_count_[ type ] ) {
			//   tr.Debug << " acc. target: " << accumulated_target_[ type ] << " acc. count: " << accumulated_count_[ type ] << std::endl;
			core::Real average_target = accumulated_target_[ type ] / accumulated_count_[ type ];
			core::Real median;
			if ( use_median ) {
				target_values_[ type ].sort();
				TargetValues::const_iterator it = target_values_[type].begin();
				for ( Size i=1; i<= target_values_[ type ].size()/2; ++i ) {
					++it;
				}
				median = *it;
				average_target = median;
				accumulated_target_[ type ] = median * accumulated_count_[ type ];
			}
			if ( utility::isnan(average_target) || utility::isinf(average_target) ) continue;
			if ( target_sign_* ( average_target - target_) < -tolerance_ ) {
				interpolate_too_small( type );
				finished = false;
			} else if ( target_sign_* ( average_target - target_ ) > tolerance_ ) {
				interpolate_too_big( type );
				finished = false;
			}
		}
	}
	calibration_constant_[ METHYL ] = 3.0 * calibration_constant_[ BACKBONE ];
	calibration_constant_[ SIDECHAIN ] = 1.5 * calibration_constant_[ BETA_NON_METHYL ];
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
	for ( core::Size type = BACKBONE; type < max_type_; ++type ) {
		if ( types.test( type ) ) {
			accumulated_count_[ type ] += 1;
			accumulated_target_[ type ] += target;
			target_values_[ type ].push_back( target );
		}
	}
	// tr.Debug << "acc. " << accumulated_target_[ BACKBONE ] << std::endl;
}

// void PeakCalibrator::reset_statistics() {
//  for ( Size type=BACKBONE; type < max_type_; ++type ) {
//   accumulated_count_[ type ] = 0;
//   accumulated_dist_[ type ] = 0;
//  }
// }

void PeakCalibrator::reset_calibration_constants() {
	for ( Size type=BACKBONE; type < max_type_; ++type ) {
		calibration_constant_[ type ] = 1e10;
		calibration_constant_lows_[ type ] = 1;
		calibration_constant_highs_[ type ] = 1e20;
	}
}

void PeakCalibrator::do_calibration() {
	bool finished = false;
	Size max_cycles = 50;

	reset_calibration_constants();

	tr.Info << "Calibration .... for " << peaks_.size() << " crosspeaks " << std::endl;

	// Q_backbone_ = calibration_constant_[ BACKBONE ];
	tr.Info << "value   target ";
	for ( core::Size type=BACKBONE; type< (Size) max_type_; type++ ) {
		tr.Info << " " << CALIBRATOR_TYPE_NAMES[type];
	}
	tr.Info << std::endl;
	while ( !finished && max_cycles ) {
		--max_cycles;
		reset_statistics();
		set_new_upper_bounds(); //compute statistics about upper bounds
		//  show_statistics( tr.Info );
		tr.Info << target_;
		tr.Info << " ";
		for ( core::Size type=BACKBONE; type< (Size) max_type_; type++ ) {
			if ( accumulated_count_[ type ] ) {
				tr.Info << " " << accumulated_target_[ type ] / accumulated_count_[ type ];
			} else tr.Info << " -1    ";
			tr.Info << " " << calibration_constant_[ type ];
		}
		tr.Info << std::endl;
		finished = interpolate_on_statistics();
	}
}

// void PeakCalibrator::interpolate( PeakCalibrator const& cal1, PeakCalibrator const& cal2 ) {
//   Q_backbone_ = exp( 0.5*( log( cal1.Q_backbone_) + log( cal2.Q_backbone_ ) ));
//  Q_nonmethyl_beta_ = exp( 0.5*( log(cal1.Q_nonmethyl_beta_) + log( cal2.Q_nonmethyl_beta_) ));
//   Q_methyl_ = 3.0 * Q_backbone_;
//   Q_nonmethyl_sidechain_ = 1.5 * Q_nonmethyl_beta_;
// }

void PeakCalibrator::add_peak( CrossPeakOP peak ) {
	peaks_.push_back( peak );
}

void PeakCalibrator::set_new_upper_bounds() {
	core::Size ct( 1 );
	for ( auto it = peaks_.begin(); it != peaks_.end(); ++it, ++ct ) {
		TypeCumulator types;
		(*it)->calibrate( *this, types );
		core::Real dist( (*it)->distance_bound() );
		//  if ( dist <= 0.01 ) continue;
		if ( dist == 0 ) continue;
		collect_upperbound_statistics( ct, types );
	}
}


CALIBRATION_ATOM_TYPE PeakCalibrator::atom_type( core::id::NamedAtomID const& atom, core::chemical::AA aa ) {
	using namespace core::chemical; //for AA
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	if ( !params.atom_dependent_calibration_ ) return BACKBONE;
	std::string const& name = atom.atom();
	if ( name == "HA" || name == "H" || name == "1H" || name == "2H" || name == "QA" || name =="1HA" || name=="2HA" ) {
		return BACKBONE;
	}
	if ( name.find("B") != std::string::npos ) {
		if ( aa==aa_ala ) return METHYL;
		return BETA_NON_METHYL;
	}
	if ( name.find("G") != std::string::npos ) {
		if ( aa==aa_val || aa==aa_thr ) return METHYL;
	}
	if ( name.find("G2") != std::string::npos && aa==aa_ile ) {
		return METHYL;
	}
	if ( name.find("D") != std::string::npos && ( aa==aa_leu || aa==aa_ile ) ) {
		return METHYL;
	}
	if ( name.find("E") != std::string::npos && ( aa==aa_met ) ) {
		return METHYL;
	}
	return SIDECHAIN;
}


}
}
