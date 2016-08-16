// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PeakCalibratorList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_StructureDependentPeakCalibrator_hh
#define INCLUDED_protocols_noesy_assign_StructureDependentPeakCalibrator_hh


// Unit Header
#include <protocols/noesy_assign/PeakCalibrator.hh>


// Package Headers
#include <protocols/noesy_assign/CrossPeak.hh>
#include <protocols/noesy_assign/CrossPeakList.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

//// C++ headers
//#include <map>
//#include <bitset>

namespace protocols {
namespace noesy_assign {

class StructureDependentPeakCalibrator : public PeakCalibrator {
public:
	typedef utility::vector1< core::pose::PoseOP > PoseVector;

	StructureDependentPeakCalibrator( PoseVector const& structures, core::Real dcalibrate )
	:    PeakCalibrator( -1 /*reverse sign*/ ),
		structures_( structures ),
		dcalibrate_( dcalibrate )
	{}

	virtual PeakCalibratorOP fresh_instance() {
		return PeakCalibratorOP( new StructureDependentPeakCalibrator( structures_, dcalibrate_ ) );
	}


	//  virtual void reset_statistics();
	//  virtual bool interpolate_on_statistics();
	virtual void collect_upperbound_statistics( core::Size /*peak*/, TypeCumulator const& /*types*/ );
	virtual void init_calibrator(); //to create constraints for example
	void generate_constraints();

	// this is a service that has nothing to do with calibration,
	// however, after calibration we have all the  constraints already generated, that are necessary to do this task
	void eliminate_violated_constraints();

private:
	//  core::Real accumulated_viol_percentage_[ MAX_TYPE ];
	//  core::Size accumulated_count_[ MAX_TYPE ];
	PoseVector structures_;
	utility::vector1< core::scoring::constraints::ConstraintOP > constraints_;
	core::Real dcalibrate_;
};


}
}

#endif
