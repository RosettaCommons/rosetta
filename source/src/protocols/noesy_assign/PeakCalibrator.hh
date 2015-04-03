// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PeakCalibratorList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_PeakCalibrator_hh
#define INCLUDED_protocols_noesy_assign_PeakCalibrator_hh


// Unit Header
//#include <devel/NoesyAssign/PeakCalibrator.fwd.hh>

// Package Headers
// #include <devel/NoesyAssign/PeakCalibratorInfo.hh>
// #include <devel/NoesyAssign/PeakAssignment.hh>
// #include <devel/NoesyAssign/ResonanceList.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/NamedAtomID.fwd.hh>
//#include <core/chemical/AA.hh>
#include <protocols/noesy_assign/CrossPeak.hh>
#include <protocols/noesy_assign/CrossPeakList.fwd.hh>

// Utility headers
//#include <utility/exit.hh>
// #include <utility/excn/Exceptions.hh>
//#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
// #include <basic/prof.hh>
//#include <basic/Tracer.hh>
// #include <basic/options/option.hh>
// #include <basic/options/keys/abinitio.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
//#include <cstdlib>
// #include <string>
// #include <list>
#include <map>
#include <bitset>

#include <protocols/noesy_assign/PeakCalibrator.fwd.hh>
#include <utility/vector1.hh>
#include <list>

namespace protocols {
namespace noesy_assign {

class PeakCalibrator : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~PeakCalibrator();

  typedef std::bitset< MAX_TYPE > TypeCumulator;
  PeakCalibrator( int target_sign );

  virtual PeakCalibratorOP fresh_instance() = 0;

  void reset_statistics();
  void collect_target_statistics( core::Real, TypeCumulator const& );
  bool interpolate_on_statistics();

  virtual void collect_upperbound_statistics( core::Size /*peak*/, TypeCumulator const& /*types*/ ) = 0;
  //  virtual void show_statistics( std::ostream& );
  virtual void init_calibrator() {}; //to create constraints for example

  core::Real operator() ( CALIBRATION_ATOM_TYPE atom ) const { return calibration_constant( atom ); }
  core::Real calibration_constant( CALIBRATION_ATOM_TYPE type ) const { return calibration_constant_[ type ]; };
  //  void interpolate( PeakCalibrator const& cal1, PeakCalibrator const& cal2 );

  void interpolate_too_small( core::Size type );
  void interpolate_too_big( core::Size type );

  void reset_calibration_constants();

  void add_peak( CrossPeakOP peak );
  void set_new_upper_bounds();

  void do_calibration();

  void set_target_and_tolerance( core::Real target, core::Real tolerance );

  static CALIBRATION_ATOM_TYPE atom_type( core::id::NamedAtomID const& atom, core::chemical::AA aa );

  virtual void eliminate_violated_constraints() {};

protected:
  utility::vector1< CrossPeakOP > const& peaks() { return peaks_; }

private:
  core::Size max_type_direct_;
  core::Real accumulated_target_[ MAX_TYPE ];
  typedef std::list< core::Real > TargetValues;
  TargetValues target_values_[ MAX_TYPE ];
  core::Size accumulated_count_[ MAX_TYPE ];

  core::Real calibration_constant_[ MAX_TYPE ];
  core::Real calibration_constant_lows_[ MAX_TYPE ];
  core::Real calibration_constant_highs_[ MAX_TYPE ];
  core::Size max_type_;
  //  core::Real Q_backbone_;
  //core::Real Q_methyl_;
  //core::Real Q_nonmethyl_sidechain_;
  //core::Real Q_nonmethyl_beta_;
  utility::vector1< CrossPeakOP > peaks_;

  core::Real tolerance_;
  core::Real target_;
  int target_sign_; //false --> >tolerance constants too big  / true --> >tolerance if constants too small
};


class PeakCalibratorMap : public utility::pointer::ReferenceCount {
  typedef std::map < std::string, PeakCalibratorOP > CalibratorMap;
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~PeakCalibratorMap();
  PeakCalibratorMap( CrossPeakList&, PeakCalibratorOP );
  void set_new_upper_bounds();
  void do_calibration();
  void set_target_and_tolerance( core::Real target, core::Real tolerance );
  void eliminate_violated_constraints();
private:
  CalibratorMap calibrators_;
};


}
}

#endif
