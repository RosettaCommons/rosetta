// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PeakCalibratorList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_StructureIndependentPeakCalibrator_hh
#define INCLUDED_protocols_noesy_assign_StructureIndependentPeakCalibrator_hh


// Unit Header
//#include <devel/NoesyAssign/StructureIndependentPeakCalibrator.fwd.hh>
#include <protocols/noesy_assign/PeakCalibrator.hh>

// Package Headers
// #include <devel/NoesyAssign/StructureIndependentPeakCalibratorInfo.hh>
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

#include <utility/vector1.hh>


namespace protocols {
namespace noesy_assign {

class StructureIndependentPeakCalibrator : public PeakCalibrator {
public:

  StructureIndependentPeakCalibrator() : PeakCalibrator( 1 /*normal sign*/ ) {};
  virtual PeakCalibratorOP fresh_instance() {
    return PeakCalibratorOP( new StructureIndependentPeakCalibrator() );
  }

  virtual void collect_upperbound_statistics( core::Size  /*peak*/, TypeCumulator const& /*types*/ );
private:
};





}
}

#endif
