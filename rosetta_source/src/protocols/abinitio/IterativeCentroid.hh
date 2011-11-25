// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file AbrelaxMover
/// @brief  this class will be handled to a SampleProtocol as a control instance
/// @detailed responsibilities:
///           know which chainbreaks to penalize and close
///           know which jumps to use during sampling, which (if any) to keep after loop-closing
///           supply a JumpMover if jumps should be moved
///           supply a MoveMap
///           supply a "StrictMoveMap": the protocol should not move anything that is dissallowed in strict_movemap(),
///                      it should try to move just stuff in movemap()
/// should this class also know how to ramp score terms ?
/// handle the titration of constraints ?
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_IterativeCentroid_hh
#define INCLUDED_protocols_abinitio_IterativeCentroid_hh

// Unit Headers
//#include <protocols/abinitio/IterativeCentroid.fwd.hh>

// Package Headers
//#include <protocols/jd2/archive/ArchiveBase.hh>
//#include <protocols/jd2/archive/EvaluatedArchive.hh>
#include <protocols/abinitio/IterativeBase.hh>
#include <protocols/abinitio/IterativeFullatom.hh>

// Project Headers
// AUTO-REMOVED #include <protocols/abinitio/PairingStatistics.hh>

//#include <core/types.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>  // Needed so PyRosetta can generate copy constructor

//#include <protocols/loops/Loops.hh>
#include <protocols/jd2/archive/ArchiveManager.fwd.hh>

// ObjexxFCL Headers
//#include <ObjexxFCL/FArray1D.hh>
//#include <ObjexxFCL/FArray2D.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
//#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
//#include <cstdlib>
//#include <string>

#include <utility/vector1.hh>

//Auto Headers
//#include <core/scoring/constraints/Constraint.hh>


namespace protocols {
namespace abinitio {

class IterativeCentroid : public IterativeBase {
	typedef IterativeBase Parent;
public:
  IterativeCentroid(IterativeFullatom* fullatom_pool_ptr ) :
		IterativeBase( "centroid_pool" ),
		fullatom_pool_ptr_( fullatom_pool_ptr ) {};

	virtual void gen_evaluation_output( jd2::archive::Batch& batch, bool fullatom = false );

private:
	IterativeFullatom* fullatom_pool_ptr_;
};



}
}

#endif
