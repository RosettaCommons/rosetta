// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file IterativeAbrelax
/// @brief iterative protocol starting with abinitio and getting progressively more concerned with full-atom relaxed structures
/// @detailed
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/abinitio/IterativeCentroid.hh>
// AUTO-REMOVED #include <protocols/jd2/archive/ArchiveManager.hh>

// Package Headers

// Project Headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/io/silent/SilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>

// AUTO-REMOVED #include <core/fragment/ConstantLengthFragSet.hh>
// AUTO-REMOVED #include <core/fragment/FragmentIO.hh>
// AUTO-REMOVED #include <core/fragment/util.hh>

// AUTO-REMOVED #include <protocols/toolbox/DecoySetEvaluation.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>

// #include <core/scoring/ScoreType.hh>
// //#include <core/kinematics/MoveMap.hh>
// #include <core/types.hh>
// #include <core/scoring/rms_util.hh> //for ConvergenceCheck
// //#include <core/pack/task/PackerTask.fwd.hh>
// #include <core/scoring/constraints/ConstraintSet.hh>
// //only needed because of temporary output_debug_structure ...

// #include <core/io/silent/SilentStructFactory.hh>
// #include <basic/options/keys/out.OptionKeys.gen.hh>
// #include <protocols/jd2/util.hh>
// //#include <protocols/moves/Mover.hh>
// //#include <protocols/moves/MoverContainer.hh>
// #include <protocols/moves/TrialMover.hh>
// #include <protocols/moves/RepeatMover.hh>
// //#include <protocols/moves/WhileMover.hh>
// #include <protocols/checkpoint/Checkpoint.hh>

// ObjexxFCL Headers
// //#include <ObjexxFCL/string.functions.hh>

// Utility headers
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/file/FileName.hh>

// #include <utility/exit.hh>
// #include <utility/vector1.fwd.hh>
// #include <utility/pointer/ReferenceCount.hh>
// #include <utility/file/file_sys_util.hh>
// #include <numeric/numeric.functions.hh>
// #include <basic/prof.hh>
#include <basic/Tracer.hh>
// #include <basic/options/option.hh>

// Option Headers
// AUTO-REMOVED #include <basic/options/keys/abinitio.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>
// //#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
// AUTO-REMOVED #include <cstdlib>
#include <string>
// AUTO-REMOVED #include <ctime>
// AUTO-REMOVED #include <iterator>

// Utility headers
// AUTO-REMOVED #include <basic/options/option_macros.hh>

#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>


static basic::Tracer tr("protocols.iterative");

using core::Real;
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;


namespace protocols {
namespace abinitio {
using namespace jd2::archive;

void IterativeCentroid::gen_evaluation_output( jd2::archive::Batch& batch, bool fullatom ) {
	if ( fullatom ) fullatom_pool_ptr_->gen_evaluation_output( batch, fullatom );
	else Parent::gen_evaluation_output( batch, fullatom );
}

} //abinitio
} //protocols
