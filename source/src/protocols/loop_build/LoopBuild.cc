// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief demo program for implementing loop relax + FA relax
/// @author Srivatsan Raman
/// @author James Thompson
/// @author Mike Tyka
/// @author Daniel J. Mandell
/// @author Roland A. Pache, PhD

//unit headers
#include <protocols/loop_build/LoopBuild.hh>
#include <protocols/loop_build/LoopBuildMover.hh>

// basic headers
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <basic/Tracer.hh>


#include <core/chemical/ChemicalManager.hh>
#include <core/fragment/FragSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/types.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>

#include <numeric/random/random.hh>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace loop_build {

int
LoopBuild_main( bool  ) {
	basic::Tracer TR( "protocols.loop_build.LoopBuild" );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::chemical;
	using namespace core::id;
	using namespace jobdist;

	//using namespace basic::resource_manager;

	std::string remodel            ( option[ OptionKeys::loops::remodel ]() );
	std::string const intermedrelax( option[ OptionKeys::loops::intermedrelax ]() );
	std::string const refine       ( option[ OptionKeys::loops::refine ]() );
	std::string const relax        ( option[ OptionKeys::loops::relax ]() );
	//bool const keep_time      ( option[ OptionKeys::loops::timer ]() );

	TR << "==== Loop protocol: ================================================="
		<< std::endl;
	TR << " remodel        " << remodel        << std::endl;
	TR << " intermedrelax  " << intermedrelax  << std::endl;
	TR << " refine         " << refine         << std::endl;
	TR << " relax          " << relax          << std::endl;


	// fragment initialization
	utility::vector1< core::fragment::FragSetOP > frag_libs;
	if ( option[ OptionKeys::loops::frag_files ].user()) {
        loops::read_loop_fragments( frag_libs );
	}

	//setup of looprelax_mover
	comparative_modeling::LoopRelaxMover looprelax_mover;
	looprelax_mover.frag_libs( frag_libs );
	looprelax_mover.relax( relax );
	looprelax_mover.refine( refine );
	looprelax_mover.remodel( remodel );
	looprelax_mover.intermedrelax( intermedrelax );

	LoopBuildMoverOP loopbuild_mover( new protocols::loop_build::LoopBuildMover(looprelax_mover) );


	protocols::jd2::JobDistributor::get_instance()->go(loopbuild_mover);

	return 0;
} // Looprelax_main

} // namespace loop_build
} // namespace protocols
