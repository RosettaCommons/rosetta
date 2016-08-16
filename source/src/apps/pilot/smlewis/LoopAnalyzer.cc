// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/smlewis/LoopAnalyzer.cc
/// @brief Q&D protocol to run LoopAnalyzerMover as protocol
/// @author Steven Lewis

// Unit Headers
#include <protocols/analysis/LoopAnalyzerMover.hh>

// Project Headers
#include <protocols/loops/Loops.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.smlewis.LoopAnalyzerMover" );

class hackLAMMover : public protocols::moves::Mover {
public:
	hackLAMMover(protocols::loops::Loops const & loops) : LAM(new protocols::analysis::LoopAnalyzerMover(loops, true)){}

	virtual
	void
	apply(core::pose::Pose & pose ){
		LAM->apply(pose);
		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY); //hacky way to prevent output
		return;
	}

	virtual
	std::string
	get_name() const { return "hackLAMMover"; }

	protocols::analysis::LoopAnalyzerMoverOP LAM;
};

int
main( int argc, char* argv[] )
{

	try {

	devel::init(argc, argv);

	// read loops file
	protocols::loops::Loops loops( basic::options::option[ basic::options::OptionKeys::loops::loop_file ].value()[1] );

	TR << "initial loops " << loops << std::flush;
	TR << "Read in files" << std::endl;

	protocols::jd2::JobDistributor::get_instance()->go(protocols::moves::MoverOP(new hackLAMMover(loops)));

	TR << "************************d**o**n**e**************************************" << std::endl;

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
