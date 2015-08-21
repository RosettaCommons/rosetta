// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_build/LoopmodelWrapper.hh>
#include <protocols/loop_build/LoopmodelWrapperCreator.hh>
#include <protocols/loop_build/LoopBuildMover.hh>

// Protocol headers
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/loops/loops_main.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/vector1.hh>

// RosettaScripts headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace loop_build {

using namespace std;
using core::Size;
using core::Real;

protocols::moves::MoverOP LoopmodelWrapperCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopmodelWrapper );
}

std::string LoopmodelWrapperCreator::keyname() const {
	return "LoopmodelWrapper";
}

protocols::moves::MoverOP LoopmodelWrapper::clone() const {
	return protocols::moves::MoverOP( new LoopmodelWrapper );
}

void LoopmodelWrapper::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &) {

	using namespace basic::options;
	using utility::options::StringVectorOption;
	using utility::options::BooleanOption;

	// The loopmodel application makes it impossible to specify which loop file
	// to use without setting a user-specified command-line option.  (Really
	// impossible; there's no quick hack that would loosen this restriction.)
	// Since I need to be able to specify a loops file from an XML script, I had
	// to hack the options system.  This is a pretty brutal hack for a number of
	// reasons, not least of which because it breaks the "don't modify the
	// options data structure" convention.  The consolation is that loopmodel
	// will be obsolete soon.

	string loops_file = tag->getOption<string>("loops_file");
	StringVectorOption & loop_files = option[OptionKeys::loops::loop_file];
	loop_files.clear();
	loop_files.value(loops_file);

	// I also want to be able to set run:test_cycles from a rosetta script.  This
	// is not as critical, but still useful.

	bool go_fast = tag->getOption<bool>("fast");
	BooleanOption & test_cycles = option[OptionKeys::run::test_cycles];
	test_cycles.value(go_fast);
}

void LoopmodelWrapper::apply(core::pose::Pose & pose) {
	// This entire method was copied from LoopBuild_main.  The only difference is
	// that this method calls apply() at the end, while LoopBuild_main calls
	// jd2::go().  I couldn't think of any way to avoid copying this code.

	basic::Tracer TR("protocols.loop_build.LoopBuild");

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::chemical;
	using namespace core::id;
	using namespace jobdist;

	std::string remodel            ( option[ OptionKeys::loops::remodel ]() );
	std::string const intermedrelax( option[ OptionKeys::loops::intermedrelax ]() );
	std::string const refine       ( option[ OptionKeys::loops::refine ]() );
	std::string const relax        ( option[ OptionKeys::loops::relax ]() );

	TR << "==== Loop protocol: ================================================="
		<< std::endl;
	TR << " remodel        " << remodel        << std::endl;
	TR << " intermedrelax  " << intermedrelax  << std::endl;
	TR << " refine         " << refine         << std::endl;
	TR << " relax          " << relax          << std::endl;

	utility::vector1< core::fragment::FragSetOP > frag_libs;
	if ( option[ OptionKeys::loops::frag_files ].user() ) {
		loops::read_loop_fragments( frag_libs );
	}

	comparative_modeling::LoopRelaxMover looprelax_mover;
	looprelax_mover.frag_libs( frag_libs );
	looprelax_mover.relax( relax );
	looprelax_mover.refine( refine );
	looprelax_mover.remodel( remodel );
	looprelax_mover.intermedrelax( intermedrelax );

	LoopBuildMoverOP loopbuild_mover( new protocols::loop_build::LoopBuildMover(looprelax_mover) );

	loopbuild_mover->apply(pose);
}

}
}
