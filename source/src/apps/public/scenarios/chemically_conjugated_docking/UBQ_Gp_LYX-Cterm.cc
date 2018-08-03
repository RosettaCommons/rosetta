// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/public/scenarios/chemically_conjugated_docking/UBQ_Gp_LYX-Cterm.cc
/// @brief  this application is a one-shot for modeling a ubiquitinated G-protein; this version uses a native lysine to ubiquitin linkage via LYX
/// @author Steven Lewis and Hope Anderson

// Project Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

//movers
#include <protocols/chemically_conjugated_docking/UBQ_GTPaseMover.hh>

#include <protocols/moves/MoverFactory.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/SingletonBase.hh>
#include <basic/prof.hh>

// C++ headers
#include <string>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/chemically_conjugated_docking.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

//tracers
using basic::Error;
using basic::Warning;
static basic::Tracer TR("apps.public.scenarios.chemically_conjugated_docking.UBQ_Gp_LYX-Cterm");

int main( int argc, char* argv[] )
{
	try {
		//initialize options
		devel::init(argc, argv);

		// got rid of this so that the default value for GTPasepdb is not attempted to be used.
		//basic::prof_reset();

		protocols::chemically_conjugated_docking::UBQ_GTPaseMoverOP ubq_mov(new protocols::chemically_conjugated_docking::UBQ_GTPaseMover());

		// setting ubq file path
		std::string const & ubqFile = basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::UBQpdb].value();
		ubq_mov->set_UBQpdb(ubqFile);

		// setting residue selector from lysine index
		core::Size const lysNum = basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::GTPase_residue].value();

		//core::select::residue_selector::ResidueSelectorOP res_op(new core::select::residue_selector::ResidueIndexSelector(lysNum));
		ubq_mov->make_index_selector(lysNum);

		//ubq_mov->set_selector(utility::pointer::static_pointer_cast< core::select::residue_selector::ResidueIndexSelector >(res_op));

		// determining if there are extra bodies
		if ( basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::extra_bodies].user() ) {
			ubq_mov->set_extra_bodies(true);
		}

		// setting n_tail_res
		core::Size const nTailRes = basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::n_tail_res].value();
		ubq_mov->set_n_tail_res(nTailRes);

		// setting filters
		core::Real const scorefilter = basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::scorefilter].value();
		ubq_mov->set_scorefilter(scorefilter);

		core::Real const SASAfilter = basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::SASAfilter].value();
		ubq_mov->set_SASAfilter(SASAfilter);

		protocols::jd2::JobDistributor::get_instance()->go(ubq_mov);

		//basic::prof_show();
		TR << "************************d**o**n**e**************************************" << std::endl;
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
