// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/public/revert_design_to_native.cc
/// @brief For protein-interface design, reverts residues at the interface to their wildtype identities and if binding energy
/// isn't adversely affected applies the change. Produces a report at the end of the run of all of the changes.

/// Documentation is available as part of RosettaScripts in the RosettaCommons documentation

// Project headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/Pose.hh>

#include <basic/options/keys/OptionKeys.hh>
#include <devel/init.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>

#include <basic/Tracer.hh>

// Unit headers
#include <protocols/protein_interface_design/design_utils.hh>

// C++ headers
#include <iostream>

// option key includes

#include <basic/options/option.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace protocols::protein_interface_design;

using namespace basic::options;

namespace revert_app
{
	StringOptionKey wt( "revert_app:wt" );
	StringOptionKey design( "revert_app:design" );
	IntegerOptionKey ddg_cycles( "revert_app:ddg_cycles" );
	RealOptionKey threshold( "revert_app:threshold" );
	BooleanOptionKey post_repack( "revert_app:post_repack" );
}

static thread_local basic::Tracer TR( "apps.public.rosetta_scripts.revert_design_to_native" );

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
	using namespace core::pose;
	using namespace core::scoring;
	using namespace conformation;
	using namespace protocols::protein_interface_design;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add( revert_app::wt, "The wt file name" );
	option.add( revert_app::design, "The design file name");
	option.add( revert_app::ddg_cycles, "how many ddg cycles to compute (more leads to higher numerical stability)" ).def( 5 );
	option.add( revert_app::threshold, "ddg threshold for acceptance of reversion" ).def( 0.5 );
	option.add( revert_app::post_repack, "attempt repacking of the reverted structure prior to output" ).def( 0 );

	// setup random numbers and options
	devel::init(argc, argv);

	// read the pose
	pose::Pose pose_wt, pose_des;
	std::string const wt_fname(  option[ revert_app::wt ] );
	std::string const des_fname(  option[ revert_app::design ] );
	core::import_pose::pose_from_pdb( pose_wt, wt_fname );
	core::import_pose::pose_from_pdb( pose_des, des_fname );
	pose::Pose pose_ref(pose_des);

	ScoreFunctionOP scorefxn( core::scoring::get_score_function() ); // defaults to sc12
	core::Real const ddg_thres( option[ revert_app::threshold ] );
	Revert rev( scorefxn, ddg_thres, option[ revert_app::ddg_cycles ] );
	rev.apply( pose_wt, pose_des );

	if( option[ revert_app::post_repack ] ) {
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose_des ) );
		task->initialize_from_command_line().or_include_current( true );
		task->restrict_to_repacking();
		utility::vector1< bool > packable(pose_des.n_residue(),0);
		//Only repack positions which changed identity
		for(core::Size ii(1); ii <= pose_des.n_residue(); ++ii) {
			if( pose_des.aa(ii) != pose_ref.aa(ii) ) { packable[ii] = 1; }
		}
		task->restrict_to_residues(packable);
		runtime_assert(!task->design_any());
		utility::vector1< bool > torepack(task->repacking_residues());
		TR.Debug << "Post reversion repacking of residues ";
		for(core::Size ii(1); ii <= torepack.size(); ++ii) {
			if(torepack[ii]) { TR.Debug << ii << "+";}
		}
		TR.Debug << std::endl;
		core::pack::pack_rotamers( pose_des, *scorefxn, task );
	}

	pose_des.dump_scored_pdb( des_fname+".revert.pdb", *scorefxn );//, *scorefxn );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

