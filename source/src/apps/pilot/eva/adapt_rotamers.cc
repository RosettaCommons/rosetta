// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @brief
/// takes 2 input poses and allows one pose to adopt a given residue of the other pose at a given residue.
/// currently only uses rosetta numbering
/// @author Eva

#include <core/types.hh>
#include <core/conformation/Residue.hh>
//#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.hh>

#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/util.hh>
#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>

// Unit headers
#include <protocols/protein_interface_design/design_utils.hh>

// C++ headers
#include <iostream>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


using basic::T;
using basic::Error;
using basic::Warning;
static numeric::random::RandomGenerator RG( 25051978 ); // <- Magic number, do not change it!!!

static basic::Tracer TR( "apps.pilot.eva.adapt_rotamers" );

using namespace core;
using namespace protocols::protein_interface_design;

using namespace basic::options;

namespace adapt_rotamers
{
	StringOptionKey ori( "adapt_rotamers:ori" );
	StringOptionKey adapting( "adapt_rotamers:adapting" );
//	StringOptionKey residue ("adapt_rotamers:residue");
//	StringOptionKey chain1 ("adapt_rotamer:wt chain");
//	StringOptionKey chain2 ("adapt_rotamers:adapting chain");

	}

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

	option.add( adapt_rotamers::ori, "The ori file name" );
	option.add( adapt_rotamers::adapting, "The file name of the pdb that adapts certain rotamers");
//	option.add( adapt_rotamers::residue, "position of that should be taken over");
//	option.add( adapt_rotamers::chain1, "wt chain ID");
//	option.add( adapt_rotamers::chain2. "chain ID of the structurs that is supposed to take over the rotamer");

	// setup random numbers and options
	devel::init(argc, argv);

	// read the pose
	pose::Pose pose_ori;
	pose::Pose pose_adpt;
	core::import_pose::pose_from_pdb( pose_ori, option[ adapt_rotamers::ori ] );
	core::import_pose::pose_from_pdb( pose_adpt, option[ adapt_rotamers::adapting ] );

	typedef conformation::Residue Residue;

	for( core::Size res=1; res<=pose_adpt.total_residue(); ++res ) {

		//	TR<<"this is res "<<res <<" this is pose_adpt.residue (res)" <<pose_adpt.residue(res) << "\n";

		if (res==16)
		{
			pose_adpt.replace_residue( res, pose_ori.residue( 15 ), true/*orient_backbone*/ );
			//update_residue_coordinates(16, true);
			pose_adpt.conformation().update_polymeric_connection( res );
		}
		else if (res==18)
		{
			pose_adpt.replace_residue( res, pose_ori.residue( 17 ), true/*orient_backbone*/ );
			pose_adpt.conformation().update_polymeric_connection( res );
		}
		else if (res==19)
		{
			pose_adpt.replace_residue( res, pose_ori.residue( 18 ), true/*orient_backbone*/ );
			pose_adpt.conformation().update_polymeric_connection( res );
		}
		else if (res==192)
		{
			pose_adpt.replace_residue( res, pose_ori.residue( 191 ), true/*orient_backbone*/ );
			pose_adpt.conformation().update_polymeric_connection( res );
	}
		else if (res==200)
		{
			pose_adpt.replace_residue( res, pose_ori.residue( 199 ), true/*orient_backbone*/ );
			pose_adpt.conformation().update_polymeric_connection( res );
		}

	}
	TR<<"total residue number adpt "<< pose_adpt.total_residue() <<" and ori : " << pose_ori.total_residue()<<std::endl;

	pose_adpt.dump_pdb( "adapted.pdb");
	pose_ori.dump_pdb( "ori.pdb");
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
}

