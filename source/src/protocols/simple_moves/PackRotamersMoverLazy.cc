// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Hetu Kamisetty

// Unit headers
#include <protocols/simple_moves/PackRotamersMoverLazy.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/PackRotamersMoverLazyCreator.hh>


#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

// option key includes

#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>


namespace protocols {
namespace simple_moves {

using namespace core;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;

using basic::Warning;
using basic::t_warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.PackRotamersMoverLazy" );

PackRotamersMoverLazy::PackRotamersMoverLazy(
	ScoreFunctionCOP scorefxn,
	PackerTaskCOP task,
	core::Size nloop
) : protocols::simple_moves::PackRotamersMover(scorefxn,task,nloop)
{}

PackRotamersMoverLazy::PackRotamersMoverLazy() : protocols::simple_moves::PackRotamersMover( std::string("PackRotamersMoverLazy") )
{}

PackRotamersMoverLazy::~PackRotamersMoverLazy(){}

void
PackRotamersMoverLazy::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & pose
)
{
	if ( tag->getName() != "PackRotamersMoverLazy" ) {
		TR(t_warning) << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( tag->hasOption("nloop") ) {
		nloop(tag->getOption<Size>("nloop",1));
		runtime_assert( nloop() > 0 );
	}
	parse_score_function( tag, datamap, filters, movers, pose );
	parse_task_operations( tag, datamap, filters, movers, pose );
}

/// PackRotamersMoverLazy


std::string
PackRotamersMoverLazyCreator::mover_name()
{
	return "PackRotamersMoverLazy";
}

void
PackRotamersMoverLazy::call_setup( Pose & pose)
{
	this->setup(pose);//to bypass the protected status of setup()
}
void
PackRotamersMoverLazy::apply_to_rotpack( Pose & pose , utility::vector0< int > rot_to_pack)
{
	//assume setup has been called once.

	core::PackerEnergy best_energy(0.);
	Pose best_pose;
	best_pose = pose;
	for ( Size run(1); run <= nloop(); ++run ) {
		// run SimAnnealer
		core::PackerEnergy packer_energy( this->run( pose, rot_to_pack) );
		//  Real const score( scorefxn_( pose ) ); another option for deciding which is the 'best' result
		if ( run == 1 || packer_energy < best_energy ) {
			best_pose = pose;
			best_energy = packer_energy;
		}
	}
	if ( nloop() > 1 ) pose = best_pose;

	ScoreFunctionCOP scorefxn_ = score_function();
	(*scorefxn_)(pose);
}
/*
core::PackerEnergy PackRotamersMoverLazy::run_with_ig( Pose & pose, utility::vector0< int > rot_to_pack, InteractionGraphBaseOP ig) const
{
return pack_rotamers_run( pose, task(), rotamer_sets(), ig, rot_to_pack );
}
*/

}//moves
}//protocols
