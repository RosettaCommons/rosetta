// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief protocols for folding into density
/// @detailed
/// @author Frank DiMaio

#include <protocols/rbsegment_relax/MakeStarTopology.hh>
#include <protocols/rbsegment_relax/MakeStarTopologyCreator.hh>

#include <protocols/loops/loops_main.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Remarks.hh>

#include <protocols/rbsegment_relax/util.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace rbsegment_relax {

static thread_local basic::Tracer TR( "protocols.cryst.cryst_movers" );

using namespace protocols;
using namespace core;


//////////////////
//////////////////
/// creator


std::string
MakeStarTopologyMoverCreator::keyname() const {
	return MakeStarTopologyMoverCreator::mover_name();
}

protocols::moves::MoverOP
MakeStarTopologyMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MakeStarTopologyMover );
}

std::string
MakeStarTopologyMoverCreator::mover_name() {
	return "MakeStarTopology";
}


//////////////////
//////////////////
/// mover

void MakeStarTopologyMover::apply( core::pose::Pose & pose ) {
	using namespace protocols::rbsegment_relax;

	if (restore_) {
		pose.fold_tree( *ft_restore_ );
		// ??? cutpoint variants
		protocols::loops::remove_cutpoint_variants( pose );
	} else {
		*ft_restore_ = pose.fold_tree();
		if (mode_ == "disconnected")
			setup_disconnected( pose );
		else
			//default
			setup_star_topology( pose );
	}
}

void MakeStarTopologyMover::parse_my_tag( 
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & data,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & pose )
{
	mode_ = tag->getOption<std::string>("mode", "default");

	//
	restore_ = tag->getOption<bool>("restore", false);
	tag_ = tag->getOption<std::string>("tag", "nulltag");

	// look for tag on datamap
	if( data.has( "foldtrees", tag_ ) ){
		ft_restore_ = data.get_ptr<core::kinematics::FoldTree>( "foldtrees", tag_ );
		TR << "Found foldtree " << tag_ << " on datamap" << std::endl;
	} else {
		ft_restore_ = core::kinematics::FoldTreeOP( new core::kinematics::FoldTree( pose.fold_tree() ) );
			// make a copy of current pose ... if restore is called before 'set' this could be oddly behaved
		data.add( "foldtrees", tag_, ft_restore_ );
		TR << "Adding foldtrees " << tag_ << " to datamap" << std::endl;
	}
}

}
}
