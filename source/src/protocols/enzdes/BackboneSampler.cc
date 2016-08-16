// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (C) 199x-27 University of Washington
// (C) 199x-27 University of California Santa Cruz
// (C) 199x-27 University of California San Francisco
// (C) 199x-27 Johns Hopkins University
// (C) 199x-27 University of North Carolina, Chapel Hill
// (C) 199x-27 Vanderbilt University

/// @file
/// @brief
/// @author

// Unit headers
#include <protocols/enzdes/BackboneSampler.hh>

// Package headers
#include <protocols/enzdes/BackboneSamplerCreator.hh>
#include <protocols/enzdes/EnzdesTaskOperations.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>

// Project Headers
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// Utility Headers
#include <numeric/random/random.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// C++ headers
//#include <map>
#include <set>

//basic includes
#include <basic/Tracer.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>


namespace protocols {
namespace enzdes {

static THREAD_LOCAL basic::Tracer TR( "protocols.enzdes.BackboneSampler" );

std::string
BackboneSamplerCreator::keyname() const
{
	return BackboneSamplerCreator::mover_name();
}

protocols::moves::MoverOP
BackboneSamplerCreator::create_mover() const {
	return protocols::moves::MoverOP( new BackboneSampler() );
}

std::string
BackboneSamplerCreator::mover_name()
{
	return "BackboneSampler";
}

BackboneSampler::BackboneSampler() :
	LigandBaseProtocol(),
	bb_moves_( 1000 ),
	mc_kt_( 0.6 ) {}

using namespace std;
using namespace core;
using namespace core::scoring;
using namespace pack::task;
using namespace protocols::moves;

// all constants for the backrub mover were taken from Colin's backrub.cc application
BackboneSampler::BackboneSampler
( ScoreFunctionCOP scorefxn,
	core::Size const bb_moves,
	core::Real const mc_kt
)
: protocols::ligand_docking::LigandBaseProtocol()
{
	bb_moves_ = bb_moves;
	mc_kt_ = mc_kt;
	scorefxn_repack_ = scorefxn->clone();
}

BackboneSampler::~BackboneSampler() {}

protocols::moves::MoverOP
BackboneSampler::clone() const {
	return( protocols::moves::MoverOP( new BackboneSampler( *this ) ) );
}

void
BackboneSampler::apply( Pose & pose )
{
	//  pose.dump_pdb("pose.pdb");
	kinematics::FoldTree const ft_in = pose.fold_tree();

	//checks wheather there are any segments set from cacheable observers
	//makes them obey insertions and deletions through lenght observer
	std::set< core::Size > segments;
	DetectProteinLigandInterface::add_observer_cache_segments_to_set( pose, segments );

	//makes protein segements movable
	utility::vector1< bool > mobile_bb( pose.n_residue(), false );
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	//movemap -> clear();
	if ( segments.empty() ) {
		for ( core::Size ii=1; ii <= pose.n_residue(); ++ii ) {
			if ( pose.residue(ii).is_protein() ) {
				movemap -> set_bb( ii, true );
				mobile_bb[ii] = true;
			} else {
				//C-term res in the protein is frozen
				movemap -> set_bb( ii-1, false );
				mobile_bb[ii-1] = false;
			}
		}
	} else {
		for ( core::Size ii=1; ii <= pose.n_residue(); ++ii ) {
			if ( pose.residue(ii).is_protein() && ( segments.find(ii) != segments.end() ) ) {
				movemap -> set_bb( ii, true );
				mobile_bb[ii] = true;
			}
		}
	}

	//make fold tree that complies with movemap
	//reorder_foldtree_around_mobile_regions from LigandBaseProtocol
	TR << "Original fold tree: " << pose.fold_tree() << std::endl;
	core::Size lig_jump = get_ligand_jump_id(pose);
	core::Size lig_id = get_ligand_id(pose);
	reorder_foldtree_around_mobile_regions( pose, lig_jump,  mobile_bb, lig_id);
	protocols::loops::add_cutpoint_variants( pose );
	TR << "Modified fold tree: " << pose.fold_tree() << std::endl;

	protocols::simple_moves::BBG8T3AMover bbg8t3amover;
	bbg8t3amover.movemap( movemap );
	bbg8t3amover.factorA( 0.5 ); // values suggested by Yuan
	bbg8t3amover.factorB( 10.0 );

	//core::pose::PoseCOP pose_copy = new core::pose::Pose( pose );
	//bbg8t3amover.set_input_pose( pose_copy );
	//bbg8t3amover.set_native_pose( pose_copy );

	TR << "Score After PDB Load:" << std::endl;
	scorefxn_repack_->show(TR, pose);
	TR.flush();

	protocols::moves::MonteCarlo mc( pose, *scorefxn_repack_, mc_kt_ );
	mc.reset( pose );

	TR << "Running " << bb_moves_ << " trials..." << std::endl;
	for ( core::Size ii =1; ii != bb_moves_; ++ii ) {
		bbg8t3amover.apply( pose );
		mc.boltzmann( pose, bbg8t3amover.type() );
	}
	mc.show_counters();

	TR << "Lowest score: \n";
	pose = mc.lowest_score_pose();
	scorefxn_repack_->show(TR, pose);
	TR.flush();

	protocols::loops::remove_cutpoint_variants( pose );
	pose.fold_tree( ft_in );
	TR << "Reinstated original fold tree: " << pose.fold_tree() << std::endl;
	//  pose.dump_pdb("BBG_pose.pdb");
}

std::string
BackboneSampler::get_name() const {
	return BackboneSamplerCreator::mover_name();
}

void BackboneSampler::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & )
{
	bb_moves_ = tag->getOption<core::Size>( "moves", 1000 );
	scorefxn_repack_ = protocols::rosetta_scripts::parse_score_function(tag, data)->clone();
}


} //enzdes
} //protocols
