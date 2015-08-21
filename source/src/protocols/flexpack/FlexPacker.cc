// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/flexpack/FlexPacker.cc
/// @brief
/// @author Florian Richter (floric@u.washington.edu), oct 08


// Unit Headers
#include <protocols/flexpack/FlexPacker.hh>

// Package headers
#include <protocols/flexpack/annealer/FlexbbSimAnnealer.hh>
#include <protocols/flexpack/interaction_graph/FlexbbIGFactory.hh>
#include <protocols/flexpack/interaction_graph/FlexbbInteractionGraph.hh>
#include <protocols/flexpack/rotamer_set/FlexbbRotamerSets.hh>
#include <protocols/flexpack/OtherContextScoreFunction.hh>
#include <core/fragment/Frame.hh>


// Project headers
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// tracer
#include <basic/Tracer.hh>

// ObjexxFCL
#include <ObjexxFCL/FArray1D.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace flexpack {

static thread_local basic::Tracer tr( "protocols.flexpack.FlexPacker" );

typedef core::PackerEnergy PackerEnergy;


FlexPacker::FlexPacker(
	core::pack::task::PackerTaskCOP task, // SHOULD TAKE A TASK FACTORY, NOT A TASK COP!
	utility::vector1< core::fragment::FrameOP > const & frames,
	core::scoring::ScoreFunctionCOP scorefxn
) : task_(task), scorefxn_(scorefxn)
{
	this->set_frames( frames );
}

FlexPacker::~FlexPacker() {}

void
FlexPacker::set_sfxn( core::scoring::ScoreFunctionCOP sfxn ){
	scorefxn_ = sfxn;
}

void
FlexPacker::set_task( core::pack::task::PackerTaskCOP task ){
	task_ = task;
}

void
FlexPacker::set_task_factory( core::pack::task::TaskFactoryOP factory ){
	factory_ = factory;
}


void
FlexPacker::apply(
	core::pose::Pose & pose
)
{
	using namespace ObjexxFCL;

	//std::cout << "Flexpacker start score: " << (*scorefxn_)( pose ) << std::endl;

	scorefxn_->setup_for_packing( pose, task_->repacking_residues(), task_->designing_residues() );

	rotamer_set::FlexbbRotamerSetsOP flex_rotsets( new rotamer_set::FlexbbRotamerSets( task_ ) );

	flex_rotsets->set_frames( pose, frames_ );

	/// DEBUG HACK!
	core::scoring::methods::EnergyMethodOptions opts = scorefxn_->energy_method_options();
	opts.hbond_options().decompose_bb_hb_into_pair_energies( true );
	(const_cast< core::scoring::ScoreFunction & > (*scorefxn_)).set_energy_method_options( opts );


	core::graph::GraphOP flex_neighbor_graph = flex_rotsets->flexpack_neighbor_graph( pose, *scorefxn_, task_ );

	flex_rotsets->build_rotamers( pose, *scorefxn_, *flex_neighbor_graph );
	tr << "Flexxbb RotamerSet contains a total of " << flex_rotsets->nrotamers() << " rotamers." << std::endl;

	interaction_graph::FlexbbInteractionGraphOP flex_ig =
		interaction_graph::FlexbbIGFactory::create_flexbb_interaction_graph( *task_, *flex_rotsets, pose, *scorefxn_ );

	flex_ig->initialize( *flex_rotsets );
	flex_rotsets->precompute_energies( pose, *scorefxn_, flex_neighbor_graph, *flex_ig );

	FArray1D_int bestrotamer_at_seqpos( pose.total_residue(), 0 );
	FArray1D< PackerEnergy > rot_freq( flex_rotsets->nrotamers() );
	FArray1D_int current_rot_index( pose.total_residue(), 1 );
	//for ( Size ii = 1; ii <= flex_rotsets->nmoltenres(); ++ii ) {
	// current_rot_index( flex_rotsets->moltenres_2_resid( ii ) ) = 1 + flex_rotsets->nrotamer_offset_for_moltenres( ii );
	//}
	tr << "Current rot index: ";
	for ( Size ii = 1; ii <= flex_rotsets->nmoltenres(); ++ii ) {
		tr << current_rot_index( ii ) << " ";
	}
	tr << std::endl;
	PackerEnergy bestE (0);

	annealer::FlexbbSimAnnealerOP annealer( new annealer::FlexbbSimAnnealer(
		bestrotamer_at_seqpos, bestE, false /*start_with_current*/, flex_ig,
		flex_rotsets, current_rot_index, false, rot_freq) );

	tr << "FlexbbIG Memory Use: " << flex_ig->getTotalMemoryUsage() << " bytes " << std::endl;
	OtherContextScoreFunction oc_sfxn( pose );
	using namespace core::scoring;
	oc_sfxn.set_energy_method_options( scorefxn_->energy_method_options() );
	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		ScoreType iist = (ScoreType) ii;
		if ( scorefxn_->weights()[ iist ] != 0.0 ) {
			oc_sfxn.set_weight( iist, scorefxn_->weights()[ iist ] );
		}
	}

	pose.energies().clear();
	oc_sfxn.pre_scoring();
	tr << "Other Context Score: " << oc_sfxn( pose ) << std::endl;
	pose.energies().total_energies().show_weighted( tr, oc_sfxn.weights() );
	tr << "OC weights: ";
	oc_sfxn.weights().show_nonzero( tr );
	tr << std::endl;
	pose.energies().clear();
	tr << std::endl << "Correct Context Score: " << (*scorefxn_)( pose ) << std::endl;
	pose.energies().total_energies().show_weighted( tr, scorefxn_->weights() );
	tr << std::endl;
	tr << "scorefxn_ weights: ";
	scorefxn_->weights().show_nonzero( tr );
	tr << std::endl;


	std::cerr << "Start score: " << oc_sfxn( pose ) << std::endl;
	annealer->run();

	tr << "Final energy: " << bestE << std::endl;
	for ( Size ii = 1; ii <= flex_rotsets->nmoltenres(); ++ii ) {
		Size iiresid =flex_rotsets->moltenres_2_resid( ii );
		//tr << "Replacing residue " << iiresid << " mres: " << ii << std::endl;
		pose.replace_residue(
			iiresid,
			*(flex_rotsets->rotamer( bestrotamer_at_seqpos( iiresid ))),
			false
		);
	}

	tr << "The final assigned backbone fragments are: ";
	for ( Size jj = 1; jj<= flex_rotsets->nflexible_segments(); ++jj ) {
		Size representative_seqpos = flex_rotsets->flexsegment_start_resid( jj );
		Size representative_moltenres = flex_rotsets->resid_2_moltenres( representative_seqpos );
		Size rep_rotamer = bestrotamer_at_seqpos( representative_seqpos ) - flex_rotsets->nrotamer_offset_for_moltenres( representative_moltenres );

		tr << "for segment " << jj << " is " << flex_ig->get_bb_for_state( representative_moltenres , rep_rotamer ) ;
		tr << "; ";
	}
	tr << std::endl;


	pose.energies().clear();
	tr << "Other Context Score: " << oc_sfxn( pose ) << std::endl;
	pose.energies().total_energies().show_weighted( tr, oc_sfxn.weights() );
	tr << "OC weights: ";
	oc_sfxn.weights().show_nonzero( tr );
	tr << std::endl;
	pose.energies().clear();
	tr << std::endl << "Correct Context Score: " << (*scorefxn_)( pose ) << std::endl;
	pose.energies().total_energies().show_weighted( tr, scorefxn_->weights() );
	tr << std::endl;
	tr << "scorefxn_ weights: ";
	scorefxn_->weights().show_nonzero( tr );
	tr << std::endl;

	//pose.dump_pdb( "flexpacked.pdb" );


} //apply function

std::string
FlexPacker::get_name() const {
	return "FlexPacker";
}

void
FlexPacker::set_frames(
	utility::vector1< core::fragment::FrameOP > const & frames
)
{

	for ( utility::vector1< core::fragment::FrameOP >::const_iterator frame_it = frames.begin(); frame_it != frames.end(); ++frame_it ) {

		frames_.push_back( *frame_it );
	}

} //set_frames

} //namespace flexpack
} //namespace protocols
