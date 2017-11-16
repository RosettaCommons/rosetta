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
/// @author Yifan Song
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_loop_grower_SheetSampler_hh
#define INCLUDED_protocols_loop_grower_SheetSampler_hh

#include <protocols/loop_grower/SheetSampler.fwd.hh>

#include <iostream>
#include <fstream>

#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FrameList.hh>
#include <map>

#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/FoldTree.hh>

#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>

#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>

//possibily duplicate includes here

#include <basic/database/open.hh>

#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/JumpingFrame.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <protocols/jumping/PairingLibrary.hh>



namespace protocols {
namespace loop_grower {

static basic::Tracer TR("protocols.loop_grower.SheetSampler");

class SheetSampler: public protocols::moves::Mover {
	protocols::jumping::PairingLibraryOP instance_;
public:
	SheetSampler():
		clashtolerance_(5),ideal_sheets_(true) {}

	SheetSampler(core::Size start, core::Size end,  core::Size frag){
		instance_ = protocols::jumping::PairingLibraryOP ( new protocols::jumping::PairingLibrary );
		instance_->read_from_file( basic::database::full_name("scoring/score_functions/jump_templates_SSpairs_v2.dat") );
		start_ = start;
		end_ = end;
		frag_ = frag;
		ideal_sheets_ = true;
		runtime_assert( end_ >= start_  );

		//set default score function settings use set_sf to set custom scorefunction
		sf_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
		sf_->set_weight( core::scoring::hbond_sr_bb, 0);
		sf_->set_weight( core::scoring::hbond_lr_bb, 0);
		sf_->set_weight( core::scoring::vdw, 1);
		sf_->set_weight( core::scoring::elec_dens_fast, 0);
		clashtolerance_ = 5;
	}

	// orientation: 1 = parallel, 2 = antiparallel
	void
	generate_jump_frags( core::pose::Pose const & pose, core::Size jump_nr, int orientation, int pleating, core::fragment::FrameList& all_frames ) {

		int const startpos( pose.fold_tree().downstream_jump_residue( jump_nr ) );
		int const endpos( pose.fold_tree().upstream_jump_residue( jump_nr ) );

		// create storage for jump fragments
		core::fragment::FragDataOPs frag_data;
		core::Size const length( 2 ); // incl endpoint torsions
		core::fragment::JumpingFrameOP frame( new core::fragment::JumpingFrame( startpos, endpos, length ) );
		frame->set_pos( 1, startpos ); frame->set_pos( 2, endpos );
		//frame->set_pos( 3, endpos );   frame->set_pos( 4, endpos );

		instance_->create_jump_fragments( orientation, pleating, false, frag_data );
		frame->add_fragment( frag_data );
		all_frames.push_back( frame );
	}

	void
	alignPerfectCA( core::pose::Pose & pose, core::Size moving, core::Size ref );

	void
	alignStrand( core::pose::Pose & pose, core::Size ref, core::Size moving, core::Size strand_start, core::Size strand_size );

	core::Real
	sheethbonds( core::pose::Pose & pose, core::Size lower, core::Size upper );

	void apply(core::pose::Pose & pose);

	void
	set_sf( core::scoring::ScoreFunctionOP newsf ){ sf_ = newsf; }

	void set_ideal_sheets( bool ideal_sheets ){ ideal_sheets_ = ideal_sheets; }

	std::string get_name() const { return "SheetSampler"; }



private:
	core::Size start_, end_, frag_;
	core::scoring::ScoreFunctionOP sf_;
	core::Real clashtolerance_;
	bool ideal_sheets_;
	//protocols::jumping::PairingLibraryOP instance_ = protocols::jumping::PairingLibraryOP ( new protocols::jumping::PairingLibrary );

};



} //loop_grower
} //protocols

#endif
