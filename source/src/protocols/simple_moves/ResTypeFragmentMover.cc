// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @brief  swaps both fragment xyz and amino acid
/// @author TJ Brunette (tjbrunette@gmail.com)

// Unit Headers
#include <protocols/simple_moves/ResTypeFragmentMover.hh>

// Package Headers

// Project Headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/types.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <protocols/moves/Mover.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/random/random.hh>

//// C++ headers
#include <cstdlib>
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

using namespace core;
using namespace fragment;
using namespace basic;

static thread_local basic::Tracer tr( "protocols.simple_moves.ResTypeFragmentMover" );

ResTypeFragmentMover::ResTypeFragmentMover(core::fragment::FragSetCOP fragset):ClassicFragmentMover(fragset, "ResTypeFragmentMover"){}

ResTypeFragmentMover::ResTypeFragmentMover(core::fragment::FragSetCOP fragset,	core::kinematics::MoveMapCOP movemap):ClassicFragmentMover( fragset, movemap, "ResTypeFragmentMover"){}


ResTypeFragmentMover::~ResTypeFragmentMover() {}

bool ResTypeFragmentMover::apply_frames( pose::Pose &pose, core::fragment::FrameList const& frames ) const {
	Size frame_num;
	Size frag_num;
	bool success( false );
	if ( !choose_fragment( frames, pose, frame_num /*output*/, frag_num /*output*/ ) ) return false;

	if ( tr.Trace.visible() )
		tr.Trace<< "frag (" << frames[ frame_num ]->start() << ","
															<< frag_num << ","
															<< frames[ frame_num ]->nr_res_affected( *movemap_ )
															<< ")" << std::endl;
	if (tr.Trace.visible() )
		tr.Trace <<"seq " << frames[frame_num]->fragment(frag_num).sequence() << std::endl;

	if ( !check_ss() ){
		std::string swap_sequence = frames[frame_num]->fragment(frag_num).sequence();
		swap_residue_types(pose,swap_sequence,frames[ frame_num ]->start());
		return apply_fragment( *frames[ frame_num ], frag_num, *movemap_, pose );
	}
	// now do the ss-check!
	//	tr.Trace << "now do the ss-check!"<< std::endl;
	// get actual ss from pose
	std::string proposed_ss;
	proposed_ss.reserve( pose.total_residue() );
	proposed_ss = pose.secstruct();

	std::string old_ss = proposed_ss;

	// check if old ss is valid
	bool valid = !valid_ss( old_ss ); // if old_ss is not valid we can apply the fragment anyway

	// if old ss was fine ---> check fragments effect on ss
	if ( !valid ) { // if old_ss was valid we check if proposed_ss is still valid.
		frames[ frame_num ]->apply_ss( *movemap_, frag_num, proposed_ss );
		//		tr.Trace << !valid << " old_ss: " << old_ss << std::endl;
		valid = valid_ss( proposed_ss );
		//		tr.Trace << valid << "new_ss: " << proposed_ss << std::endl;
	}
	//	tr.Trace << "finished the ss-check! : " << valid << std::endl;
	if ( valid ) {
		std::string swap_sequence = frames[frame_num]->fragment(frag_num).sequence();
		swap_residue_types(pose,swap_sequence,frames[ frame_num ]->start());
		success = apply_fragment( *frames[ frame_num ], frag_num, *movemap_, pose );
	} else {
		//		tr.Trace << "dissallow insertion due to short helix/strand " << std::endl;
	}
	return success;
}

void ResTypeFragmentMover::swap_residue_types( pose::Pose &pose, std::string const sequence, Size const startSeqPos ) const {
	using namespace chemical;
	for(Size ii=0; ii<sequence.size(); ++ii){
		char aa = sequence[ii];
		AA my_aa = aa_from_oneletter_code( aa );
		ResidueTypeSetCAP const &residue_set(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID ));
		ResidueTypeCOPs const & rsd_type_list( residue_set->aa_map( my_aa ) );
		Size best_index = 1;
		ResidueType const & rsd_type( *(rsd_type_list[ best_index ]) );
		replace_pose_residue_copying_existing_coordinates(pose,startSeqPos+ii,rsd_type);
	}
}

} // simple_moves
} // protocols
