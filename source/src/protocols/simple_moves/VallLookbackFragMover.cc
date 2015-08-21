// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  makes sure all fragments fall within X distance of a naturally occuring fragment during sampling.
/// @author TJ Brunette
/// Also saves the value so a score function can be used.

// Unit Headers
#include <protocols/simple_moves/VallLookbackFragMover.hh>

// history
#include <core/scoring/methods/vall_lookback/VallLookbackPotential.hh>
#include <core/scoring/methods/vall_lookback/VallLookbackData.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
// Package Headers


// Project Headers
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>

#include <core/indexed_structure_store/FragmentStore.hh>
#include <core/indexed_structure_store/FragmentLookup.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/ScoringManager.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>
#include <core/sequence/ABEGOManager.hh>

#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <protocols/moves/Mover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/indexed_structure_store.OptionKeys.gen.hh>
// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/random/random.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <map>


namespace protocols {
namespace simple_moves {


using namespace core;
using namespace core::indexed_structure_store;
static basic::Tracer tr("protocols.simple_moves.VallLookbackFragMover");

VallLookbackFragMover::VallLookbackFragMover(
	core::fragment::FragSetCOP fragset, core::kinematics::MoveMapCOP movemap) :
	ClassicFragmentMover( fragset, movemap,"VallLookbackFragMover")
{
}

VallLookbackFragMover::~VallLookbackFragMover() {}

std::string
VallLookbackFragMover::get_name() const {
	return "VallLookbackFragMover";
}
//DB stored in VallLookbackPotential. That way only 1 db is used.
Real VallLookbackFragMover::lookbackMaxRmsd(pose::Pose & pose, Size resStart,Size resEnd) const{
	using namespace core::scoring::methods;
	VallLookbackPotential potential_ = scoring::ScoringManager::get_instance()->get_vallLookbackPotential();
	return(potential_.lookbackRegion(pose, resStart, resEnd));
}

void VallLookbackFragMover::markChangedResid(pose::Pose & pose,Size resStart,Size resEnd) const{
	using namespace core::pose::datacache;
	using namespace core::scoring::methods;
	if ( !pose.data().has( CacheableDataType::VALL_LOOKBACK_DATA) ) {
		VallLookbackDataOP history = core::scoring::methods::VallLookbackDataOP(new VallLookbackData(pose));
		pose.data().set( CacheableDataType::VALL_LOOKBACK_DATA,history);
	} else {
		VallLookbackDataOP history_ = utility::pointer::static_pointer_cast<core::scoring::methods::VallLookbackData >(pose.data().get_ptr( CacheableDataType::VALL_LOOKBACK_DATA ));
		for ( Size ii=resStart; ii<=resEnd; ++ii ) {
			history_->set_res_changed(ii,true);
		}
	}
}

bool VallLookbackFragMover::apply_frames( pose::Pose &pose, FrameList const& frames ) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::indexed_structure_store;
	using namespace OptionKeys::indexed_structure_store;
	using namespace core::scoring::methods;
	using namespace core::scoring;
	using namespace core::pose::datacache;
	Size frame_num;
	Size frag_num;
	if ( !choose_fragment( frames, pose, frame_num /*output*/, frag_num /*output*/ ) ) return false;
	if ( tr.Trace.visible() ) {
		tr.Trace
			<< "frag (" << frames[ frame_num ]->start() << ","
			<< frag_num << ","
			<< frames[ frame_num ]->nr_res_affected( *movemap_ )
			<< ")" << std::endl;
	}
	//count number of positions above threshold
	VallLookbackPotential potential_ = ScoringManager::get_instance()->get_vallLookbackPotential();
	Size fragmentStore_fragment_length = potential_.fragLengthInDB();
	Size fragBeingSwapped= frames[frame_num]->start();
	//Check to see if ABEGO warrants checking all the fragments
	Size fragment_length = frames[frame_num]->length();
	Size adjRes = fragmentStore_fragment_length;
	Size resStart = 999;
	Size resEnd = 999;
	Real thresholdDistance = option[fragment_threshold_distance]();
	if ( (int)fragBeingSwapped-(int)adjRes >= 1 ) {
		resStart = fragBeingSwapped-adjRes;
	} else {
		resStart = 1;
	}
	if ( fragBeingSwapped+fragment_length+adjRes<= pose.total_residue() ) {
		resEnd = fragBeingSwapped+adjRes;
	} else {
		resEnd = pose.total_residue()-fragmentStore_fragment_length;
	}
	Real startLookbackRmsd= lookbackMaxRmsd(pose,resStart,resEnd);
	//cache fragment to be changed from initial pose so you can apply it later if the move fails.
	utility::vector1< Real >phi;
	utility::vector1< Real >psi;
	utility::vector1< Real >omega;
	utility::vector1< char > secstruct;
	utility::vector1< Real > rmsd_cache;
	VallLookbackDataOP history_ = utility::pointer::static_pointer_cast<core::scoring::methods::VallLookbackData >(pose.data().get_ptr( CacheableDataType::VALL_LOOKBACK_DATA ));
	for ( Size ii=fragBeingSwapped; ii<=fragBeingSwapped+fragment_length; ++ii ) {
		phi.push_back(pose.phi(ii));
		psi.push_back(pose.psi(ii));
		omega.push_back(pose.omega(ii));
		secstruct.push_back(pose.secstruct(ii));
	}
	for ( Size ii=resStart; ii<=resEnd; ++ii ) {
		rmsd_cache.push_back(history_->get_rmsd(ii));
	}
	bool valid = apply_fragment( *frames[ frame_num ], frag_num, *movemap_, pose );
	if ( valid ) {
		markChangedResid(pose,resStart,resEnd);
		Real endLookbackRmsd = lookbackMaxRmsd(pose,resStart,resEnd);
		//std::cout << "end:" << endLookbackRmsd <<" ,start:" << startLookbackRmsd << " ,thresh" << thresholdDistance << std::endl;
		if ( endLookbackRmsd <= startLookbackRmsd || endLookbackRmsd <= thresholdDistance ) {
			return true;
		} else { //report false & go back to cached pose
			for ( Size ii=1; ii<=fragment_length; ++ii ) {
				pose.set_phi(fragBeingSwapped+ii-1,phi.at(ii));
				pose.set_psi(fragBeingSwapped+ii-1,psi.at(ii));
				pose.set_omega(fragBeingSwapped+ii-1,omega.at(ii));
				pose.set_secstruct(fragBeingSwapped+ii-1,secstruct.at(ii));
				history_->set_rmsd(fragBeingSwapped+ii-1,rmsd_cache.at(ii));
			}
			for ( Size ii=resStart; ii<=resEnd; ++ii ) {
				history_->set_rmsd(ii,rmsd_cache.at(ii-resStart+1));
			}
			return false;
		}
	} else {
		return false; //return from valid check up above.
	}
}
} // simple_moves
} // protocols
