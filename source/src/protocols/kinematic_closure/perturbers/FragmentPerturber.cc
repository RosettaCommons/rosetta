// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Roland A. Pache, PhD

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/perturbers/FragmentPerturber.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/SingleResidueFragData.hh>
#include <core/fragment/BBTorsionSRFD.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/kinematic_closure/vector.hh>
#include <boost/foreach.hpp>


namespace protocols {
namespace kinematic_closure {
namespace perturbers {

using namespace std;
using numeric::conversions::DEGREES;

FragmentPerturber::FragmentPerturber(utility::vector1< core::fragment::FragSetCOP > const & frag_libs) {
	frag_libs_=frag_libs;
}

void FragmentPerturber::perturb_subset(Pose const &, IndexList const & residues, ClosureProblemOP problem) {
	
	//init variables
	Real phi, psi, omega;
	bool success (false);
	Size nfail = 0;
	const Size max_nfail = 100;//limit number of sampling attempts
	const Size num_frag_libs = frag_libs_.size();
	Size frag_lib_index;
	core::fragment::FragSetCOP frag_lib;
	core::fragment::FrameList overlapping_frames;
	const Size region_start=residues.front();
	const Size region_end=residues.back();
	basic::Tracer TR("protocols.looprelax.FragmentPerturber");
	TR << endl << "residues to sample: " << region_start << " to " << region_end << endl;
	core::kinematics::MoveMap move_map;
	//allow all motions in the movemap, since the pose itself is not changed by this method
	move_map.set_bb(true);
	move_map.set_chi(true);
	move_map.set_jump(true);
	const Size min_overlap = 1;
	const Size min_length = 1;
	Size num_overlapping_frames;
	Size overlapping_frame_index;
	core::fragment::FrameCOP overlapping_frame;
	Size num_fragments;
	Size fragment_number;
	core::fragment::FragDataCOP fragment;
	Size frag_pos;
	core::fragment::SingleResidueFragDataCOP fragment_residue;
	core::fragment::BBTorsionSRFDCOP fragment_data;
	
	//try fragment selection max_nfail times
	while (!success && (nfail < max_nfail)) {
		TR << "trial " << nfail << endl;
		
		//randomly choose a fragment library to sample from (i.e. uniform sampling of the fragment length from all available fragment lengths)
		frag_lib_index = 1 + static_cast<int>(numeric::random::uniform() * num_frag_libs);//adding 1, since the vector starts at index 1
		frag_lib = frag_libs_[frag_lib_index];
		TR << "chose random frag_lib #" << frag_lib_index << endl;
		
		//search that library for fragment alignment frames that overlap with the given residues region and contain at least one valid fragment
		if (!frag_lib->overlapping_with_region(move_map, region_start, region_end, min_overlap, min_length, overlapping_frames)) {
			nfail++;
			continue;//try again
		}
		num_overlapping_frames = overlapping_frames.size();
		/*//start debug
		 for (Size i=1; i<=num_overlapping_frames; ++i) {
		 TR << overlapping_frames[i]->type() << " " << overlapping_frames[i]->length() << " " << overlapping_frames[i]->nr_frags() << endl;
		 //overlapping_frames[i]->show(std::cout);
		 }
		 //end debug */
		TR << "number of overlapping frames found: " << num_overlapping_frames << endl;
		
		//randomly choose an alignment frame to sample from (i.e. uniform sampling of the alignment frame from all overlapping frames)
		overlapping_frame_index = 1 + static_cast<int>(numeric::random::uniform() * num_overlapping_frames);//adding 1, since the vector starts at index 1
		overlapping_frame = overlapping_frames[overlapping_frame_index];
		num_fragments = overlapping_frame->nr_frags();
		TR << "chose random overlapping frame #" << overlapping_frame_index << ", going from residue " << overlapping_frame->start() << " to " << overlapping_frame->end() << endl;
		/*//start debug
		 for (Size i=1; i<=num_fragments; ++i) {
		 TR << i << " " << (overlapping_frame->fragment_ptr(i))->sequence() << " " << (overlapping_frame->fragment_ptr(i))->secstruct() << " " << (overlapping_frame->fragment_ptr(i))->score() << endl;
		 }
		 //end debug */
		TR << "this frame is " << overlapping_frame->end()-overlapping_frame->start() << " residues long and contains " << num_fragments << " fragments" << endl;
		
		//randomly choose a fragment from the selected alignment frame (i.e. uniform fragment sampling from all fragments of the selected frame)
		fragment_number = 1 + static_cast<int>(numeric::random::uniform() * num_fragments);//adding 1, since the vector starts at index 1
		fragment = overlapping_frame->fragment_ptr(fragment_number);
		TR << "chose random fragment #" << fragment_number << " " << fragment->sequence() << " " << fragment->secstruct() << " " << fragment->score() << endl;
		success=true;
	}
	if (!success) {
		TR.Error << "couldn't find any valid fragment!" << endl;
		return;
	}
	
	BOOST_FOREACH(Size residue, residues) {
		//check if residue is part of the given frame
		if (overlapping_frame->contains_seqpos(residue)) {
			//convert sequence position to internal fragment position
			frag_pos = residue - overlapping_frame->start() + 1;//adding 1, since the vector starts at index 1
			TR << "residue: " << residue << ", frag_pos: " << frag_pos << endl;
			fragment_residue = fragment->get_residue(frag_pos);
			//dynamically cast the base class fragment data pointer to the specialized BBTorsion subtype
			//(dynamic cast required, since different fragment types can exist in parallel at runtime)
			fragment_data = core::fragment::BBTorsionSRFDCOP(dynamic_cast<core::fragment::BBTorsionSRFD const*>(fragment_residue.get()));
			if (fragment_data) {
				//if the dynamic cast was successful, extract backbone torsion angles
				phi = fragment_data->torsion(1);
				psi = fragment_data->torsion(2);
				omega = fragment_data->torsion(3);
				TR << "phi: " << phi << ", psi: " << psi << ", omega: " << omega << endl;
				//update KIC matrix with new backbone torsion angles
				problem->perturb_phi(residue, phi, DEGREES);
				problem->perturb_psi(residue, psi, DEGREES);
				problem->perturb_omega(residue, omega, DEGREES);
			}
		}
	}
}

}
}
}
