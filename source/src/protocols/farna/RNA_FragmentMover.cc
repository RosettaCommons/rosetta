// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to RNA_FragmentMover
/// @author Rhiju Das


// Unit headers
#include <protocols/farna/RNA_FragmentMover.hh>
#include <protocols/farna/RNA_Fragments.hh>
#include <protocols/farna/RNA_MatchType.hh>

// Project headers
#include <core/types.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <protocols/toolbox/AllowInsert.hh>
#include <protocols/toolbox/AllowInsert.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>

//C++ headers
#include <vector>
#include <string>
#include <sstream>


using namespace core;
using basic::T;


static THREAD_LOCAL basic::Tracer TR( "protocols.rna.RNA_FragmentMover" );

namespace protocols {
namespace farna {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
RNA_FragmentMover::RNA_FragmentMover(
	RNA_Fragments const & rna_fragments,
	protocols::toolbox::AllowInsertCOP allow_insert
) : Mover(),
	rna_fragments_( rna_fragments ),
	allow_insert_( allow_insert ),
	num_insertable_residues_( 0 ),
	insert_map_frag_size_( 0 ),
	frag_size_( 0 )
{
	Mover::type("RNA_FragmentMover");
}

// Copy constructor
RNA_FragmentMover::RNA_FragmentMover(RNA_FragmentMover const & object_to_copy) : Mover(object_to_copy),
	rna_fragments_(object_to_copy.rna_fragments_),
	allow_insert_(object_to_copy.allow_insert_),
	insert_map_(object_to_copy.insert_map_),
	num_insertable_residues_(object_to_copy.num_insertable_residues_),
	insert_map_frag_size_(object_to_copy.insert_map_frag_size_),
	frag_size_(object_to_copy.frag_size_)
{}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
RNA_FragmentMover::~RNA_FragmentMover()
{
}

std::string
RNA_FragmentMover::get_name() const {
	return "RNA_FragmentMover";
}

protocols::moves::MoverOP
RNA_FragmentMover::clone() const
{
	return protocols::moves::MoverOP( new RNA_FragmentMover(*this) );
}

// protocols::moves::MoverOP
// RNA_FragmentMover::fresh_instance() const
// {
// 	return protocols::moves::MoverOP( new RNA_FragmentMover() );
// }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMover::apply(
	core::pose::Pose & pose
)
{
	//Note, its better to call random_fragment_insertion directly...
	random_fragment_insertion( pose, frag_size_ );
}

///////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMover::update_insert_map( pose::Pose const & pose )
{

	// Maybe we don't need to do any updating.
	if ( frag_size_ == insert_map_frag_size_ ) return;

	// OK we do.
	num_insertable_residues_ = 0;
	insert_map_.clear();

	for ( Size i = 1; i <= pose.total_residue() - frag_size_ + 1; i++ ) {

		// Look to see if frame has *any* insertable residues.
		// This is different from the past. Now we have a atom-resolution
		// allow_insert map, so when we do the insertion, we can
		// avoid changing atom positions that should be fixed.
		//bool frame_ok = true;

		//  if ( !rna_fragments_->is_fullatom() ){
		bool frame_ok = false;
		for ( Size offset = 1; offset <= frag_size_; offset++ ) {
			Size const n = i + offset - 1;
			if ( allow_insert_->get( n ) ) { //sucka!
				frame_ok = true;
				break;
			}
		}

		//  } else {

		//   frame_ok = true;
		//   //Old school. Default for full-atom fragments.
		//   for (Size offset = 1; offset <= frag_size_; offset++ ){
		//    if ( !allow_insert_->get( i + offset - 1 ) ) { //sucka!
		//     frame_ok = false;
		//     break;
		//    }
		//   }
		//  }

		if ( !frame_ok ) continue;

		// Must make sure the whole frame is RNA, of course.

		for ( Size offset = 1; offset <= frag_size_; offset++ ) {
			if ( !pose.residue_type( i + offset - 1 ).is_RNA() ) {
				frame_ok = false; break;
			}
		}

		// Check for cutpoints that interrupt frame. Wait. why?
		//  for (Size offset = 1; offset <= frag_size_; offset++ ){
		//   if ( offset < frag_size_ &&
		//      pose.fold_tree().is_cutpoint( i + offset - 1) &&
		//      !( pose.residue_type( i+offset-1).has_variant_type( chemical::CUTPOINT_LOWER ) &&
		//       pose.residue_type( i+offset  ).has_variant_type( chemical::CUTPOINT_UPPER ) ) ) {
		//    frame_ok = false; break;
		//   }
		//  }

		if ( !frame_ok ) continue;  //  for (Size offset = 1; offset <= frag_size_; offset++ ){
		//   if ( offset < frag_size_ &&
		//      pose.fold_tree().is_cutpoint( i + offset - 1) &&
		//      !( pose.residue_type( i+offset-1).has_variant_type( chemical::CUTPOINT_LOWER ) &&
		//       pose.residue_type( i+offset  ).has_variant_type( chemical::CUTPOINT_UPPER ) ) ) {
		//    frame_ok = false; break;
		//   }


		num_insertable_residues_++;
		insert_map_[ num_insertable_residues_ ] = i;

	}

	insert_map_frag_size_ = frag_size_; //up to date!

}

////////////////////////////////////////////////////////////////////////////////////////
Size
RNA_FragmentMover::random_fragment_insertion(
	core::pose::Pose & pose,
	Size const & frag_size
)
{
	frag_size_ = frag_size;

	//Size const type = protocols::farna::MATCH_YR;
	Size const type = protocols::farna::MATCH_EXACT;

	// Make this insertion stuff a class before checking in?
	update_insert_map( pose );
	if ( num_insertable_residues_ == 0 ) return 0; // nothing to do
	Size const position_index = static_cast <int> ( numeric::random::rg().uniform() * num_insertable_residues_ ) + 1;
	Size const position = insert_map_[ position_index ];

	// std::cout << " --- Trying fragment! at " << position << std::endl;

	rna_fragments_.apply_random_fragment( pose, position, frag_size_, type, allow_insert_ );

	return position;
}


void
RNA_FragmentMover::set_frag_size(
	Size const fragment_size
)
{
	frag_size_ = fragment_size;
}

} //farna
} //protocols
