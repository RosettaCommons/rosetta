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
/// @detailed
/// @author Rhiju Das


#include <protocols/rna/RNA_FragmentMover.hh>
#include <protocols/rna/RNA_Fragments.hh>
#include <protocols/rna/AllowInsert.hh>
#include <protocols/rna/AllowInsert.fwd.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_Util.hh>
#include <core/pose/Pose.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray3D.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

// External library headers


//C++ headers
#include <vector>
#include <string>
#include <sstream>
// AUTO-REMOVED #include <fstream>
// AUTO-REMOVED #include <ctime>

#include <core/kinematics/Jump.hh>
#include <protocols/rna/RNA_MatchType.hh>
#include <utility/vector1.hh>

using namespace core;
using basic::T;

static numeric::random::RandomGenerator RG(22220);  // <- Magic number, do not change it!

static basic::Tracer TR( "protocols.rna.rna_fragment_mover" ) ;

namespace protocols {
namespace rna {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	RNA_FragmentMover::RNA_FragmentMover(
																		 RNA_FragmentsOP rna_fragments,
																		 protocols::rna::AllowInsertOP allow_insert ):
  Mover(),
	rna_fragments_( rna_fragments ),
	allow_insert_( allow_insert ),
	num_insertable_residues_( 0 ),
	insert_map_frag_size_( 0 ),
	frag_size_( 0 )
{
	Mover::type("RNA_FragmentMover");
}

	// This constructor is not actually used anymore -- better to use AllowInsert object above.
	RNA_FragmentMover::RNA_FragmentMover( RNA_FragmentsOP all_rna_fragments,
																				ObjexxFCL::FArray1D<bool> const & allow_insert_in,
																				pose::Pose const & pose ):
  Mover(),
	rna_fragments_( all_rna_fragments ),
	num_insertable_residues_( 0 ),
	insert_map_frag_size_( 0 ),
	frag_size_( 0 )
{
	Mover::type("RNA_FragmentMover");

	allow_insert_ = new AllowInsert( pose );
	allow_insert_->set( false );
	for ( Size i = 1; i <= allow_insert_in.size(); i++ ){
		if ( pose.residue_type( i ).is_RNA() && allow_insert_in[ i ] ) allow_insert_->set( i, true );
	}

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
RNA_FragmentMover::~RNA_FragmentMover()
{
}

std::string
RNA_FragmentMover::get_name() const {
	return "RNA_FragmentMover";
}

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
	if (frag_size_ == insert_map_frag_size_) return;

	// OK we do.
	num_insertable_residues_ = 0;
	insert_map_.clear();

	for (Size i = 1; i <= pose.total_residue() - frag_size_ + 1; i++ ) {

		// Look to see if frame has *any* insertable residues.
		// This is different from the past. Now we have a atom-resolution
		// allow_insert map, so when we do the insertion, we can
		// avoid changing atom positions that should be fixed.
		bool frame_ok = true;

		//		if ( !rna_fragments_->is_fullatom() ){
		frame_ok = false;
		for (Size offset = 1; offset <= frag_size_; offset++ ){
			if ( allow_insert_->get( i + offset - 1 ) ) { //sucka!
				frame_ok = true;
				break;
			}
		}

		//		} else {

		//			frame_ok = true;
		//			//Old school. Default for full-atom fragments.
		//			for (Size offset = 1; offset <= frag_size_; offset++ ){
		//				if ( !allow_insert_->get( i + offset - 1 ) ) { //sucka!
		//					frame_ok = false;
		//					break;
		//				}
		//			}
		//		}

		if ( !frame_ok ) continue;

		// Check for cutpoints that interrupt frame. Wait. why?
		//		for (Size offset = 1; offset <= frag_size_; offset++ ){
		//			if ( offset < frag_size_ &&
		//					 pose.fold_tree().is_cutpoint( i + offset - 1) &&
		//					 !( pose.residue_type( i+offset-1).has_variant_type( chemical::CUTPOINT_LOWER ) &&
		//							pose.residue_type( i+offset  ).has_variant_type( chemical::CUTPOINT_UPPER ) ) ) {
		//				frame_ok = false; break;
		//			}
		//		}

		if ( !frame_ok ) continue;

		num_insertable_residues_++;
		insert_map_[ num_insertable_residues_ ] = i;

	}

	// std::cout << "NUM_INSERTABEL_RESIDUES: " << num_insertable_residues_ << " for FRAG SIZE " << frag_size_ << std::endl;

	// std::cout << "ALLOW INSERT! ALLOW INSERT! ALLOW INSERT!" << std::endl;
	// for (Size i = 1; i <= pose.total_residue(); i++ ){
	// 	std::cout << allow_insert_->get( i );
	// }
	// std::cout << std::endl;
	// std::cout << "ALLOW INSERT! ALLOW INSERT! ALLOW INSERT!"<< std::endl;



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

	//Size const type = protocols::rna::MATCH_YR;
	Size const type = protocols::rna::MATCH_EXACT;

	// Make this insertion stuff a class before checking in?
	update_insert_map( pose );
	if ( num_insertable_residues_ == 0) return 0; // nothing to do
	Size const position_index = static_cast <int> ( RG.uniform() * num_insertable_residues_ ) + 1;
	Size const position = insert_map_[ position_index ];

	//	std::cout << " --- Trying fragment! at " << position << std::endl;

	rna_fragments_->apply_random_fragment( pose, position, frag_size_, type, allow_insert_ );

	return position;
}


void
RNA_FragmentMover::set_frag_size(
	 Size const fragment_size
)
{
	frag_size_ = fragment_size;
}

} // namespace rna
} // namespace protocols
