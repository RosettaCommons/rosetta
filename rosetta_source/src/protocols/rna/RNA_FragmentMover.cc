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
#include <protocols/rna/RNA_FragmentsClasses.hh>
// AUTO-REMOVED #include <protocols/rna/RNA_ProtocolUtil.hh>
// AUTO-REMOVED #include <protocols/rna/RNA_SecStructInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/rna/RNA_Util.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.hh>
#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <basic/options/util.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray3D.hh>

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

//Auto Headers
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <utility/options/keys/BooleanOptionKey.hh>


using namespace core;
using basic::T;

static numeric::random::RandomGenerator RG(22220);  // <- Magic number, do not change it!

static basic::Tracer TR( "protocols.rna.rna_fragment_mover" ) ;

namespace protocols {
namespace rna {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
RNA_FragmentMover::RNA_FragmentMover( RNA_FragmentsOP & all_rna_fragments, FArray1D_bool const & allow_insert ):
  Mover(),
	all_rna_fragments_( all_rna_fragments ),
	allow_insert_( allow_insert ),
	num_insertable_residues_( 0 ),
	insert_map_frag_size_( 0 ),
	frag_size_( 0 )
{
	Mover::type("RNA_FragmentMover");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
RNA_FragmentMover::~RNA_FragmentMover()
{
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

std::string
RNA_FragmentMover::get_name() const {
	return "RNA_FragmentMover";
}

void
RNA_FragmentMover::set_frag_size(
	 Size const fragment_size
)
{
	frag_size_ = fragment_size;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMover::insert_fragment(
	 core::pose::Pose & pose,
	 Size const position,
	 protocols::rna::TorsionSet const torsion_set
){

	using namespace protocols::rna;
	using namespace scoring::rna;

	Size const size = torsion_set.get_size();

	for (Size offset = 0; offset < size; offset++){

		Size const position_offset = position + offset;

		pose.set_secstruct( position_offset, torsion_set.secstruct( offset ) );

		for (Size j = 1; j <= NUM_RNA_TORSIONS; j++) {
			id::TorsionID rna_torsion_id( position_offset, id::BB, j );
			if ( j > NUM_RNA_MAINCHAIN_TORSIONS) rna_torsion_id = id::TorsionID( position_offset, id::CHI, j - NUM_RNA_MAINCHAIN_TORSIONS );
			pose.set_torsion( rna_torsion_id,
												torsion_set.torsions( j, offset ) );
		}

		//			pose.set_name( position_offset, torsion_set.torsion_source_name( offset ) );

		//TR << std::endl;
	}

	//////////////////////////////////////////////////////////////
	if ( torsion_set.non_main_chain_sugar_coords_defined ) {

		initialize_non_main_chain_sugar_atoms();

		//Force one refold.
		pose.residue(1).xyz( 1 );
		pose::Pose const & reference_pose( pose ); //This will avoid lots of refolds. I think.

		for (Size offset = 0; offset < size; offset++){

			Size const position_offset = position + offset;
			utility::vector1< Vector > vecs;
			for (Size n = 1; n <= non_main_chain_sugar_atoms.size(); n++  ) {
				Vector v( torsion_set.non_main_chain_sugar_coords( offset, n, 1) ,
									torsion_set.non_main_chain_sugar_coords( offset, n, 2) ,
									torsion_set.non_main_chain_sugar_coords( offset, n, 3) ) ;
				vecs.push_back( v );
			}
			apply_non_main_chain_sugar_coords( vecs, pose, reference_pose, position_offset );

		}
	}



}

// ///////////////////////////////////////////////////////////////////////////////
// // non_main_chain_sugar_atoms should be C2*, C1*, O4* (check out RNA_FragmentsClasses.h)
// void
// RNA_FragmentMover::apply_non_main_chain_sugar_coords_fragment(
//     FArray3D <Real > const & non_main_chain_sugar_coords,
// 		pose::Pose & pose,
// 		Size const position,
// 		Size const frag_size)
// {

// 	using namespace core::scoring::rna;
// 	using namespace core::id;
// 	initialize_non_main_chain_sugar_atoms();

// 	for (Size offset = 0; offset < frag_size; offset++){

// 		utility::vector1< Vector > vecs;

// 		for (Size n = 1; n <= non_main_chain_sugar_atoms.size(); n++  ) {
// 			Vector v( non_main_chain_sugar_coords( offset, n, 1) ,
// 								non_main_chain_sugar_coords( offset, n, 2) ,
// 								non_main_chain_sugar_coords( offset, n, 3) ) ;
// 			vecs.push_back( v );
// 		}

// 		Size const i = position+offset;

// 		utility::vector1< Real > rna_torsions;

// 		apply_non_main_chain_sugar_coords( vecs, pose, reference_pose, i );
// 	}

// }

///////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMover::update_insert_map( pose::Pose const & pose )
{
	// Maybe we don't need to do any updating.
	if (frag_size_ == insert_map_frag_size_) return;

	// OK we do.
	//	TR << "Updating insert map for frag size " << frag_size_ << std::endl;
	num_insertable_residues_ = 0;
	insert_map_.clear();

	for (Size i = 1; i <= allow_insert_.size() - frag_size_ + 1; i++ ) {

		bool frame_ok = true;
		for (Size offset = 1; offset <= frag_size_; offset++ ){
			if ( !allow_insert_( i + offset - 1 ) ) { //sucka!
				frame_ok = false; break;
			}
			if ( offset < frag_size_ && pose.fold_tree().is_cutpoint( i + offset - 1) ){
				frame_ok = false; break;
			}
		}

		if (frame_ok) {
			num_insertable_residues_++;
			insert_map_[ num_insertable_residues_ ] = i;
			//			TR << " INSERT_MAP: " << num_insertable_residues_ << " ==> " << i << std::endl;
		}

	}

	insert_map_frag_size_ = frag_size_; //up to date!

}

/////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMover::random_fragment_insertion(
	 core::pose::Pose & pose,
	 Size const & frag_size
)
{
	using namespace protocols::rna;

	frag_size_ = frag_size;

	Size const type = MATCH_EXACT;
	//Size const type = MATCH_YR;
	TorsionSet torsion_set( frag_size_ );

	// Generalize to move maps. (in fact, maybe use a move map???)
	// Size position = static_cast <int> ( RG.uniform() * ( nres - frag_size + 1) ) + 1;

	// Make this insertion stuff a class before checking in?
	update_insert_map( pose );
	if ( num_insertable_residues_ == 0) return; // nothing to do
	Size const position_index = static_cast <int> ( RG.uniform() * num_insertable_residues_ ) + 1;
	Size const position = insert_map_[ position_index ];
	all_rna_fragments_->pick_random_fragment( torsion_set, pose, position, frag_size_, type  );
	insert_fragment( pose, position, torsion_set );

}


} // namespace rna
} // namespace protocols
