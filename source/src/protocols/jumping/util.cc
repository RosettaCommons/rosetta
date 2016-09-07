// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
/// @author Oliver Lange


// Unit Headers
#include <protocols/jumping/util.hh>
#include <core/scoring/dssp/Dssp.hh>

// Package Headers
// Project Headers
#include <core/pose/Pose.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>
#include <core/fragment/JumpingFrame.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/JumpSRFD.hh>
#include <core/conformation/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymDof.hh>

#include <protocols/loops/loop_closure/ccd/SlidingWindowLoopClosure.hh>
#include <protocols/loops/Exceptions.hh>

#include <protocols/checkpoint/CheckPointer.hh>

#ifdef WIN32
#include <core/scoring/dssp/PairingsList.hh>
#endif


// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

//numeric headers

//// C++ headers
// #include <string>
#include <list>

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>


namespace protocols {
namespace jumping {
using namespace core;
using namespace pose;
using namespace kinematics;


static THREAD_LOCAL basic::Tracer tr( "protocols.jumping" );


// if false pose does still contain jumps in the fold-tree
void
close_chainbreaks(
	loops::loop_closure::ccd::SlidingWindowLoopClosureOP closure_protocol,
	Pose& open_pose,
	protocols::checkpoint::CheckPointer &checkpoint,
	const std::string &decoy_tag,
	core::kinematics::FoldTree const& final_fold_tree_in  /*default empty fold-tree */
)
{
	core::kinematics::FoldTree final_fold_tree;
	kinematics::MoveMap const& movemap = closure_protocol->movemap();
	if ( !final_fold_tree_in.nres() ) { //empty fold-tree
		final_fold_tree.simple_tree( open_pose.total_residue() );
	} else if ( core::pose::symmetry::is_symmetric( open_pose ) ) {
		final_fold_tree = core::conformation::symmetry::replaced_symmetric_foldtree_with_new_monomer(
			open_pose.fold_tree(), *core::pose::symmetry::symmetry_info( open_pose ),
			final_fold_tree_in );
	} else {
		final_fold_tree = final_fold_tree_in;
	}
	//some sanity checks
	runtime_assert( closure_protocol != nullptr );
	//runtime_assert  that all cut-points in final_fold_tree are also contained the actual fold_tree
	for ( Size ncut = 1; ncut <= (Size)final_fold_tree.num_cutpoint(); ncut++ ) {
		if ( !open_pose.fold_tree().is_cutpoint( final_fold_tree.cutpoint( ncut ) ) ) {
			throw( loops::EXCN_Loop_not_closed( "Foldtree mismatch." ) );
		}
	}
	tr.Debug << "close chainbreaks until final fold-tree is reached\n";
	simple_visualize_fold_tree( final_fold_tree, tr.Debug );
	//close cuts that have a lot of freedom to move first
	//since removal of some jumps might make other cuts removable
	//hence we first get a list sorted by their max-loop-length
	//after each cut removal we recompute this list.
	//do this until all cuts that are not in final_fold_tree are removed
	Pose pose = open_pose;

	Size close_count = 0;
	while ( pose.fold_tree().num_cutpoint() > final_fold_tree.num_cutpoint() ) {
		close_count++;
		//make list of all cutpoints that need removal, list of std::pairs: first: max-loop-size; second: cutpoint
		std::list< std::pair< Size, Size > > cuts;
		for ( Size ncut = 1; ncut <= (Size) pose.fold_tree().num_cutpoint(); ncut++ ) {
			Size cutpoint = pose.fold_tree().cutpoint( ncut );
			if ( !final_fold_tree.is_cutpoint( cutpoint  ) ) {

				// compute  max_loop_size...
				// extend loop from cutpoint away until either a jump-residue or an unmovable bb-torsion is found
				Size min_loop_begin ( cutpoint + 1 );
				Size max_loop_end  ( cutpoint );
				while (
						min_loop_begin > 1
						&& !pose.fold_tree().is_jump_point( min_loop_begin - 1 )
						&& movemap.get_bb( min_loop_begin - 1 )
						) --min_loop_begin;
				while (
						max_loop_end < pose.total_residue()
						&& !pose.fold_tree().is_jump_point( max_loop_end + 1 )
						&& movemap.get_bb( max_loop_end + 1)
						) ++max_loop_end;
				if ( max_loop_end > min_loop_begin ) {    // put in list
					cuts.push_back( std::make_pair(  max_loop_end - min_loop_begin, cutpoint ) );
				}
			}
		}
		cuts.sort();

		if ( tr.Debug.visible() ) {
			for (auto & cut : cuts) {
				tr.Debug << "cut " << cut.second << " max_loop_size " << cut.first << std::endl;
			}
		}
		if ( cuts.empty() ) { //size() == 0 ) {
			throw( loops::EXCN_Loop_not_closed( "no moveable piece to close loop" ) );
		}

		Size const cutpoint( cuts.back().second ); //largest is last... so take last element
		tr.Info << "close chainbreak at position " << cutpoint << "..." << std::endl;

		// TODO need to also save foldtree here ! Use binary silent file format!
		const bool checkpoint_foldtree=true;


		if ( ! checkpoint.recover_checkpoint(
				open_pose, nullptr, decoy_tag,  "_lc_" + ObjexxFCL::string_of( close_count ), false, checkpoint_foldtree ) ) {
			if ( remove_cut( cutpoint, pose, final_fold_tree ) ) {
				if ( tr.Debug.visible() ) {
					tr.Debug << "current fold-trees " << std::endl;
					simple_visualize_fold_tree_and_movemap( open_pose.fold_tree(), movemap, tr.Debug );
					simple_visualize_fold_tree_and_movemap( pose.fold_tree(), movemap, tr.Debug );
					tr.Debug << std::endl;
				}
				closure_protocol->set_loop( closure_protocol->determine_loop( open_pose, pose ) );
				closure_protocol->apply( open_pose, pose );
				if ( closure_protocol->bIdealLoopClosing() ) {
					open_pose = pose;
				} else {
					//the best closing fragment is already applied to open_pose, now fix the fold-tree leaving an unideal residue at the chainbreak.
					open_pose.fold_tree( pose.fold_tree() );
				}
			}
			checkpoint.checkpoint( open_pose, nullptr, decoy_tag,  "_lc_" + ObjexxFCL::string_of( close_count ), checkpoint_foldtree );
		} else {
			// make sure the two poses have the same foldtree brefore removing one of the cuts in
			// the closed pose. (this is not necessary when checkpointing is not active, but since only the
			// open_pose is recovered, this is necessary here.
			pose = open_pose;
		}
	}
}

bool
remove_cut( Size cutpoint, Pose& pose, FoldTree const& final_fold_tree /*default empty fold-tree */) {
	FoldTree new_fold_tree = pose.fold_tree();
	bool success( remove_cut( cutpoint, new_fold_tree, final_fold_tree ) );
	if ( success ) {
		//  pose.dump_pdb("tt_old_f.pdb");
		tr.Debug << "set new fold-tree " << new_fold_tree << std::endl;
		pose.fold_tree( new_fold_tree );
		tr.Debug << "idealize positions " << std::endl;
		conformation::idealize_position( cutpoint, pose.conformation() );
		conformation::idealize_position( cutpoint+1, pose.conformation() );
	}
	return success;
}

bool
remove_cut( Size cutpoint, FoldTree &new_fold_tree, FoldTree const& final_fold_tree /*default empty fold-tree */ )
{
	tr.Info << "close-loops: remove cuts until fold-tree is : " << final_fold_tree << std::endl;
	// construct the new tree formed when we glue this cutpoint and
	// delete a single jump
	FoldTree old_fold_tree( new_fold_tree );
	FoldTree const& f ( old_fold_tree );
	runtime_assert( f.is_cutpoint( cutpoint ) );

	//find enclosing jump points
	Size jump_pos1( cutpoint );
	Size jump_pos2( cutpoint + 1 );
	while ( jump_pos1 > 1 && !f.is_jump_point( jump_pos1 ) )
			--jump_pos1;
	while ( jump_pos2 < f.nres() && !f.is_jump_point( jump_pos2 ) )
			++jump_pos2;
	runtime_assert( jump_pos1 <= jump_pos2 );

	tr.Info << "remove cutpoint: " << cutpoint << " between " << jump_pos1 << " " << jump_pos2 << " in " << std::endl;
	new_fold_tree.reorder( 1 );
	tr.Info << new_fold_tree << std::endl;


	if ( jump_pos1 < cutpoint ) {
		new_fold_tree.delete_unordered_edge( jump_pos1, cutpoint, -1);
	}
	if ( cutpoint+1 < jump_pos2 ) {
		new_fold_tree.delete_unordered_edge( cutpoint+1, jump_pos2, -1);
	}
	new_fold_tree.add_edge( jump_pos1, jump_pos2, -1 );
	// I think there may be more than one jump which
	// we could delete. This just chooses the first one
	for ( Size i=1; i<= f.num_jump(); ++i ) {
		// mark as "deleted" for the purposes of connectivity checking
		bool in_final( false );
		for ( Size nf = 1; !in_final && nf<= final_fold_tree.num_jump(); nf ++ ) {
			in_final =  ( final_fold_tree.jump_point(1, nf) == f.jump_point(1,i) && final_fold_tree.jump_point(2, nf) == f.jump_point(2,i) )
				|| ( final_fold_tree.jump_point(1, nf) == f.jump_point(2,i) && final_fold_tree.jump_point(2, nf) == f.jump_point(1,i) );
		}
		if ( !in_final ) {
			new_fold_tree.update_edge_label( f.jump_point(1,i), f.jump_point(2,i), i,  0 ); // label with 0
			if ( new_fold_tree.connected() ) {
				// safe to cut this edge
				new_fold_tree.delete_unordered_edge( f.jump_point(1,i), f.jump_point(2,i), 0 );
			} else {
				// dont cut this edge after all
				new_fold_tree.update_edge_label( f.jump_point(1,i), f.jump_point(2,i), 0, i );
			}
		}
	}
	tr.Debug << "new_fold_tree: " << new_fold_tree << std::endl;
	if ( new_fold_tree.connected() ) {
		// we deleted a jump so the jump label numbers may be screwed up
		// note that this will screw up the jump_transform mapping
		new_fold_tree.renumber_jumps();
		tr.Debug << "new_fold_tree renumbered: " << new_fold_tree << std::endl;
		// plus, there may be vertices in the tree which are not jumps and
		// not at a break
		new_fold_tree.delete_extra_vertices();
		tr.Debug << "new_fold_tree deleted_extra_vertices: " << new_fold_tree << std::endl;
		new_fold_tree.reorder( 1 );
		tr.Debug << "new_fold_tree reordered to 1 " << new_fold_tree << std::endl;
		// fold_tree = new_fold_tree;
		Size new_root = 1;
		if ( new_fold_tree.num_jump() > 0 ) {
			new_root = new_fold_tree.jump_point( 1, 1 );
		}
		new_fold_tree.reorder( new_root );
		tr.Debug << "new_fold_tree reordered to " << new_root << " " << new_fold_tree << std::endl;
		return true;
	}
	return false;
}


void safe_secstruct( pose::Pose& pose ) {
	kinematics::FoldTree const& fold_tree( pose.fold_tree() );

	Size const num_jump ( fold_tree.num_jump() );
	Size const nres( pose.total_residue() );

	for ( Size i = 1; i <= num_jump; ++i ) {
		for ( Size j = 1; j <= 2; ++j ) {
			Size const pos ( j==1 ? fold_tree.jump_edge( i ).start() : fold_tree.jump_edge( i ).stop() );
			char const ss( pose.secstruct( pos ) );
			if ( ss != 'L' ) {
				for ( Size k = std::max( (Size) 1, pos-2 ), ke = std::min( nres, pos+2 );
						k <= ke; ++k ) {
					if ( pose.secstruct(k) != ss ) {
						pose.set_secstruct( k, ss );
					}
				}
			}
		}
	}
}

core::fragment::JumpingFrameOP generate_empty_jump_frame( Size startpos, Size endpos, Size length ) {
	using namespace core::fragment;
	JumpingFrameOP frame( new JumpingFrame( startpos, endpos, length ) );
	if ( length <= 1 || length > 4 ) utility_exit_with_message("called generate_jump_frame with inappropriate length argument");
	Size pos = 1;
	if ( length == 4 ) frame->set_pos( pos++, startpos );
	frame->set_pos( pos++, startpos );
	frame->set_pos( pos++, endpos );
	if ( length >=3 ) frame->set_pos( pos, endpos );
	return frame;
}

core::fragment::JumpingFrameOP generate_jump_frame( Size startpos, Size endpos, bool bWithTorsion ) {
	using namespace core::fragment;
	FragDataOP frag_data( new FragData );
	if ( bWithTorsion ) {
		BBTorsionSRFDOP start( new BBTorsionSRFD( 3, 'E', 'X' ) );
		frag_data->add_residue( start );
	}
	frag_data->add_residue( SingleResidueFragDataOP( new UpJumpSRFD ) );
	frag_data->add_residue( SingleResidueFragDataOP( new DownJumpSRFD ) );

	if ( bWithTorsion ) {
		BBTorsionSRFDOP stop( new BBTorsionSRFD( 3, 'E', 'X' ) );
		frag_data->add_residue( stop );
	}
	JumpingFrameOP frame = generate_empty_jump_frame( startpos, endpos, frag_data->size() ); //see above
	frame->add_fragment( frag_data );
	return frame;
}

///////////////////////////////////////////////////////////////////////////////
//// some legacy code
///
typedef utility::vector1< PointPosition > PointList;

void
get_CA_vectors(
	PointList const & ca1, // pass by reference, so no tricks:: 3x3
	PointList const & ca2, // pass by reference, so no tricks:: 3x3
	Vector & a,
	Vector & b,
	Vector & c
)
{
	/*       a goes from c-alpha #1 to c-alpha #3 */
	a = ca1[3]-ca1[1];
	a.normalize();

	/*       b gives direction of pleat for ca1 c-alphas */
	Vector b1 = ca1[2] - ca1[1];
	Vector b2 = ca1[2] - ca1[3];
	b = b1 + b2 ;
	b.normalize();

	/*       c goes from ca1 triple to ca2 triple (central alpha-carbons) */
	c = ca2[2] - ca1[2];
	c.normalize();
}

void
get_pairing_geometry(
	pose::Pose const& pose,
	Size const res1,
	Size const res2,
	Real& orientation,
	Real& pleating1,
	Real& pleating2
)
{
	runtime_assert( res1>1 && res1 < pose.total_residue() &&
		res2 > res1 && res2 < pose.total_residue() );

	chemical::ResidueType const& rt1 ( pose.residue_type ( res1 ) );
	chemical::ResidueType const& rt2 ( pose.residue_type ( res2 ) );

	PointList pCA1(3); //space for 3 CAs positions
	PointList pCA2(3);

	// read CAs of 3 consecutive residues
	for ( int i = -1 ; i<=1 ; i++ ) {
		id::AtomID CA1( rt1.atom_index ("CA") , res1+i );
		id::AtomID CA2( rt2.atom_index ("CA") , res2+i );
		pCA1[ i+2 ] = pose.xyz( CA1 );
		pCA2[ i+2 ] = pose.xyz( CA2 );
	};

	Vector dCaCa = pCA1[ 2 ] - pCA2[ 2 ];
	if ( dCaCa.length() > 10.5 ) {
		tr.Warning << "the CA-CA distance for pairing " << res1 << " " << res2 << " is " << dCaCa.length() << std::endl;
		// I put this in because I had this case once, and such an exit would have saved me time.
		//If you have longer distances you probably don't want to choose pleating and orientation based on the native structure
		// --- which is a temporary convenience hack anyway.
		//  utility_exit_with_message(" the CA-CA distance for the chosen pairing is more than 10.5 A check your pairings!!! ");
	}


	Vector a1,b1,c1;
	Vector a2,b2,c2;
	get_CA_vectors( pCA1, pCA2, a1,b1,c1 );
	get_CA_vectors( pCA2, pCA1, a2,b2,c2 );

	orientation = dot_product( a1, a2 );

	Vector ab1,ab2;
	cross(a1,b1,ab1);
	cross(a2,b2,ab2);

	Real const d1 = dot_product(ab1,c1);
	Real const d2 = dot_product(ab2,c2);

	pleating1 = d1;

	if ( orientation < 0 ) {
		pleating2 =  d2; // antiparallel
	} else {
		pleating2 = -d2;
	}
}

void
get_pleating(
	pose::Pose const& pose,
	Size const pos1,
	Size const pos2,
	Size &orientation,
	Size &pleating
)
{

	//Why did this have to get so complicated?
	// Its actually a pretty simple concept!
	//
	// For some reasons, get_pairing_geometry flips
	// pleating2 depending on the orientation --
	// in ideal strand pairs, pleating1 and pleating2 then have the same sign.
	//
	// But for some twisted strand pairs (see, e.g, 22,48 in 1brs.pdb),
	// the numbers get a little crazy...


	Real forientation, pleating1, pleating2;
	if ( pos1 < pos2 ) {
		get_pairing_geometry( pose, pos1, pos2, forientation,pleating1,pleating2);

		//This isn't always quite true...
		//  runtime_assert( pleating1 * pleating2 > 0.0 );
		pleating = ( (pleating1+pleating2) < 0 ? 1 : 2 );
	} else {
		get_pairing_geometry( pose, pos2, pos1, forientation,pleating1,pleating2);


		//This isn't always quite true...
		//  runtime_assert( pleating1 * pleating2 > 0.0 );

		if ( forientation < 0 ) {
			// pleating for anti-parallel pairings is preserved when we
			// interchange positions
			pleating = ( (pleating1+pleating2) < 0 ? 1 : 2 );
		} else {
			// pleating for parallel pairings is reversed when we
			// interchange positions
			pleating = ( (pleating1+pleating2) < 0 ? 2 : 1 );
		}
	}
	tr.Debug << " orientation " << forientation << " pleating " << pleating << std::endl;
	orientation = forientation < 0 ? 1 : 2;
}

void assign_ss_dssp( core::pose::Pose & pose ) {
	core::scoring::dssp::Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );
}

} // jumping
} // protocols

