// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  some utils for fragments
/// @author Oliver Lange

// Unit Headers
#include <core/fragment/util.hh>

// Package Headers
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/JumpSRFD.hh>

#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>

// Project Headers
#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh> //for reading in fragsets from a pdb

#include <core/pose/annotated_sequence.hh>
#include <core/conformation/util.hh>

#include <core/types.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>

//#include <utility/pointer/ReferenceCount.hh>
//#include <numeric/numeric.functions.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/xyz.functions.hh>

//#include <core/pose/util.hh>
#include <utility/io/ozstream.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/scoring/ScoreFunction.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <fstream>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>

namespace core {
namespace fragment {

static THREAD_LOCAL basic::Tracer tr( "core.fragment" );

using namespace ObjexxFCL::format;

void retain_top(core::Size k, FragSetOP fragments) {
	for ( FrameIterator i = fragments->nonconst_begin(); i != fragments->nonconst_end(); ++i ) {
		FrameOP existing_frame = *i;

		// create a new frame containing only the desired fragments
		Frame new_frame(existing_frame->start(),
			existing_frame->stop(),
			existing_frame->length());

		for ( core::Size j = 1; j <= std::min(k, existing_frame->nr_frags()); ++j ) {
			new_frame.add_fragment(existing_frame->fragment_ptr(j));
		}

		*existing_frame = new_frame;
	}
}

void steal_constant_length_frag_set_from_pose ( pose::Pose const& pose_in, ConstantLengthFragSet& fragset ) {
	//Size nbb ( 3 ); // three backbone torsions for Protein
	Size len = fragset.max_frag_length();
	runtime_assert( len > 0 );
	FrameOP frame;
	pose::Pose pose = pose_in;
	pose::set_ss_from_phipsi( pose );
	for ( Size pos = 1; pos <= pose.size() - len + 1; ++pos ) {

		//don't want to use fragments that go across a cutpoint... certainly non-ideal geometry there
		// check for cutpoints
		bool free_of_cut = true;
		for ( Size icut = 1; icut <= (Size) pose.fold_tree().num_cutpoint(); icut++ ) {
			Size const cut ( pose.fold_tree().cutpoint( icut ) );
			if ( cut >= pos && cut < pos+len-1 ) {
				free_of_cut = false;
				break;
			}
		}

		//steal backbone torsion fragment
		if ( free_of_cut ) {
			frame = FrameOP( new Frame( pos, FragDataCOP( FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), len ) ) ) ) );
			frame->steal( pose );
			fragset.add( frame );
		}
	}
}

void steal_frag_set_from_pose ( pose::Pose const& pose_in, FragSet& fragset, core::fragment::FragDataCOP frag_type ) {
	//Size nbb ( 3 ); // three backbone torsions for Protein
	Size const len( frag_type->size() );
	runtime_assert( len > 0 );
	FrameOP frame;
	pose::Pose pose = pose_in;
	pose::set_ss_from_phipsi( pose );
	for ( Size pos = 1; pos <= pose.size() - len + 1; ++pos ) {
		frame = FrameOP( new Frame( pos, frag_type ) );
		frame->steal( pose );
		fragset.add( frame );
	}
}

void steal_frag_set_from_pose ( pose::Pose const& pose_in, Size const begin, Size const end, FragSet& fragset, core::fragment::FragDataCOP frag_type ) {

	Size const len( frag_type->size() );
	runtime_assert( len > 0 );
	FrameOP frame;
	pose::Pose pose = pose_in;
	pose::set_ss_from_phipsi( pose );
	for ( Size pos = begin; pos <= end - len + 1; ++pos ) {
		frame = FrameOP( new Frame( pos, frag_type ) );
		frame->steal( pose );
		fragset.add( frame );
	}
}

void steal_frag_set_from_pose (
	pose::Pose const& pose_in,
	FragSet& fragset,
	core::fragment::FragDataCOP frag_type,
	std::set< core::Size > const& selected_residues )
{
	//Size nbb ( 3 ); // three backbone torsions for Protein
	Size const len( frag_type->size() );
	runtime_assert( len == 1 ); // for this type of residue-by-residue we only support 1-mers
	FrameOP frame;
	pose::Pose pose = pose_in;
	pose::set_ss_from_phipsi( pose );
	for ( core::Size selected_residue : selected_residues ) {
		frame = FrameOP( new Frame( selected_residue, frag_type ) );
		frame->steal( pose );
		fragset.add( frame );
	}
}

// chop fragments into sub-fragments
void chop_fragments( core::fragment::FragSet& source, core::fragment::FragSet& dest ) {
	Size slen( source.max_frag_length() );
	Size tlen( dest.max_frag_length() );
	runtime_assert( tlen < slen );
	FrameList dest_frames;
	for ( Size pos = 1; pos <= source.max_pos() + slen - tlen; pos++ ) {
		dest_frames.push_back( core::fragment::FrameOP( new Frame( pos, tlen ) ) );
	}
	for ( ConstFrameIterator it=source.begin(), eit=source.end(); it!=eit; ++it ) {
		Frame const& fr( **it );
		for ( Size pos = fr.start(); pos<= fr.end() - tlen + 1; pos++ ) {
			Frame& dest_fr( *dest_frames[ pos ] );
			for ( Size nr = 1; nr <= fr.nr_frags(); ++nr ) {
				FragData const& long_frag( fr.fragment( nr ) );
				dest_fr.add_fragment( long_frag.generate_sub_fragment( pos-fr.start()+1, pos-fr.start()+tlen ) );
			}
		}
	}
	for ( FrameList::const_iterator it = dest_frames.begin(), eit = dest_frames.end(); it!=eit; ++it ) {
		if ( (*it)->nr_frags() ) {
			dest.add( *it );
		}
	}
}

void compute_per_residue_coverage( core::fragment::FragSet const& _frags, utility::vector1< core::Size > &nr_frags ) {
	if ( nr_frags.size() < _frags.max_pos() ) {
		nr_frags.resize( _frags.max_pos(), 0 );
	}

	for ( ConstFrameIterator it=_frags.begin(), eit=_frags.end(); it!=eit; ++it ) {
		// so far this implementation doesn't work for non-continous frames
		runtime_assert( it->is_continuous() );

		for ( Size i = it->start(); i<=it->stop(); i++ ) {
			nr_frags[ i ] += it->nr_frags();
		}
	}
}

void flatten_list( FrameList& frames, FragID_List& frag_ids ) {
	frag_ids.reserve( frag_ids.size() + frames.flat_size() );
	for ( FragID_Iterator it = frames.begin(), eit=frames.end();
			it!=eit; ++it ) {
		frag_ids.push_back( *it );
	}
}

FragSetOP
merge_frags( FragSet const& good_frags, FragSet const& filling, Size min_nr_frags, bool bRandom ) {

	// all good_frags go into final set
	FragSetOP merged_frags = good_frags.clone();

	utility::vector1< Size > nr_frags( filling.max_pos(), 0 );
	compute_per_residue_coverage( good_frags, nr_frags );
	tr.Info << "fragment coverage:";
	for ( Size pos = 1; pos <= nr_frags.size(); pos++ ) {
		tr.Info << "    " << pos << " " << nr_frags[ pos ];
	}
	tr.Info << std::endl;
	//fill up with filling fragments:
	for ( Size pos = 1; pos<=filling.max_pos(); pos++ ) {
		if ( nr_frags[ pos ] < min_nr_frags ) {
			FrameList fill_frames;
			filling.frames( pos, fill_frames );
			Size nr_fill( min_nr_frags - nr_frags[ pos ] );
			tr.Info << nr_frags[ pos ] << " fragments at pos " << pos << ". required: " << min_nr_frags << std::endl;
			tr.Info << "attempt to fill up with " << nr_fill << " frags at position " << pos << " ... ";
			// select randomly from filling?
			// generate random sequence from 1 .. N
			FragID_List frag_ids;
			flatten_list( fill_frames, frag_ids );
			if ( bRandom ) {
				numeric::random::random_permutation( frag_ids, numeric::random::rg() );
				numeric::random::random_permutation( frag_ids, numeric::random::rg() ); //playing safe
			}

			for ( auto it = frag_ids.begin(), eit = frag_ids.end();
					it != eit && nr_fill; ++it, --nr_fill ) {
				merged_frags->add( *it );
				// book keeping: raise counter for affected residues
				runtime_assert( it->frame().is_continuous() );
				for ( Size p = it->frame().start(); p<=it->frame().stop(); p++ ) nr_frags[ p ]++;

			}
			if ( nr_fill ) {
				tr.Info << nr_fill << " fragments short " << std::endl;
			} else {
				tr.Info << "succeeded! " << std::endl;
			}
		} //pos needs fill up
	}// loop over pos
	return merged_frags;
}


void
apply_best_scoring_fragdata(
	pose::Pose & pose,
	Frame const & frame,
	scoring::ScoreFunction const & sfxn
)
{

	if ( frame.nr_frags() < 1 ) return;

	frame.apply( 1, pose );

	core::Size best_frag(1);

	core::Real best_score( sfxn( pose ) );

	//tr << "frag 1 has score " << best_score << std::endl;

	for ( core::Size i = 2; i <= frame.nr_frags(); ++i ) {

		frame.apply( i, pose );

		core::Real cur_score( sfxn( pose ) );
		//tr << "frag " << i << " has score " << cur_score << std::endl;

		if ( cur_score < best_score ) {
			best_frag = i;
			best_score = cur_score;
		}

	}

	frame.apply( best_frag, pose );
	//tr << "applying frag " << best_frag << std::endl;

	//tr << "pose score before exit is " << sfxn( pose ) << std::endl;

} //apply_best_scoring_fragdata

void dump_frames_as_pdb(
	pose::Pose const & pose,
	utility::vector1< FrameOP > const & frames,
	std::string const & filename,
	Size const start_frag
)
{

	//we need to make a copy of the pose to muck around with
	pose::Pose frame_pose = pose;
	Size model_count(1), atom_counter(0);

	std::ofstream outfile( filename.c_str() );
	outfile << "REMARK  666  Fragment set for outtag pose \n";

	//now let's go through every fragment of every frame and dump to pdb
	for ( auto const & frame : frames ) {

		for ( Size frag = start_frag; frag <= frame->nr_frags(); ++frag ) {

			frame->apply( frag, frame_pose );

			outfile << "MODEL" << I(9, model_count) << "\n";

			for ( Size rescount = frame->start(); rescount <= frame->end(); ++rescount ) {

				core::io::pdb::dump_pdb_residue( frame_pose.residue( rescount ), atom_counter, outfile );

			} //loop over residues of fragment

			outfile << "ENDMDL\n";

			model_count++;

		} //loop over all fragments to apply

	} //iterator over all frames

	outfile.close();

} //dump_frames_as_pdb


/// @details this is a little tricky: this should support functionality for both creating a frameset entirely from
/// @details the input pdb (i.e. with no prior information about frag data length or srfds ), but it should also
/// @details  be possible to pass in non_empty frames such that the information in the pdb will generate FragData
/// @details  objects that are compatible to the ones already in the passed in frames. hmpf
/// @details
bool fill_template_frames_from_pdb(
	pose::Pose const & pose,
	utility::vector1< FrameOP > const & template_frames,
	std::string const & filename
)
{

	//we rely on the multimodel pdb reader for file processing
	utility::vector1< pose::Pose > in_poses;
	core::import_pose::pose_from_file( in_poses, filename , core::import_pose::PDB_file);

	Size frame_counter(1);
	bool return_val( true );
	pose::Pose frag_pose = pose;

	for ( utility::vector1< pose::Pose >::const_iterator pose_it = in_poses.begin(); pose_it != in_poses.end(); ++pose_it ) {

		//first we need to figure out which residues we are dealing with
		Size pdb_start_res( pose_it->pdb_info()->number( 1 ) );
		Size pdb_stop_res( pose_it->pdb_info()->number( pose_it->size() ) );
		char pdb_res_chain( pose_it->pdb_info()->chain( 1 ) );
		char pdb_start_res_icode( pose_it->pdb_info()->icode( 1 ) );
		char pdb_stop_res_icode( pose_it->pdb_info()->icode( pose_it->size() ) );

		//safety check
		if ( pose_it->pdb_info()->chain( pose_it->size() ) != pdb_res_chain ) {
			utility_exit_with_message("PDB file containing fragments is corrupted, one model contains multiple chains");
		}

		Size frame_start = pose.pdb_info()->pdb2pose( pdb_res_chain, pdb_start_res, pdb_start_res_icode );
		Size frame_stop = pose.pdb_info()->pdb2pose( pdb_res_chain, pdb_stop_res, pdb_stop_res_icode );


		while ( pose_it->size() != template_frames[ frame_counter ]->length()
				&& frame_start != template_frames[ frame_counter ]->start()
				&& frame_stop != template_frames[ frame_counter]->stop() )
				{
			frame_counter++;
			if ( frame_counter > template_frames.size() ) return false; //means there were input poses not compatible with the requested frames.
		}

		frag_pose.copy_segment( template_frames[ frame_counter ]->length(), *pose_it, template_frames[ frame_counter ]->start(), 1 );
		template_frames[ frame_counter ]->steal( frag_pose );

	} //iterator over input poses


	return return_val;

} //get_frames_from_pdb

void read_std_frags_from_cmd( FragSetOP& fragset_large, FragSetOP& fragset_small ) { //needed for setup_broker_from_cmdline
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string frag_large_file( "NoFile" );
	std::string frag_small_file( "NoFile" );
	if ( option[ in::file::frag9 ].user() ) {
		frag_large_file  = option[ in::file::frag9 ]();
	}

	if ( option[ in::file::frag3 ].user() ) {
		frag_small_file  = option[ in::file::frag3 ]();
	}

	if ( frag_large_file != "NoFile" ) {
		// fragset_large_ = FragmentIO().read( frag_large_file );
		fragset_large = FragmentIO(
			option[ OptionKeys::abinitio::number_9mer_frags ](),
			option[ OptionKeys::frags::nr_large_copies ](),
			option[ OptionKeys::frags::annotate ]()
			).read_data( frag_large_file );
	}

	if ( frag_small_file != "NoFile" ) {
		fragset_small = FragmentIO(
			option[ OptionKeys::abinitio::number_3mer_frags ],
			1, //nr_copies
			option[ OptionKeys::frags::annotate ]
			).read_data( frag_small_file );
	}

	if ( option[ OptionKeys::abinitio::steal_3mers ]() || option[ OptionKeys::abinitio::steal_9mers ]() ) {
		// read native pose to get sequence
		pose::PoseOP native_pose( new pose::Pose );
		if ( option[ in::file::native ].user() ) {
			core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
			pose::set_ss_from_phipsi( *native_pose );
		} else {
			utility_exit_with_message(" can't steal natie fragments without in:file:native " );
		}
		tr.Info << " stealing fragments from native pose: ATTENTION: native pose has to be IDEALIZED!!! " << std::endl;
		//  utility_exit_with_message(" stealing fragments from pose: currently not supported! ask Oliver " );
		if ( option[ OptionKeys::abinitio::steal_9mers ]() ) {
			if ( !fragset_large ) fragset_large = FragSetOP( new ConstantLengthFragSet( 9 ) );
			steal_frag_set_from_pose( *native_pose, *fragset_large, core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), fragset_large->max_frag_length() ) ) ) );
		}
		if ( option[ OptionKeys::abinitio::steal_3mers ]() ) {
			if ( !fragset_small ) fragset_small = FragSetOP( new ConstantLengthFragSet( 3 ) );
			steal_frag_set_from_pose( *native_pose, *fragset_small, core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), fragset_small->max_frag_length() ) ) ) );
		}
	}

	if ( fragset_small && option[ OptionKeys::abinitio::dump_frags ]() ) { //diagnosis
		utility::io::ozstream dump_frag_small( "fragset_small.dump" );
		for ( ConstFrameIterator it=fragset_small->begin(), eit=fragset_small->end(); it!=eit; ++it ) {
			(*it)->show( dump_frag_small );
		}
	}
	if ( fragset_large && option[ OptionKeys::abinitio::dump_frags ]() ) { //diagnosis
		utility::io::ozstream dump_frag_large( "fragset_large.dump" );
		for ( ConstFrameIterator it=fragset_large->begin(), eit=fragset_large->end(); it!=eit; ++it ) {
			(*it)->show( dump_frag_large );
		}
	}
}


/// @brief given a JumpFrame with Up and DownJumpSRFDs as LAST SRFDs this will make a fold-tree compatible with the
/// Frame...   this is NOT GOOD for sampling, since it introduces cut-points outside of fragments
/// later for sampling: one could probably write a routine that looks if it can move existing Jumps in Fold-tree to
/// fit the FRAME ... if not it returns failure...

/// one little assumption: create frames always like this:
/// contigues piece JUMP contig. piec JUMP contig. piece.
/// then good candidates for  cutpoints are the last contig.piece residue before a jump
void make_simple_fold_tree_from_jump_frame( Frame const& frame, Size total_residue, kinematics::FoldTree& new_fold_tree ) {
	/// @brief how many actual jumps are in this Frame?
	utility::vector1< core::Size > ups;
	utility::vector1< core::Size > downs;
	utility::vector1< core::Size > cuts;

	runtime_assert( frame.nr_frags() >= 1 );
	FragData const& fragdata( frame.fragment( 1 ) );
	for ( Size i=1; i<=frame.length(); ++i ) {
		// KAB 2015-06-10: Fixing warnings:
		// error: expression with side effects will be evaluated despite being used as an operand to 'typeid' [-Werror,-Wpotentially-evaluated-expression]
		// typeid() doesn't seem to like owning pointers (because they are a class),
		//   so creating a reference
		const SingleResidueFragData& fragdata_residue_i = *fragdata.get_residue( i );

		if ( typeid( fragdata_residue_i ) == typeid( UpJumpSRFD ) ) {
			const SingleResidueFragData& fragdata_residue_iplus1 = *fragdata.get_residue( i + 1 );
			runtime_assert( i + 1 <= frame.length() && typeid( fragdata_residue_iplus1 ) == typeid( DownJumpSRFD ) );
			ups.push_back( frame.seqpos( i ) );
			downs.push_back( frame.seqpos( i+1 ) );

			//the following configuration is not impossible, but if it should be allowed we need to change our heuristic to find cut-points
			runtime_assert( i >= 2 && typeid( fragdata_residue_i ) != typeid( DownJumpSRFD ) );
			cuts.push_back( frame.seqpos( i-1 ) );

			// we have used up two positions...
			++i;
		}
	}
	ObjexxFCL::FArray2D_int jump_point( 2, ups.size(), 0 );
	ObjexxFCL::FArray1D_int cut_point( ups.size() );
	for ( Size i = 1; i <= ups.size() ; ++i ) {
		jump_point( 1, i ) = ups[ i ];
		jump_point( 2, i ) = downs[ i ];
		cut_point( i ) = cuts[ i ];
	}
	new_fold_tree.tree_from_jumps_and_cuts( total_residue, ups.size(), jump_point, cut_point );
}

void fragment_set_slice ( ConstantLengthFragSetOP & fragset, Size const & min_res, Size const & max_res ){

	utility::vector1< Size > slice_res;
	for ( Size n = min_res; n <= max_res; n++ ) slice_res.push_back( n );
	fragment_set_slice( fragset, slice_res );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void
fragment_set_slice( core::fragment::ConstantLengthFragSetOP & fragset,
	utility::vector1< core::Size > const & slice_res ) {

	using namespace core::fragment;

	Size const len( fragset->max_frag_length() );

	ConstantLengthFragSetOP fragset_new( new ConstantLengthFragSet );

	for ( Size n = 1; n <= (slice_res.size() - len + 1); n++ ) {

		Size const & pos = slice_res[ n ];

		FrameList frames;

		if ( pos > (fragset->max_pos()-len+1) ) {
			tr.Warning << "NO FRAGS FOR POSITION " << pos << std::endl;
			continue;
		}

		fragset->frames( pos, frames );

		// CURRENTLY ONLY WORKS FOR CONST FRAG LENGTH SETS!!!! ASSUMES ONE FRAME!!!
		debug_assert( frames.size() == 1 );

		FrameOP & frame( frames[1] );
		FrameOP frame_new( new Frame( n, len ) );

		for ( Size k = 1; k <= frame->nr_frags(); k++ ) {
			frame_new->add_fragment( frame->fragment_ptr( k ) );
		}

		fragset_new->add( frame_new );

	}

	fragset = fragset_new;

}


/// @brief Finds the fold tree boundaries to the left and right of <pos>.
void FindBoundaries(const core::kinematics::FoldTree& tree, core::Size pos, core::Size* left, core::Size* right) {
	using core::Size;
	debug_assert(left);
	debug_assert(right);

	Size lower_cut = 0;
	Size upper_cut = tree.nres();
	Size num_cutpoints = tree.num_cutpoint();
	for ( Size i = 1; i <= num_cutpoints; ++i ) {
		Size cutpoint = tree.cutpoint(i);

		// find the upper boundary (inclusive)
		if ( cutpoint >= pos && cutpoint < upper_cut ) {
			upper_cut = cutpoint;
		}

		// find the lower boundary (exclusive)
		if ( cutpoint < pos && cutpoint > lower_cut ) {
			lower_cut = cutpoint;
		}
	}

	// set output parameters
	*left = lower_cut + 1;
	*right = upper_cut;
}

core::kinematics::Stub getxform(numeric::xyzVector<core::Real> m1,
	numeric::xyzVector<core::Real> m2,
	numeric::xyzVector<core::Real> m3,
	numeric::xyzVector<core::Real> f1,
	numeric::xyzVector<core::Real> f2,
	numeric::xyzVector<core::Real> f3 ) {
	core::kinematics::Stub s;
	s.M = alignVectorSets(m1-m2, m3-m2, f1-f2, f3-f2);
	s.v = f2-s.M*m2;
	return s;
}

void xform_pose(core::pose::Pose& pose,
	const core::kinematics::Stub& s,
	core::Size sres,
	core::Size eres ) {
	using core::id::AtomID;
	if ( eres==0 ) eres = pose.size();
	for ( core::Size ir = sres; ir <= eres; ++ir ) {
		for ( core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia ) {
			AtomID aid(AtomID(ia,ir));
			pose.set_xyz(aid, s.local2global(pose.xyz(aid)));
		}
	}
}

// frags must be ordered by sequence position
void make_pose_from_frags( pose::Pose & pose, std::string sequence, utility::vector1<FragDataCOP> frags, bool chains ) {
	// clear the pose
	pose.clear();
	// generate pose
	core::pose::make_pose_from_sequence(pose, sequence, core::chemical::CENTROID );
	// idealize and extend pose
	for ( Size pos = 1; pos<=pose.size(); pos++ ) {
		core::conformation::idealize_position( pos, pose.conformation() );
	}
	Real const init_phi  ( -150.0 );
	Real const init_psi  (  150.0 );
	Real const init_omega(  180.0 );
	for ( Size pos = 1; pos<=pose.size(); pos++ ) {
		if ( pos != 1 ) pose.set_phi( pos,  init_phi );
		if ( pos != pose.size() ) pose.set_psi( pos,  init_psi );
		if ( ( pos != 1 ) && ( pos != pose.size() ) ) pose.set_omega( pos,  init_omega );
	}
	// insert the torsions from the fragments
	Size insert_pos = 1;
	for ( Size i = 1; i <= frags.size(); ++i ) {
		FragDataCOP frag = frags[i];
		frag->apply(pose, insert_pos, insert_pos + frag->size() - 1);
		insert_pos += frag->size();
	}

	// if single fragment, reorientation is not necessary so just return
	if ( frags.size() == 1 ) return;

	// construct fold tree
	core::kinematics::FoldTree tree(pose.fold_tree());
	Size total_size = 0;
	for ( Size i = 1; i < frags.size(); ++i ) {
		FragDataCOP frag_1 = frags[i];
		FragDataCOP frag_2 = frags[i+1];
		Size midpoint_1 = static_cast<Size>(ceil(frag_1->size() / 2.0)) + total_size;
		total_size += frag_1->size();
		Size midpoint_2 = static_cast<Size>(ceil(frag_2->size() / 2.0)) + total_size + 1;
		int jump_id = tree.new_jump(midpoint_1, midpoint_2, total_size);
		tree.set_jump_atoms(jump_id, "CA", "CA");
	}
	pose.fold_tree(tree);

	// Ensure that the FoldTree is left in a consistent state
	debug_assert(pose.fold_tree().check_fold_tree());

	// orient the fragments in the pose based on CA fragment data
	total_size = 0;
	for ( Size i = 1; i <= frags.size(); ++i ) {
		FragDataCOP frag = frags[i];
		// Construct stubs from the 3 central CA atoms
		Size central_residue = static_cast< Size >(ceil(frag->size() / 2.0));
		BBTorsionSRFDCOP f1 = utility::pointer::static_pointer_cast< BBTorsionSRFD const >(frag->get_residue(central_residue - 1));
		BBTorsionSRFDCOP f2 = utility::pointer::static_pointer_cast< BBTorsionSRFD const >(frag->get_residue(central_residue    ));
		BBTorsionSRFDCOP f3 = utility::pointer::static_pointer_cast< BBTorsionSRFD const >(frag->get_residue(central_residue + 1));
		numeric::xyzVector<Real> fa1(f1->x(), f1->y(), f1->z());
		numeric::xyzVector<Real> fa2(f2->x(), f2->y(), f2->z());
		numeric::xyzVector<Real> fa3(f3->x(), f3->y(), f3->z());
		Size midpoint_pose = central_residue + total_size;
		total_size += frag->size();
		numeric::xyzVector<Real> m1 = pose.residue(midpoint_pose - 1).xyz("CA");
		numeric::xyzVector<Real> m2 = pose.residue(midpoint_pose).xyz("CA");
		numeric::xyzVector<Real> m3 = pose.residue(midpoint_pose + 1).xyz("CA");
		// Compute the transform
		core::kinematics::Stub x = getxform(m1, m2, m3, fa1, fa2, fa3);
		// Starting at the midpoint of the insertion point in the pose, propagate
		// the change to the left and right until we reach either the end of the
		// chain or a cutpoint
		core::Size region_start, region_stop;
		FindBoundaries(tree, midpoint_pose, &region_start, &region_stop);
		xform_pose(pose, x, region_start, region_stop);
		if ( chains && i < frags.size() ) pose.conformation().insert_chain_ending( total_size );
	}
}

} //fragment
} //core
