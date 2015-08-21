// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/loops_main.cc
/// @brief loop building tools
/// @author Mike Tyka
/// @author Chu Wang
/// @author Daniel J. Mandell

// Unit headers
#include <protocols/loops/loops_main.hh>

// Package headers
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Project headers
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/util.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/util.hh>
#include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/types.hh>

// Basic headers
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>

// C++ Headers
#include <iostream>
#include <map>
#include <string>

// ObjexxFCL includes
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace protocols {
namespace loops {

///////////////////////////////////////////////////////////////////////////////
using namespace core;

static thread_local basic::Tracer tt( "protocols.loops.loops_main" );


void read_loop_fragments(
	std::vector< core::fragment::FragSetOP > & frag_libs
) {
	using namespace basic::options;
	using namespace utility::file;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::fragment;

	utility::vector1<int> frag_sizes( option[ OptionKeys::loops::frag_sizes ] );
	FileVectorOption      frag_files( option[ OptionKeys::loops::frag_files ] );

	if ( frag_sizes.size() != frag_files.size() ) {
		utility_exit_with_message( "You must specify as many fragment sizes as fragment file names " );
	}

	for ( Size i = 1; i <= frag_sizes.size(); ++i ) {
		Size const frag_size = Size(frag_sizes[i]);

		FragSetOP frag_lib_op( new ConstantLengthFragSet( frag_size ) );
		//protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
		tt.Error << "Frag libraries debug " << frag_files[i] << " " << frag_size << std::endl;

		if ( frag_files[i] != std::string("none") ) {
			//frag_lib_op->read_fragment_file( frag_files[i]  );
			frag_lib_op = (FragmentIO().read_data( frag_files[i] ));
		}
		frag_libs.push_back(frag_lib_op);
	}

	Size prev_size(10000);
	FragSetOP prev_lib_op(0);

	// Loop over the temporary map to generate missing/noninitialized fragment
	// libraries
	for ( std::vector< FragSetOP >::const_iterator
			it = frag_libs.begin(),
			it_end = frag_libs.end();
			it != it_end; ++it
			) {
		Size const frag_size( (*it)->max_frag_length() );
		Size const n_frags( (*it)->size() );

		if ( frag_size > prev_size ) {
			//std::cerr << frag_size << "  " << prev_size << std::endl;
			std::string msg;
			msg += " frag_size = " + string_of( frag_size );
			msg += " prev_size = " + string_of( prev_size );
			msg += "\nFragment size must be given in order !!\n";
			utility_exit_with_message( msg );
		}

		if ( (n_frags == 0) && prev_lib_op && (prev_lib_op->size() != 0) ) {
			tt.Info << "Set up " << frag_size << "-mer library from " << prev_size << "-mer library" << std::endl;

			chop_fragments( *prev_lib_op, **it );
		}
		prev_size = frag_size;
		prev_lib_op = *it;
	}


	// Steal fragments as requested

	if ( option[ OptionKeys::loops::stealfrags ].user() ) {
		utility::vector1< FileName > pdbfiles = option[  OptionKeys::loops::stealfrags ]();
		utility::vector1< FileName >::iterator file_it = pdbfiles.begin(), file_it_end = pdbfiles.end();
		// Loop over all the islent input files
		for ( ; file_it != file_it_end; ++file_it ) {
			std::string infile  = *file_it;

			core::pose::Pose stealpose;
			core::import_pose::centroid_pose_from_pdb( stealpose, infile );

			for ( std::vector< FragSetOP >::const_reverse_iterator
					it = frag_libs.rbegin(),
					it_end = frag_libs.rend();
					it != it_end; ++it ) {
				tt.Info << "Stealing fragments from " << infile << "  "
					<< option[  OptionKeys::loops::stealfrags_times ]() << " times" << std::endl;

				for ( int c=0; c< option[  OptionKeys::loops::stealfrags_times ](); c++ ) {
					//steal_constant_length_frag_set_from_pose ( stealpose, **it );
					steal_frag_set_from_pose( stealpose, **it, core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( core::fragment::SingleResidueFragDataOP( new BBTorsionSRFD ), (*it)->max_frag_length() ) ) ) );
				}
			}
		} // loop over input files
	} // if stealfrags

	for ( std::vector< FragSetOP >::const_reverse_iterator
			it = frag_libs.rbegin(),
			it_end = frag_libs.rend();
			it != it_end; ++it ) {

		Size const frag_size( (*it)->max_frag_length() );
		Size const n_frags( (*it)->size() );
		tt.Info << "Fragment libraries: " << frag_size << "   " << n_frags
			<< std::endl;
	}
}


void read_loop_fragments(
	utility::vector1< core::fragment::FragSetOP > & frag_libs
) {
	using std::vector;
	using utility::vector1;

	std::vector< core::fragment::FragSetOP > temp_libs;
	read_loop_fragments( temp_libs );
	frag_libs.resize( temp_libs.size() );
	for ( Size i = 1; i <= temp_libs.size(); ++i ) {
		frag_libs[i] = temp_libs[i-1];
	}
}


//////////////////////////////////////////////////////////////////////////////////////////
/// @details for each loop defined, add a fixed jump from the residue before loop_start and
/// the residue after the loop_end. The cutpoint is from loop_cut. Support terminal loops, i.e.,
/// starting at residue 1 or ending at total_residue
////////////////////////////////////////////////////////////////////////////////////////////
void
fold_tree_from_loops(
	core::pose::Pose const & pose,
	Loops const & loops,
	kinematics::FoldTree & f,
	bool terminal_cutpoint
)
{
	using namespace kinematics;

	f.clear();

	// "symmetry-safe" version
	FoldTree const &f_in = core::conformation::symmetry::get_asymm_unit_fold_tree( pose.conformation() );

	// nres points to last protein residue;
	Size totres = f_in.nres();
	Size nres = totres;
	if ( nres != pose.total_residue() ) nres--;   // only true if pose is symm. ...  asymm foldtree is then rooted on VRT

	// following residues (e.g. ligands) will be attached by jumps
	while ( ! ( pose.residue( nres ).is_protein() || pose.residue( nres ).is_carbohydrate() ) ) {
		nres -= 1;
	}

	Loops tmp_loops = loops;
	tmp_loops.sequential_order();

	Size jump_num = 0;
	Size prev_interchain_jump = 0;
	for ( Loops::const_iterator it=tmp_loops.begin(), it_end=tmp_loops.end(), it_next; it != it_end; ++it ) {
		it_next = it;
		it_next++;

		bool is_lower_term = pose.residue( it->start() ).is_lower_terminus();
		bool is_upper_term = pose.residue( it->stop() ).is_upper_terminus();

		Size jump_start = it->start()-1;
		Size jump_stop  = it->stop()+1;
		Size const jump_cut   = it->cut();

		if ( jump_cut == 0 ) {
			utility_exit_with_message("Can't build a fold tree from a loop with an unspecified cut point.");
		}

		Size const jump_next_start = ( it_next == it_end ) ? nres : it_next->start() - 1;

		if ( jump_stop <= nres && jump_stop < jump_next_start ) f.add_edge(jump_stop, jump_next_start, Edge::PEPTIDE);

		// If this is a terminal loop and the terminus should be used for the cut instead of whatever
		// the loop's cutpoint is set to (i.e. terminal_cutpoint is false, which is the default by the
		// way), just make an Edge that spans the loop.
		if ( (is_lower_term || is_upper_term) && !terminal_cutpoint ) {
			if ( is_lower_term && is_upper_term ) {
				utility_exit_with_message("Trying to make a loop that covers an entire chain. No anchor position possible!");
			}

			if ( is_lower_term ) {
				f.add_edge(jump_start+1, jump_stop, Edge::PEPTIDE);
				if ( prev_interchain_jump > 0 ) {
					f.add_edge(prev_interchain_jump, jump_stop, ++jump_num);
				} else if ( jump_start>0 ) {
					f.add_edge(jump_start, jump_stop, ++jump_num);
				}
				prev_interchain_jump = 0;
			} else { // is_upper_term
				prev_interchain_jump = jump_start;
				f.add_edge(jump_start, jump_stop-1, Edge::PEPTIDE);

				// if we're not rebuilding the adjacent n-term, add a jump to the next segment
				if ( (it_next != it_end && it_next->start() != jump_stop) || ( it_next == it_end && jump_stop != nres+1) ) {
					f.add_edge(jump_start, jump_stop, ++jump_num);
				}
			}
			continue;
		}
		prev_interchain_jump = 0;

		if ( jump_start==0 ) jump_start=1;
		if ( jump_stop>nres ) jump_stop=nres;

		f.add_edge( jump_start, jump_cut,  Edge::PEPTIDE );

		// Increase the jump number and add a jump connecting the residues around the loop
		++jump_num;
		f.add_edge(jump_start, jump_stop, jump_num);

		// Add edges directed from the terminal loop residues to the cutpoint
		f.add_edge(jump_cut + 1, jump_stop, Edge::PEPTIDE);
	}

	// This is kind of dumb -- why don't we bail way eariler if there are no loops?
	if ( tmp_loops.size() > 0 ) {
		Size first_start = tmp_loops.begin()->start();
		first_start = (first_start == 1) ? first_start : first_start - 1;
		f.add_edge( 1, first_start, Edge::PEPTIDE );

		// reorder
		Size root;
		if ( tmp_loops.begin()->start() == 1 && tmp_loops.begin()->stop() != nres ) {
			root = tmp_loops.begin()->stop()+1;
		} else {
			root = 1;
		}

		// Test for successful reordering.
		if ( !f.reorder(root) ) {
			throw utility::excn::EXCN_Msg_Exception("Unable to reorder the FoldTree for this loops set!");
		}
	}

	// Attach remaining (non-protein) residues by jumps to the tree root.
	Size jump_anchor = f.root();
	while ( nres < totres ) {
		nres += 1;
		f.add_edge( jump_anchor, nres, f.num_jump()+1 );
	}

	if ( pose.residue( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) {
		// special case for fold trees rooted on a VRT (i.e. symmetry)
		if ( f_in.nres() != pose.total_residue() ) {
			f.reorder( f_in.nres() );
		} else {
			f.reorder( pose.fold_tree().root() );
		}
	}

	// symmetrize the fold tree (which is over asymm unit)
	core::pose::symmetry::symmetrize_fold_tree( pose, f );
}


//////////////////////////////////////////////////////////////////////////////////
/// @details  Make a single fold tree that brackets the loop
void set_single_loop_fold_tree(
	core::pose::Pose & pose,
	Loop const & loop
)
{

	//setup fold tree for this loop
	kinematics::FoldTree f;
	Loops loops;
	loops.add_loop(loop);
	fold_tree_from_loops(pose, loops, f);
	tt.Warning << "Pose fold tree " << f << std::endl;
	pose.fold_tree(f);

}


/// @details  Slide a loop cutpoint to a (potentially) new position
/// @note  Updates the pose's foldtree, either moving or adding a loop cutpoint
/// @note  Updates CUTPOINT_UPPER and CUTPOINT_LOWER variant status of residues in loop to match new cutpoint location
void
set_loop_cutpoint_in_pose_fold_tree(
	Size const new_cutpoint,
	pose::Pose & pose,
	Size const loop_begin,
	Size const loop_end
)
{
	using namespace chemical;
	kinematics::FoldTree f( pose.fold_tree() );
	if ( f.is_cutpoint( new_cutpoint ) ) return;

	// find the current cutpoint
	Size cut(0);
	for ( Size i=loop_begin-1; i<= loop_end; ++i ) {
		if ( f.is_cutpoint( i ) ) {
			if ( cut ) utility_exit_with_message( "multiple cutpoints in single loop!" );
			cut = i;
		}
		// remove cutpoint variants if they're present
		core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, i   );
		core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, i+1 );
	}
	if ( !cut ) {
		tt.Warning << "set_loop_cutpoint_in_pose_fold_tree: no cutpoint in loop, so adding new jump to foldtree between " <<
			loop_begin-1 << " and " << loop_end +1 << " with cut at " << new_cutpoint << std::endl;
		f.new_jump( loop_begin-1, loop_end+1, new_cutpoint );
	} else {
		f.slide_cutpoint( cut, new_cutpoint );
	}
	pose.fold_tree( f );

	core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, new_cutpoint   );
	core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, new_cutpoint+1 );

}


//////////////////////////////////////////////////////////////////////////////////
/// @details  Remove cutpoint variants.
void remove_cutpoint_variants(
	core::pose::Pose & pose,
	bool force
)
{
	using namespace core::chemical;
	bool pose_changed( false );
	pose::Pose init_pose = pose;
	if ( force ) {
		for ( core::Size ir=1; ir<=pose.total_residue() ; ++ir ) {
			if ( pose.residue(ir).has_variant_type(CUTPOINT_LOWER) ) {
				pose_changed = true;
				core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, ir );
			}
			if ( pose.residue(ir).has_variant_type(CUTPOINT_UPPER) ) {
				core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, ir );
				pose_changed = true;
			}
		}
	} else {

		for ( int i=1; i<=pose.fold_tree().num_cutpoint() ; ++i ) {
			int cutpoint = pose.fold_tree().cutpoint(i);
			if ( pose.residue(cutpoint).has_variant_type(CUTPOINT_LOWER) ) {
				pose_changed = true;
				core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, cutpoint );
			}
			if ( pose.residue(cutpoint+1).has_variant_type(CUTPOINT_UPPER) ) {
				core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, cutpoint+1 );
				pose_changed = true;
			}
		}

		//WTF ?
		if ( pose.residue(1).has_variant_type(CUTPOINT_LOWER) ) {
			pose_changed = true;
			core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, 1 );
		}
		if ( pose.residue(2).has_variant_type(CUTPOINT_UPPER) ) {
			core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, 2 );
			pose_changed = true;
		}
		if ( pose.residue(pose.total_residue()-1).has_variant_type(CUTPOINT_LOWER) ) {
			pose_changed = true;
			core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, pose.total_residue()-1 );
		}
		if ( pose.residue(pose.total_residue()).has_variant_type(CUTPOINT_UPPER) ) {
			core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, pose.total_residue() );
			pose_changed = true;
		}
	}

	//I am pretty sure this assignment was a bug, and it was meant to be an equality operator - 4/18/11 SML
	//if ( pose_changed = true ) pose.constraint_set( init_pose.constraint_set()->remapped_clone( init_pose, pose ) );
	if ( pose_changed == true ) pose.constraint_set( init_pose.constraint_set()->remapped_clone( init_pose, pose ) );
}


//////////////////////////////////////////////////////////////////////////////////
// Add cutpoint variants.
void
add_cutpoint_variants( core::pose::Pose & pose )
{
	for ( int i = 1; i <= pose.fold_tree().num_cutpoint() ; ++i ) {
		core::uint cutpoint = pose.fold_tree().cutpoint( i );
		add_single_cutpoint_variant( pose, cutpoint );
	}
}


void
add_single_cutpoint_variant( core::pose::Pose & pose, const Loop & loop )
{
	add_single_cutpoint_variant( pose, loop.cut() );
}


void
add_single_cutpoint_variant( core::pose::Pose & pose, const core::uint cutpoint )
{
	using namespace std;
	using namespace core::chemical;
	using namespace core::pose;

	Pose init_pose = pose;
	bool pose_changed( false );

	// Cutpoint and terminus variant types are currently incompatible, so check for presence.
	// Also, the logic for the chainbreak energy method only works with peptides and carbohydrates.
	if ( ( pose.residue( cutpoint ).is_protein() || pose.residue( cutpoint ).is_carbohydrate() ) &&
			!( pose.residue( cutpoint ).is_upper_terminus() ) ) {
		add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, cutpoint );
		tt.Info << "Added cutpoint variant to residue " << cutpoint << endl;
		pose_changed = true;
	} else {
		tt.Warning << "Residue " << cutpoint <<
			" is not compatible with cutpoints; variant type not changed." << endl;
	}

	if ( ( pose.residue( cutpoint + 1 ).is_protein() || pose.residue( cutpoint + 1 ).is_carbohydrate() ) &&
			!( pose.residue( cutpoint + 1 ).is_lower_terminus() ) ) {
		add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, cutpoint + 1 );
		tt.Info << "Added cutpoint variant to residue " << cutpoint + 1 << endl;
		pose_changed = true;
	} else {
		tt.Warning << "Residue " << cutpoint + 1 <<
			" is not compatible with cutpoints; variant type not changed." << endl;
	}

	if ( pose_changed ) {
		pose.constraint_set( init_pose.constraint_set()->remapped_clone( init_pose, pose ) );
	}
}


//////////////////////////////////////////////////////////////////////////////////////
/// @details omega backbone torsion is always fixed. phi/psi backbone torsions within
/// the loop region are flexible. Depending on whether -fix_natsc flag, sidechain DOFs
/// of loop residues and/or their neighboring residues in the template will be set as
/// movable. Default neighbors are 10A CB dist from loop residues; neighbor_dist can
/// be used to further filter the neighbors. This is a wrapper function which determine
/// moveable sidechains and then call actual loop_set_move_map function to set up move
/// map properly
//////////////////////////////////////////////////////////////////////////////////////
void
loops_set_move_map(
	pose::Pose & pose,
	Loops const & loops,
	bool const fix_template_sc,
	core::kinematics::MoveMap & mm,
	Real neighbor_dist
)
{
	using namespace basic::options;
	loops_set_move_map(
		pose, loops, fix_template_sc, mm, neighbor_dist,
		option[OptionKeys::loops::allow_omega_move].user(),
		option[OptionKeys::loops::allow_takeoff_torsion_move].user());
}

void
loops_set_move_map(
	pose::Pose & pose,
	Loops const & loops,
	bool const fix_template_sc,
	core::kinematics::MoveMap & mm,
	Real neighbor_dist,
	bool const allow_omega_move,
	bool const allow_takeoff_torsion_move
)
{
	using namespace core::id;
	pose.update_residue_neighbors();

	utility::vector1<bool> allow_sc_move( pose.total_residue(), false);
	select_loop_residues( pose, loops, !fix_template_sc, allow_sc_move, neighbor_dist);
	loops_set_move_map( loops, allow_sc_move,mm, allow_omega_move, allow_takeoff_torsion_move);

	//fpd symmetric version
	if ( core::pose::symmetry::is_symmetric( pose ) )  {
		core::pose::symmetry::make_symmetric_movemap( pose, mm );
	}
}

void
loops_set_move_map(
	Loops const & loops,
	utility::vector1<bool> const & allow_sc_move,
	core::kinematics::MoveMap & mm
)
{
	using namespace basic::options;
	loops_set_move_map(
		loops, allow_sc_move, mm,
		option[OptionKeys::loops::allow_omega_move].user(),
		option[OptionKeys::loops::allow_takeoff_torsion_move].user());
}

void
loops_set_move_map(
	Loops const & loops,
	utility::vector1<bool> const & allow_sc_move,
	core::kinematics::MoveMap & mm,
	bool const allow_omega_move,
	bool const allow_takeoff_torsion_move
)
{
	using namespace core::id;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// allow chi to move
	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );
	// allow phi/psi in loops to move
	for ( Loops::const_iterator it=loops.begin(), it_end=loops.end();
			it != it_end; ++it ) {

		for ( Size i=it->start(); i<=it->stop(); ++i ) {
			mm.set_bb(i, true);
			if ( !allow_omega_move ) mm.set(TorsionID(i,BB,3), false); // omega is fixed by default
			mm.set_chi(i, true); // chi of loop residues
		}

		if ( allow_takeoff_torsion_move ) {
			mm.set( TorsionID( it->start()-1, BB, 2 ), true );
			if ( allow_omega_move ) mm.set( TorsionID( it->start()-1, BB, 3 ), true );
			mm.set( TorsionID( it->stop()+1, BB, 1 ), true );
		}

	}
	// set chi move map based on the input allow_sc array, which is filled based on fix_natsc info
	for ( Size i = 1; i <= allow_sc_move.size(); ++i ) {
		mm.set_chi(i, allow_sc_move[i] );
	}
	// chu in case we have ligand attached in fold tree and they need to move. Assuming that
	// loop jumps come first followed by lig jumps.
	if ( basic::options::option[ OptionKeys::loops::allow_lig_move ]() ) {
		mm.set_jump( true );
		for ( Size i = 1; i <= loops.num_loop(); ++i ) {
			mm.set_jump( i, false );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @details Wrapper around loops_set_move_map which works on a single loop
/// object and returns a MoveMapOP.  See the detailed documentation for
/// loops_set_move_map for more information about the generated MoveMap.
///////////////////////////////////////////////////////////////////////////////

core::kinematics::MoveMapOP
move_map_from_loops(
	core::pose::Pose & pose,
	Loops const & loops,
	bool const fix_template_sc,
	core::Real neighbor_dist,
	bool const flanking_residues
)
{
	using core::kinematics::MoveMap;
	using core::kinematics::MoveMapOP;

	MoveMapOP move_map( new MoveMap );
	loops_set_move_map(pose, loops, fix_template_sc, *move_map, neighbor_dist);

	if ( flanking_residues ) {
		add_loop_flank_residues_bb_to_movemap(loops, *move_map);
	}

	return move_map;
}

///////////////////////////////////////////////////////////////////////////////
/// @details Wrapper around loops_set_move_map which works on a single loop
/// object and returns a MoveMapOP.  See the detailed documentation for
/// loops_set_move_map for more information about the generated MoveMap.
///////////////////////////////////////////////////////////////////////////////

core::kinematics::MoveMapOP
move_map_from_loop(
	core::pose::Pose & pose,
	Loop const & loop,
	bool const fix_template_sc,
	core::Real neighbor_dist,
	bool const flanking_residues
)
{
	Loops loops = Loops();
	loops.add_loop(loop);

	return move_map_from_loops(
		pose, loops, fix_template_sc, neighbor_dist, flanking_residues);
}

//////////////////////////////////////////////////////////////////////////////////////
/// @details omega backbone torsion is always fixed. phi/psi backbone torsions within
/// the loop region are flexible. Depending on whether -fix_natsc flag, sidechain DOFs
/// of loop residues and/or their neighboring residues in the template will be set as
/// movable.
//////////////////////////////////////////////////////////////////////////////////////
void
set_move_map_for_centroid_loop(
	Loop const & loop,
	core::kinematics::MoveMap & mm
)
{
	using namespace core::id;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	bool const allow_omega_move = basic::options::option[ OptionKeys::loops::allow_omega_move ].user();
	bool const allow_takeoff_torsion_move = basic::options::option[ OptionKeys::loops::allow_takeoff_torsion_move ].user();

	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );
	// allow phi/psi in loops to move
	for ( Size i=loop.start(); i<=loop.stop(); ++i ) {
		mm.set_bb(i, true);
		if ( !allow_omega_move ) mm.set( TorsionID( i, BB, 3 ), false ); // omega is fixed
		mm.set_chi(i, true); // chi of loop residues
	}

	if ( allow_takeoff_torsion_move ) {
		mm.set( TorsionID( loop.start()-1, BB, 2 ), true );
		if ( allow_omega_move ) mm.set( TorsionID( loop.start()-1, BB, 3 ), true );
		mm.set( TorsionID( loop.stop()+1, BB, 1 ), true );
	}

	// chu in case we have ligand attached in fold tree and they need to move. Assuming that
	// loop jumps come first followed by lig jumps.
	if ( basic::options::option[ OptionKeys::loops::allow_lig_move ].user() ) {
		mm.set_jump( true );
		mm.set_jump( 1, false ); //jump for the single loop
	}


}

//////////////////////////////////////////////////////////////////////////////////////
///// @details When building loops on a homology model, the qualities of the loop stems are
///// hard to control. Implementing small/shear/ccd movers to them may be too much. One
///// may just want to minimize them when the loop refinement protocl is doing minimization.
///// This is the function to add these extra residues to the movemap.  --JQX
////////////////////////////////////////////////////////////////////////////////////////
void
add_loop_flank_residues_bb_to_movemap(
	Loops const & loops,
	core::kinematics::MoveMap & mm,
	core::Size flank_size
){

	for ( Loops::const_iterator it=loops.begin(), it_end=loops.end(); it != it_end; ++it ) {

		for ( Size i=(it->start()-flank_size); i<=(it->start()-1); i++ ) {
			mm.set_bb(i, true);
		}

		for ( Size i=(it->stop()+1); i<=(it->stop()+flank_size); i++ ) {
			mm.set_bb(i, true);
		}

	}

}


///////////////////////////////////////////////////////////////////////////////////////////////////
/// @details This function takes a pose and loop definitions and closes each loop separately by the CCD algorithm.
/// All parameters related to CCD are hard-coded.  If one wishes more control over the options, s/he should use the
/// CCDLoopClosureMover.
void
ccd_close_loops(
	pose::Pose & pose,
	Loops const & loops,
	kinematics::MoveMap const & mm )
{
	loop_closure::ccd::CCDLoopClosureMover ccd_loop_closure_mover;
	ccd_loop_closure_mover.movemap( kinematics::MoveMapCOP( kinematics::MoveMapOP( new kinematics::MoveMap( mm ) ) ) );

	for ( Loops::const_iterator it = loops.begin(), it_end = loops.end(); it != it_end; ++it ) {
		ccd_loop_closure_mover.loop( *it );
		ccd_loop_closure_mover.apply( pose );
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details use TenANeighborGraph. As input, residue_positions[i] is true for residues to be counted.
/// As output, residue_position[i] is true for all neighbor residues including orginal input residues.
/// The function is used to find all neighboring residues of the loop residues in case they need to be
/// repacked or minimized in fullatom refinement.
void get_tenA_neighbor_residues(
	pose::Pose const & pose,
	utility::vector1<bool> & residue_positions
)
{
	//make a local copy first because we will change content in residue_positions
	utility::vector1<bool> local_residue_positions = residue_positions;
	core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
	for ( Size i=1; i <= local_residue_positions.size(); ++i ) {
		if ( ! local_residue_positions[i] ) continue;
		core::graph::Node const * current_node( tenA_neighbor_graph.get_node(i)); // find neighbors for this node
		for ( core::graph::Node::EdgeListConstIter it = current_node->const_edge_list_begin();
				it != current_node->const_edge_list_end(); ++it ) {
			Size pos = (*it)->get_other_ind(i);
			if ( pose.residue(pos).type().is_disulfide_bonded() ) {
				residue_positions[ pos ] = false;
			} else {
				residue_positions[ pos ] = true;
			}
		}
	}
}
/////////////////////////////////////////////////////////////////////////////
/// @details use 10A CB distance cutoff as neighboring residue defintion. The function
///is used for conveniently setting up sidechain movable residues in loop modeling.
///The 10A residue set is further reduced if neighbor_dist < 10.0
/////////////////////////////////////////////////////////////////////////////////
void select_loop_residues(
	pose::Pose const & pose,
	Loops const & loops,
	bool const include_neighbors,
	utility::vector1<bool> & map,
	Real neighbor_dist
)
{
	for ( Loops::const_iterator it=loops.begin(), it_end=loops.end();
			it != it_end; ++it ) {
		for ( Size i=it->start(); i<=it->stop(); ++i ) {
			if ( pose.residue(i).type().is_disulfide_bonded() ) {
				map[i] = false;
			} else {
				map[i] = true;
			}
		}
	}

	if ( include_neighbors ) get_tenA_neighbor_residues( pose, map );
	// if the neighbor_dist is less than 10A, filter the 10A neighbors to that distance
	if ( neighbor_dist < 10.0 ) {
		filter_loop_neighbors_by_distance( pose, map, loops, neighbor_dist );
	}
	return;
}

/////////////////////////////////////////////////////////////////////////////
/// @details use 10A CB distance cutoff as neighboring residue defintion. The function
///is used for conveniently setting up sidechain movable residues in loop modeling.
///The 10A residue set is further reduced if neighbor_dist < 10.0
utility::vector1<bool> select_loop_residues(
	pose::Pose const & pose,
	Loops const & loops,
	bool const include_neighbors,
	Real neighbor_dist
)
{
	utility::vector1<bool> map(pose.total_residue(), false);
	select_loop_residues(pose, loops, include_neighbors, map, neighbor_dist);
	return map;
}

/////////////////////////////////////////////////////////////////////////////
/// @details for one loop only
void select_loop_residues(
	pose::Pose const & pose,
	Loop const & loop,
	bool const include_neighbors,
	utility::vector1<bool> & map,
	Real neighbor_dist
)
{
	Loops loops;
	loops.add_loop( loop );
	select_loop_residues( pose, loops, include_neighbors, map, neighbor_dist );
	return;
}

/////////////////////////////////////////////////////////////////////////////
/// @details for one loop only
utility::vector1<bool> select_loop_residues(
	pose::Pose const & pose,
	Loop const & loop,
	bool const include_neighbors,
	Real neighbor_dist
)
{
	utility::vector1<bool> map(pose.total_residue(), false);
	select_loop_residues(pose, loop, include_neighbors, map, neighbor_dist);
	return map;
}

//////////////////////////////////////////////////////////////////////////////////////
/// @details neighbors contains the set of potential neighbors to the loop residues
///given in loops. This set is reduced to only contain neighbors within dist_cutoff
///of any residue in loops.
//////////////////////////////////////////////////////////////////////////////////////
void filter_loop_neighbors_by_distance(
	pose::Pose const & pose,
	utility::vector1<bool> & map,
	Loops const & loops,
	Real & dist_cutoff
)
{
	//Size orig_neighbs=0, new_neighbs=0;
	for ( Size i=1; i<=map.size(); ++i ) {
		if ( map[i] == false ) continue; // this residue already isn't a neighbor
		//++orig_neighbs;
		map[i] = false;
		for ( Loops::const_iterator it=loops.begin(), it_end=loops.end(); it != it_end; ++it ) {
			for ( Size j=it->start(); j<=it->stop(); ++j ) {
				// Get the atom vectors for loop and scaffold CB, or CA if GLY
				numeric::xyzVector< Real > scaffold_res_vec;
				numeric::xyzVector< Real > loop_res_vec;
				scaffold_res_vec = pose.residue( i ).xyz( pose.residue( i ).nbr_atom() );
				loop_res_vec = pose.residue( j ).xyz( pose.residue( j ).nbr_atom() );
				// only keep as neighbor if dist within cutoff
				Real dist = scaffold_res_vec.distance( loop_res_vec );
				if ( dist <= dist_cutoff ) {
					map[i] = true;
				}
			}
		}
	}

	// some debug code to see the effect of this code, requires uncommenting orig_neighbs and new_neighbs above
	//for( Size z=1; z<=map.size(); ++z ) {
	// if( map[z]==1 ) ++new_neighbs;
	//}
	//tt << "Number of neighbors reduced from " << orig_neighbs << " to " << new_neighbs << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////
// if the pose sequence extends beyond the sequence from the mapping, this will extend the mapping
// semi-special-case-HACK, logic needs clarifying
//
void
extend_sequence_mapping(
	pose::Pose const & pose,
	id::SequenceMapping & mapping,
	std::string & source_seq,
	std::string & target_seq
)
{
	using namespace conformation;

	runtime_assert( mapping.size1() == source_seq.size() && mapping.size2() == target_seq.size());

	// may not represent the entire pose source sequence
	std::string const pose_seq( pose.sequence() );

	if ( source_seq != pose_seq ) {
		/// source sequence from align file does not cover entire pdb file sequence

		if ( pose_seq.find( source_seq ) == std::string::npos ) {
			tt << "src_seq:  " << source_seq << std::endl << "pose_seq: " << pose_seq << std::endl;
			utility_exit_with_message( "alignfile source sequence not contained in pose sequence" );
		}
		Size const offset( pose_seq.find( source_seq ) );

		std::string source_nterm_seq, target_nterm_seq;
		for ( Size i=1; i<= offset; ++i ) {
			// this is a residue in the input pdb that's not accounted for in the alignment
			// if its protein, unalign it
			// if its dna, align it
			Residue const & rsd( pose.residue( i ) );

			mapping.insert_source_residue( i ); // remains unaligned, ie mapping[i] == 0
			source_nterm_seq += rsd.name1();

			bool const align_this_residue( !rsd.is_protein() );

			if ( align_this_residue ) {
				target_nterm_seq += rsd.name1();
				Size const target_pos( target_nterm_seq.size() );
				mapping.insert_target_residue( target_pos );
				mapping[i] = target_pos;
			}
		}

		source_seq = source_nterm_seq + source_seq;
		target_seq = target_nterm_seq + target_seq;

		runtime_assert( pose_seq.find( source_seq ) == 0 );

		while ( source_seq.size() < pose_seq.size() ) {
			runtime_assert( mapping.size1() == source_seq.size() && mapping.size2() == target_seq.size());

			Size const pos( source_seq.size() + 1 );
			Residue const & rsd( pose.residue( pos ) );

			char const n1( rsd.name1() );
			source_seq += n1;
			mapping.push_back( 0 );

			bool const align_this_residue( !rsd.is_protein() );
			if ( align_this_residue ) {
				target_seq += n1;
				mapping.insert_target_residue( target_seq.size() );
				mapping[pos] = target_seq.size();
			}
		}

		runtime_assert( pose_seq == source_seq );
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// what can we assume about the starting fold_tree??
void
apply_sequence_mapping(
	pose::Pose & pose,
	std::string const & target_seq,
	id::SequenceMapping const & start_mapping
)
{
	using namespace conformation;
	using namespace chemical;

	id::SequenceMapping mapping( start_mapping );
	runtime_assert( mapping.size1() == pose.total_residue() && mapping.size2() == target_seq.size() );

	tt << "apply_sequence_mapping: start fold tree: " << pose.fold_tree() << std::endl;

	// try skipping this hacky rebuild of a new tree
	if ( false && pose.num_jump() ) { // check the current fold_tree
		Size const num_jump( pose.num_jump() );
		Size const nres( pose.total_residue() );

		Size segment(1);
		utility::vector1< int > seg_anchor( num_jump+1, 0 ), seg_end;
		for ( Size i=1; i<= nres; ++i ) {
			if ( mapping[i] ) {
				if ( seg_anchor[ segment ] == 0 ) {
					seg_anchor[ segment ] = i;
				}
			}
			if ( i < nres && pose.fold_tree().is_cutpoint( i ) ) {
				seg_end.push_back( i );
				++segment;
			}
		}

		FArray2D_int jumps( 2, num_jump );
		FArray1D_int cuts ( num_jump );
		for ( Size i=1; i<= num_jump; ++i ) {
			jumps(1,i) = seg_anchor[1];
			jumps(2,i) = seg_anchor[i+1];
			cuts(i) = seg_end[i];
		}

		kinematics::FoldTree f;
		bool valid_tree = f.tree_from_jumps_and_cuts( nres, num_jump, jumps, cuts );
		runtime_assert( valid_tree );
		f.reorder( jumps(1,1) );
		tt << "oldtree: " << pose.fold_tree() << std::endl << "newtree: " << f << std::endl;

		pose.fold_tree( f );
	}


	tt << "start mapping: " << std::endl;
	mapping.show();


	// alternate approach:

	// first delete all unaligned residues
	while ( !mapping.all_aligned() ) {
		for ( Size i=1; i<= mapping.size1(); ++i ) {
			if ( !mapping[i] ) {
				pose.conformation().delete_residue_slow( i );
				//pose.conformation().delete_polymer_residue( i );
				mapping.delete_source_residue( i );
				break; // because the numbering is screwed up now, start the loop again
			}
		}
	}

	// now convert sequence of aligned positions
	ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
	{
		for ( Size i=1; i<= mapping.size1(); ++i ) {
			char const new_seq( target_seq[ mapping[i]-1 ] ); // strings are 0-indexed
			if ( new_seq != pose.residue(i).name1() ) {
				// will fail if get_representative_type can't find one
				ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name1( new_seq, pose.residue(i).type().variant_types() ) );
				ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type, pose.residue(i), pose.conformation() ) );
				pose.replace_residue( i, *new_rsd, false );
			}
		}
	}

	// Add jumps+cuts at the loop positions
	{
		runtime_assert( pose.total_residue() == mapping.size1() );
		kinematics::FoldTree f( pose.fold_tree() );
		for ( Size i=1; i< pose.total_residue(); ++i ) {
			runtime_assert( mapping[i] );
			if ( mapping[i+1] != mapping[i]+1 && !f.is_cutpoint(i) ) {
				f.new_jump( i, i+1, i );
				//    } else if ( pose.chain(i) != pose.chain(i+1) ) {
				//     tt << "chain jump! " << pose.residue(i).name() << ' ' << pose.residue(i+1).name() << std::endl;
				//     f.new_jump( i, i+1, i );
			}
		}
		pose.fold_tree(f);
		tt << "FOldtree: " << f << std::endl;
		runtime_assert( f.check_fold_tree() );
	}

	// add terminal residues
	// nterm

	while ( mapping[ 1 ] != 1 ) {
		int const aligned_pos( mapping[1] - 1 );
		char const new_seq( target_seq[ aligned_pos-1 ] ); // 0-indexed
		// The representative type should have no/minimal variants added.
		ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name1( new_seq ) );
		ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
		pose.conformation().safely_prepend_polymer_residue_before_seqpos( *new_rsd, 1, true );
		pose.set_omega( 1, 180.0 );
		mapping.insert_aligned_residue( 1, aligned_pos );
	}
	// cterm
	while ( mapping[ mapping.size1() ] != mapping.size2() ) {
		int const seqpos( mapping.size1() + 1 );
		int const aligned_pos( mapping[seqpos-1] + 1 );
		char const new_seq( target_seq[ aligned_pos-1 ] ); // 0-indexed
		// The representative type should have no/minimal variants added.
		ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name1( new_seq ) );
		ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
		pose.conformation().safely_append_polymer_residue_after_seqpos( *new_rsd, seqpos-1, true );
		pose.set_omega( seqpos-1, 180.0 );
		mapping.insert_aligned_residue( seqpos, aligned_pos );
	}

	// now fill in the breaks
	{
		for ( Size cut=1; cut<= Size(pose.fold_tree().num_cutpoint()); ++cut ) {
			while ( mapping[ pose.fold_tree().cutpoint( cut )+1] != Size(pose.fold_tree().cutpoint( cut )+1 ) ) {
				Size const cutpoint( pose.fold_tree().cutpoint( cut ) );
				runtime_assert( mapping[cutpoint] == cutpoint ); // we've fixed everything up til here

				if ( pose.chain( cutpoint ) != pose.chain( cutpoint+1 ) &&
						pose.residue( cutpoint ).is_protein() && pose.residue( cutpoint+1 ).is_protein() ) {
					utility_exit_with_message( "dont know whether to add residues before or after the cutpoint between chains!" );
				}
				int const aligned_pos( cutpoint+1 );
				char const new_seq( target_seq[ aligned_pos - 1 ] ); // 0-indexed
				// The representative type should have no/minimal variants added.
				ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name1( new_seq ) );
				ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
				mapping.insert_aligned_residue( cutpoint + 1, aligned_pos );
				if ( pose.residue( cutpoint ).is_protein() ) {
					// add at the nterm of the loop
					pose.conformation().safely_append_polymer_residue_after_seqpos( *new_rsd, cutpoint, true );
					pose.set_omega( cutpoint, 180.0 );
				} else if ( pose.residue( cutpoint + 1 ).is_protein() ) {
					/// add at the cterm of the loop
					pose.conformation().safely_prepend_polymer_residue_before_seqpos( *new_rsd, cutpoint+1, true );
					pose.set_omega( cutpoint+1, 180.0 );
				} else {
					utility_exit_with_message( "dont know how to add non-protein residues!" );
				}
				tt.Trace << "added residue " << cutpoint+1 << ' ' << pose.fold_tree();
				runtime_assert( pose.fold_tree().check_fold_tree() );
			} // while missing residues

			// add cutpoint variants to allow loop modeling
			int const cutpoint( pose.fold_tree().cutpoint( cut ) );
			if ( pose.chain(cutpoint) == pose.chain(cutpoint+1) ) {
				if ( pose.residue( cutpoint ).is_protein() && pose.residue( cutpoint+1 ).is_protein() ) {
					core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, cutpoint   );
					core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, cutpoint+1 );
				} else if ( pose.residue( cutpoint ).is_protein() || pose.residue( cutpoint+1 ).is_protein() ) {
					tt.Warning << "Same-chain junction between protein and non-protein!" << std::endl;
				}
			}
		} // cut=1,ncutpoints
	} // scope

	runtime_assert( mapping.is_identity() );
	runtime_assert( pose.sequence() == target_seq );

}

/// @details  Given a sequence mapping which may have simple indels, trim back around those indels so that the
/// loops can plausibly be closed.

void
trim_back_sequence_mapping(
	//pose::Pose const & src_pose,
	id::SequenceMapping & mapping,
	std::string const & source_seq,
	std::string const & target_seq,
	Size const min_loop_size
)
{
	//Size const min_loop_size( 5 );

	std::string s1( source_seq ), s2( target_seq );

	for ( Size r=1; r<= 2; ++r ) {
		runtime_assert( s1.size() == mapping.size1() && s2.size() == mapping.size2() );

		// look for insertions that are going to be hard to close
		for ( Size i=2; i<= mapping.size1(); ++i ) {
			if ( mapping[i-1] && !mapping[i] ) {
				// i is the 1st unaligned residue
				// j is the next aligned residue
				Size j(i+1);
				for ( ; j<= mapping.size1() && !mapping[j]; ++j ) {};
				//if ( j > mapping.size1() ) break; // done scanning
				if ( j > mapping.size1() ) {
					// terminal loop
					runtime_assert( j == mapping.size1()+1 );
					Size loop_begin( i );
					Size loop_size( mapping.size1() - i + 1 );
					while ( mapping[ loop_begin-2 ] && loop_size < min_loop_size ) {
						runtime_assert( !mapping[ loop_begin   ] );
						runtime_assert(  mapping[ loop_begin-1 ] );
						--loop_begin;
						mapping[ loop_begin ] = 0;
						++loop_size;
						tt.Trace << "trimming back c-terminal loop! from " << i << " to " << loop_begin << std::endl;
					}
					break;
				}

				Size loop_begin( i ), loop_end( j-1 );
				while ( true ) {
					runtime_assert( !mapping[ loop_begin   ] && !mapping[ loop_end     ] );
					runtime_assert(  mapping[ loop_begin-1 ] &&  mapping[ loop_end + 1 ] );
					Size const size1( loop_end - loop_begin + 1 );
					Size const size2( mapping[ loop_end+1 ] - mapping[ loop_begin-1 ] - 1 );
					Size const min_size( std::min( size1, size2 ) );
					//Real const size_difference( std::abs( int(size1) - int(size2) ) );
					tt.Trace << "loopseq: " << s1.substr(loop_begin-1,size1) << std::endl;
					if ( min_size >= min_loop_size /* && size_difference / min_size <= 1.001 */ ) break;
					// trim back on one side or the other
					bool const safe_to_trim_backward( loop_begin >                  2 && mapping[ loop_begin-2 ] );
					bool const safe_to_trim_forward ( loop_end   <= mapping.size1()-2 && mapping[ loop_end+2   ] );
					if ( safe_to_trim_forward && ( numeric::random::uniform() < 0.5 || !safe_to_trim_backward ) ) {
						tt.Trace << "trimming forward " << loop_begin << ' ' << loop_end << ' ' << std::endl;
						++loop_end;
						mapping[ loop_end ] = 0;
					} else if ( safe_to_trim_backward ) {
						tt.Trace << "trimming backward " << loop_begin << ' ' << loop_end << ' ' << std::endl;
						--loop_begin;
						mapping[ loop_begin ] = 0;
					} else {
						tt.Warning << "Unable to extend loop to meet criteria! " << loop_begin << ' ' << loop_end << ' '<<
							mapping[ loop_begin-1 ] << ' ' << mapping[ loop_end+1 ] << std::endl;
						break;
					}
				}
				tt.Trace << "finalloopseq: " << s1.substr(loop_begin-1,loop_end-loop_begin+1) << ' ' << loop_begin << ' ' <<
					loop_end;
				if ( r==1 ) tt.Trace << " source_insert" << std::endl;
				else tt.Trace << " target_insert" << std::endl;
			}
		}
		mapping.reverse();
		std::string tmp( s2 );
		s2 = s1;
		s1 = tmp;
	}
}


//////////////////////////////////////////////////////////////////
/// @details read in secondary structure definition from a psipred_ss2 file and
/// set that to a Pose. A temporary setup since there is very limited ss support in
/// Pose right now.
///////////////////////////////////////////////////////////////////
bool
set_secstruct_from_psipred_ss2(
	pose::Pose & pose
)

{
	utility::vector1< char > secstructs = read_psipred_ss2_file( pose );
	// if don't have a psipred file...
	if ( secstructs.size() == 0 ) {
		// instead roughly guess at secondary structure.
		core::pose::set_ss_from_phipsi( pose );
		return false;
	}
	tt << "set_secstruct for pose: ";
	for ( Size i = 1, iend=secstructs.size(); i <= iend; ++i ) {
		pose.set_secstruct( i, secstructs[i] );
		tt << pose.secstruct(i);
	}
	tt << std::endl;

	tt << "set pose secstruct from psipred_ss2 file succesfully " << std::endl;

	return true;
}


//////////////////////////////////////////////////////////////////
/// @details read in secondary structure definition from a DSSP file and
/// set that to a Pose. (Rhiju's copy of Chu's psipred reader).
/// Chu, should these functions be moved somewhere to core, along with
/// pdb readers, etc?
///////////////////////////////////////////////////////////////////
bool
set_secstruct_from_dssp(
	pose::Pose & pose,
	std::string const & filename
)
{

	utility::io::izstream data( filename );
	if ( !data ) {
		tt << "can not open DSSP file " << filename << std::endl;
		return false;
	}

	utility::vector1< char > secstructs;
	std::string line;

	//Sorry, this is the least robust DSSP reader ever.

	//Skip ahead to the good stuff
	while ( getline( data, line ) ) {
		std::istringstream line_stream( line );
		char dummy;
		line_stream >> dummy;
		if ( dummy == '#' )  break; //yippeeee!
	}

	Size count(0);
	while ( getline( data, line ) ) {
		std::istringstream line_stream( line );
		//int dummy_int;
		char aa,sec;

		line_stream.ignore( 13 );
		line_stream.get( aa );
		line_stream.ignore( 2 );
		line_stream.get( sec );

		if ( aa == '!' ) continue; //skip a chainbreak line

		count++;
		if ( aa != oneletter_code_from_aa( pose.aa(count) ) ) {
			//std::cerr << " sequence mismatch between pose and dssp " << oneletter_code_from_aa( pose.aa(count) )
			//     << " vs " << aa << " at seqpos " << count << std::endl;
			tt.Error << " sequence mismatch between pose and dssp " << oneletter_code_from_aa( pose.aa(count) )
				<< " vs " << aa << " at seqpos " << count << std::endl;
			return false;
		}

		// Follow convention used in Rosetta++, for now.
		// Do we want to change this? E.g., why is 'G' (left-handed helix) in the
		// same category as 'H' (alpha helix)? Need to discuss
		// with Ben Blum.
		//
		if ( sec == 'E' ) {
			secstructs.push_back( 'E' );
		} else if ( sec == 'H' || sec == 'I' || sec == 'G' ) {
			secstructs.push_back( 'H' );
		} else {
			secstructs.push_back( 'L' ); //all other B, S, and T
		}
	}

	runtime_assert( secstructs.size() == pose.total_residue() );
	tt.Trace << "set_secstruct for pose: ";
	for ( Size i = 1, iend=secstructs.size(); i <= iend; ++i ) {
		pose.set_secstruct( i, secstructs[i] );
		tt.Trace << pose.secstruct(i);
	}
	tt.Trace << std::endl;

	tt.Trace << "set pose secstruct from DSSP file succesfully " << std::endl;

	return true;

}

//////////////////////////////////////////////////////////////////////////////////
/// @details   set ideal BB geometry; this must occur so that loops with missing density work.
void idealize_loop(
	core::pose::Pose & pose,
	Loop const & loop
)
{
	for ( Size i = (Size)std::max((int)1,(int)loop.start()); i <= loop.stop(); ++i ) {
		core::conformation::idealize_position(i, pose.conformation());
	}
}


//////////////////////////////////////////////////////////////////////////////////
/// @details  Set a loop to extended torsion angles.
void set_extended_torsions(
	core::pose::Pose & pose,
	Loop const & loop
)
{
	tt.Error << "USING A DEPRECATED FUNCTION!( loops::set_extended_torsions(...) ) " << std::endl;

	Real const init_phi  ( -150.0 );
	Real const init_psi  (  150.0 );
	Real const init_omega(  180.0 );

	static int counter = 0;

	tt.Debug << "Extending loop " << loop.start() << " " << loop.stop() << std::endl;

	idealize_loop(pose, loop );

	Size start_extended = std::max((Size)1,loop.start());
	Size end_extended   = std::min(pose.total_residue(),loop.stop());
	for ( Size i = start_extended; i <= end_extended; ++i ) {
		if ( i != start_extended ) pose.set_phi( i, init_phi );
		if ( i != end_extended ) pose.set_psi( i, init_psi );
		if ( ( i != start_extended ) && ( i != end_extended ) ) pose.set_omega( i, init_omega );
	}

	counter++;
}


//////////////////////////////////////////////////////////////////////////////////
/// @details  Rebuild a loop via fragment insertion + ccd closure + minimization
void remove_missing_density(
	core::pose::Pose & pose,
	Loop const & loop
)
{
	Size cutpoint = loop.cut();
	if ( cutpoint== 0 ) cutpoint =  (loop.start() + loop.stop()) / 2;
	set_single_loop_fold_tree( pose, loop );
	set_extended_torsions( pose, loop );
}


core::Real native_loop_core_CA_rmsd(
	const core::pose::Pose & native_pose,
	const core::pose::Pose & pose,
	loops::Loops loops,
	int &corelength
) {

	std::vector< core::Size > residue_exclusion;

	for ( core::Size i = 1; i <= loops.size(); i++ ) {
		for ( core::Size ir =  loops[i].start();
				ir <= loops[i].stop();
				ir ++ ) {
			residue_exclusion.push_back( ir );
		}
	}


	std::list< core::Size > residue_selection;

	for ( core::Size ir = 1; ir <= pose.total_residue(); ir ++ ) {
		if ( !pose.residue_type(ir).is_protein() ) continue;
		bool exclude = false;
		for ( core::Size p = 0; p < residue_exclusion.size(); p++ ) {
			if ( ir == residue_exclusion[p] ) {
				exclude = true;
				break;
			}
		}
		if ( !exclude ) residue_selection.push_back( ir );
	}

	corelength = residue_selection.size();
	if ( corelength == 0 ) return 0.0;
	else return core::scoring::CA_rmsd( native_pose, pose, residue_selection );
} // native_CA_rmsd


/////////////////////////////////////////////////////////////////////////////////////
/// @details pose1 is the reference and pose2 is the structure for which rmsd is calculated.
///The rmsd is calculated over four backbone atoms of all loop residues, assuming template
///regions in pose1 and pose2 are already aligned.
//////////////////////////////////////////////////////////////////////////////////////
Real
loop_rmsd_with_superimpose(
	pose::Pose const & pose1,
	pose::Pose const & pose2,
	Loops const & loops,
	bool CA_only,
	bool bb_only
)

{
	Loops core; //empty core -- take all residues for superposition
	return loop_rmsd_with_superimpose_core( pose1, pose2, loops, core, CA_only, bb_only );
}

/////////////////////////////////////////////////////////////////////////////////////
/// @details pose1 is the reference and pose2 is the structure for which rmsd is calculated.
///The rmsd is calculated over four backbone atoms of all loop residues, assuming template
///regions in pose1 and pose2 are already aligned.
//////////////////////////////////////////////////////////////////////////////////////
Real
loop_rmsd_with_superimpose_core(
	pose::Pose const & pose1,
	pose::Pose const & pose2,
	Loops const & loops,
	Loops const & core,
	bool CA_only,
	bool bb_only
)

{
	pose::Pose ncpose1 = pose1;
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, ncpose1, core::id::BOGUS_ATOM_ID );
	for ( core::Size ir=1; ir <= ncpose1.total_residue(); ++ir ) {
		if ( ncpose1.residue(ir).aa() == core::chemical::aa_vrt || pose2.residue(ir).aa() == core::chemical::aa_vrt ) continue;
		if ( !loops.is_loop_residue( ir ) && ( core.size()==0 || core.is_loop_residue( ir ) ) ) {
			id::AtomID const id1( ncpose1.residue(ir).atom_index("CA"), ir );
			id::AtomID const id2( pose2.residue(ir).atom_index("CA"), ir );
			atom_map.set(id1, id2);
		}
	}
	core::scoring::superimpose_pose( ncpose1, pose2, atom_map );

	core::Real looprms = loop_rmsd( ncpose1, pose2, loops, CA_only, bb_only );

	return looprms;
}

/////////////////////////////////////////////////////////////////////////////////////
/// @details pose1 is the reference and pose2 is the structure for which rmsd is calculated.
///The rmsd is calculated over four backbone atoms of all loop residues, assuming template
///regions in pose1 and pose2 are already aligned.
//////////////////////////////////////////////////////////////////////////////////////
Real
loop_rmsd(
	pose::Pose const & pose1,
	pose::Pose const & pose2,
	Loops const & loops,
	bool CA_only,
	bool bb_only
)

{
	if ( pose1.total_residue() != pose2.total_residue() ) {
		// strip VRTs from the end, then compare lengths
		int nnonvrt_1 = pose1.total_residue();
		while ( pose1.residue( nnonvrt_1 ).aa() == core::chemical::aa_vrt ) nnonvrt_1--;

		int nnonvrt_2 = pose2.total_residue();
		while ( pose2.residue( nnonvrt_2 ).aa() == core::chemical::aa_vrt ) nnonvrt_2--;
		if ( nnonvrt_1 != nnonvrt_2 ) {
			utility_exit_with_message( "Error in loop_rmsd: pose1.total_residue() != pose2.total_residue() ( "
				+ string_of( nnonvrt_1 ) + "!=" + string_of( nnonvrt_2 )+")" );
		}
	}

	Real rms = 0.0;
	int atom_count(0);
	for ( Loops::const_iterator it=loops.begin(), it_end=loops.end();
			it != it_end; ++it ) {
		for ( Size i = it->start(); i<=it->stop(); ++i ) {
			if ( i > pose1.total_residue() ) {
				tt.Warning <<  "[Warning]: Pose1: Loop residue " << i << "exceeds pose1 size " << pose1.total_residue() << std::endl;
				continue;
			}
			if ( i > pose2.total_residue() ) {
				tt.Warning <<  "[Warning]: Pose1: Loop residue " << i << "exceeds pose2 size " << pose2.total_residue() << std::endl;
				continue;
			}
			Size start = 1;
			Size end = bb_only ? 4 : pose1.residue(i).atoms().size();
			if ( CA_only ) { // Only include CA atoms
				start = 2;
				end = 2;
			}
			for ( Size j = start; j <= end; ++j ) {
				atom_count++;
				if ( j >  pose1.residue(i).atoms().size()  ) {
					tt.Warning <<  "[Warning]: Pose1: Atom " << j << " missing from residue " << i << std::endl;
					continue;
				}
				if ( j > pose2.residue(i).atoms().size()  ) {
					tt.Warning <<  "[Warning]: Pose2: Atom " << j << " missing from residue " << i << std::endl;
					continue;
				}
				rms += pose1.residue(i).xyz(j).distance_squared( pose2.residue(i).xyz(j) );

			} // for each bb atom
		} // for each residue
	} // for each loop
	return atom_count==0 ? 0.0 : std::sqrt(rms/atom_count);
}

/////////////////////////////////////////////////////////////////////////////////////
/// @details pose1 is the reference and pose2 is the structure for which rmsd is calculated.
///The rmsd is calculated over four backbone atoms of each loop after fitting it onto the reference
///loop and when there are multiple loops, return the mean value.
//////////////////////////////////////////////////////////////////////////////////////
Real
loop_local_rmsd(
	pose::Pose const & pose1,
	pose::Pose const & pose2,
	Loops const & loops
)
{
	runtime_assert(pose1.total_residue() == pose2.total_residue() );

	Real rms = 0.0;
	if ( loops.num_loop() == 0 ) {
		return rms;
	}

	for ( Loops::const_iterator it=loops.begin(), it_end=loops.end();
			it != it_end; ++it ) {
		int natoms = 4 * it->size() ;
		//FArray2D_double p1a(3, natoms);
		//FArray2D_double p2a(3, natoms);
		FArray2D< core::Real > p1a(3, natoms);
		FArray2D< core::Real > p2a(3, natoms);

		int atom_count(0);
		for ( Size i = it->start(); i<=it->stop(); ++i ) {
			for ( Size j = 1; j <= 4; ++j ) {
				const numeric::xyzVector< Real > & vec1( pose1.residue(i).xyz( j ) );
				const numeric::xyzVector< Real > & vec2( pose2.residue(i).xyz( j ) );
				atom_count++;
				for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
					p1a(k+1,atom_count) = vec1[k];
					p2a(k+1,atom_count) = vec2[k];
				} // k
			} // j
		} // i
		rms += numeric::model_quality::rms_wrapper( natoms, p1a, p2a );
	}

	return rms / loops.num_loop();
}

} // namespace loops
} // namespace protocols
