// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA de novo fragment assembly Structure Parameters
/// @brief User input parameters for RNA structure inference
/// @details
/// @author Rhiju Das


// Unit Headers
#include <protocols/farna/setup/RNA_DeNovoPoseInitializer.hh>
#include <protocols/farna/secstruct/RNA_SecStructLegacyInfo.hh>
#include <protocols/farna/movers/RNA_JumpMover.hh>
#include <protocols/farna/libraries/RNA_LibraryManager.hh>
#include <protocols/farna/setup/RNA_DeNovoParameters.hh>
#include <protocols/farna/base_pairs/RNA_BasePairHandler.hh>
#include <protocols/farna/util.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/rna/util.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/scoring/ScoreType.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/tools/make_vector1.hh>

#include <numeric/random/random.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <list>

#include <utility/vector1.hh>
#include <utility/stream_util.hh>

using namespace ObjexxFCL; // AUTO USING NS
using namespace ObjexxFCL::format; // AUTO USING NS
using namespace protocols::farna::secstruct;

namespace protocols {
namespace farna {
namespace setup {

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// For rna_denovo (FARNA/FARFAR) runs that start from scratch, this object handles basic pose setup.
//
// This used to be called RNA_DeNovoPoseInitializer, which was an object that had accumulated too many
//  functionalities.
//
// In the near future, this format will be deprecated in favor of direct command-line input of sequence &
//  secondary structure -- most of the code is worked out in stepwise monte carlo on RNA/proteins, but needs
//  to be ported over here. -- rhiju, 2014
//////////////////////////////////////////////////////////////////////////////////////////////////////////

static THREAD_LOCAL basic::Tracer TR( "protocols.farna.setup.RNA_DeNovoPoseInitializer" );

using namespace core;

RNA_DeNovoPoseInitializer::RNA_DeNovoPoseInitializer( RNA_DeNovoParameters const & rna_params ):
	rna_params_( rna_params ), // note that this object might be updated later...
	assume_non_stem_is_loop( false ),
	bps_moves_( false ),
	root_at_first_rigid_body_( false )
{
}

RNA_DeNovoPoseInitializer::~RNA_DeNovoPoseInitializer() {}

//////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::initialize_for_de_novo_protocol(
	core::pose::Pose & pose,
	bool const ignore_secstruct /* = false */ )
{
	initialize_secstruct( pose );
	if  ( ignore_secstruct ) override_secstruct( pose );
	if ( rna_params_.virtual_anchor_attachment_points_.size() > 0 ) append_virtual_anchor( pose );
	setup_virtual_phosphate_variants( pose );
}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::override_secstruct( core::pose::Pose & pose ){
	rna_params_.rna_secstruct_legacy_ = std::string( pose.total_residue(), 'X' );
	TR << "OVER-RIDING SECONDARY STRUCTURE WITH:   " << rna_params_.rna_secstruct_legacy_ << std::endl;
	set_rna_secstruct_legacy( pose, rna_params_.rna_secstruct_legacy_ );
}
/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::append_virtual_anchor( pose::Pose & pose )
{

	if ( rna_params_.virtual_anchor_attachment_points_.size() == 0 ) return;

	TR.Debug << "Current last residue is type: " << pose.residue( pose.total_residue() ).name3()  << std::endl;
	TR.Debug << pose.annotated_sequence() << std::endl;
	if ( pose.residue( pose.total_residue() ).name3() == "XXX" ) return; //already did virtual residue attachment.

	// Fix up the pose.
	core::chemical::ResidueTypeSetCOP residue_set = pose.residue_type(1).residue_type_set();

	// std::cout << " CHECK XXX " << residue_set.name3_map("XXX").size() << std::endl;
	// std::cout << " CHECK YYY " << residue_set.name3_map("YYY").size() << std::endl;
	// std::cout << " CHECK VRT " << residue_set.name3_map("VRT").size() << std::endl;

	core::chemical::ResidueTypeCOP rsd_type( residue_set->get_representative_type_name3("XXX") );
	core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *rsd_type ) );
	pose.append_residue_by_jump( *new_res, rna_params_.virtual_anchor_attachment_points_[1] );

	Size const virt_res = pose.total_residue();

	// Info on pairings and cutpoints.
	rna_params_.cutpoints_open_.push_back( (virt_res - 1) );

	for ( Size n = 1; n <= rna_params_.virtual_anchor_attachment_points_.size(); n++ ) {

		core::pose::rna::BasePair p;
		p.set_res1( rna_params_.virtual_anchor_attachment_points_[ n ] );
		p.set_res2( virt_res );
		rna_params_.rna_pairing_list_.push_back( p );

		utility::vector1< Size > obligate_pair;
		obligate_pair.push_back( rna_params_.rna_pairing_list_.size() );
		rna_params_.obligate_pairing_sets_.push_back( obligate_pair );
	}


}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::initialize_secstruct( core::pose::Pose & pose  )
{
	std::string & rna_secstruct_legacy = rna_params_.rna_secstruct_legacy_;

	if ( !rna_params_.secstruct_defined_ ) {

		rna_secstruct_legacy = std::string( pose.total_residue(), 'X' );

		if ( rna_params_.rna_pairing_list_.size() > 0 && assume_non_stem_is_loop ) {
			rna_secstruct_legacy = std::string( pose.total_residue(), 'L' );
		}

		for ( Size n = 1; n <= rna_params_.rna_pairing_list_.size(); n++ ) {
			core::pose::rna::BasePair const & rna_pairing( rna_params_.rna_pairing_list_[ n ] );
			if (  rna_pairing.edge1() == WATSON_CRICK &&
					rna_pairing.edge2() == WATSON_CRICK &&
					rna_pairing.orientation() == ANTIPARALLEL &&
					core::chemical::rna::possibly_canonical( pose.residue( rna_pairing.res1() ).aa(),
					pose.residue( rna_pairing.res2() ).aa() ) )  {
				rna_secstruct_legacy[ rna_pairing.res1() - 1 ] = 'H';
				rna_secstruct_legacy[ rna_pairing.res2() - 1 ] = 'H';
			}
		}
	}

	TR << "Setting desired secondary structure to: " << rna_secstruct_legacy << std::endl;

	set_rna_secstruct_legacy( pose, rna_secstruct_legacy );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::insert_base_pair_jumps( pose::Pose & pose, RNA_JumpMover const & jump_mover,  bool & success ) const
{

	Size const num_jump( pose.num_jump() );

	for ( Size i = 1; i <= num_jump; i++ ) {

		// Check that we can actually insert here. At least one of the jump partners
		// should allow moves. (I guess the other one can stay fixed).
		Size const jump_pos1( pose.fold_tree().upstream_jump_residue( i ) );
		Size const jump_pos2( pose.fold_tree().downstream_jump_residue( i ) );

		if ( moveable_jump( jump_pos1, jump_pos2, *(jump_mover.atom_level_domain_map()) ) ) jump_mover.add_new_RNA_jump( pose, i, success );

		if ( !success ) return;

	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////
// In principle, don't need to pass atom_level_domain_map into this function (its stored in rna_jump_mover),
// but I want to make sure that user known that  atom_level_domain_map will get updated -- rhiju, 2015
void
RNA_DeNovoPoseInitializer::setup_fold_tree_and_jumps_and_variants( pose::Pose & pose,
	RNA_JumpMover const & rna_jump_mover,
	protocols::toolbox::AtomLevelDomainMapOP atom_level_domain_map ) const
{
	runtime_assert( rna_jump_mover.atom_level_domain_map() == atom_level_domain_map );
	setup_jumps( pose, rna_jump_mover );
	setup_chainbreak_variants( pose, atom_level_domain_map );
	runtime_assert( rna_jump_mover.atom_level_domain_map() == atom_level_domain_map );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Only in use by cs_rosetta_rna.
// Creates a temporary JumpLibrary and AtomLevelDomainMap to do the jump setup.
void
RNA_DeNovoPoseInitializer::setup_fold_tree_and_jumps_and_variants( pose::Pose & pose ) const
{
	using namespace protocols::toolbox;
	AtomLevelDomainMapOP atom_level_domain_map( new AtomLevelDomainMap( pose ) );
	RNA_JumpMover const rna_jump_mover( RNA_LibraryManager::get_instance()->rna_jump_library_cop(), atom_level_domain_map );
	setup_fold_tree_and_jumps_and_variants( pose, rna_jump_mover, atom_level_domain_map );
}

///////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::setup_jumps( pose::Pose & pose, RNA_JumpMover const & rna_jump_mover ) const
{
	using namespace core::pose::rna;

	///////////////////////////////////////////////////////////
	// Basic setup ==> How many jumps? cuts?
	///////////////////////////////////////////////////////////
	Size const nres = pose.total_residue();
	kinematics::FoldTree f( nres );

	Size const num_cuts_closed( rna_params_.cutpoints_closed_.size() );
	Size const num_cuts_open  ( rna_params_.cutpoints_open_.size() );
	Size const num_cuts_total ( num_cuts_closed + num_cuts_open );

	////////////////////////////////////////////////////////////////////////
	utility::vector1< utility::vector1< Size > > obligate_pairing_sets,  stem_pairing_sets;
	obligate_pairing_sets = rna_params_.obligate_pairing_sets();
	if ( bps_moves_ ) { // supplement obligate_pairing_sets with stems in freely moving regions.
		for ( Size n = 1; n <= rna_params_.stem_pairing_sets_.size(); n++ ) {
			for ( Size m = 1; m <= rna_params_.stem_pairing_sets_[n].size(); m++ ) {
				Size const idx = rna_params_.stem_pairing_sets_[n][m];
				if ( rna_params_.check_in_pairing_sets( obligate_pairing_sets, rna_params_.rna_pairing_list_[ idx ] ) ) continue;
				obligate_pairing_sets.push_back( utility::tools::make_vector1(idx ) );
			}
		}
	} else {
		stem_pairing_sets = rna_params_.stem_pairing_sets_;
	}

	Size const num_pairings( rna_params_.rna_pairing_list_.size() );
	Size const num_obligate_pairing_sets( obligate_pairing_sets.size() );
	Size const num_stem_pairing_sets( stem_pairing_sets.size() );
	Size const num_chain_connections( rna_params_.chain_connections_.size() );
	// std::cout << num_stem_pairing_sets << " + " <<  num_obligate_pairing_sets << " <= " << num_pairings << std::endl;
	runtime_assert( num_stem_pairing_sets + num_obligate_pairing_sets <= num_pairings );

	//////////////////////////////////////////////////////////////////////
	// Cuts.
	//////////////////////////////////////////////////////////////////////
	utility::vector1< Size > obligate_cut_points;
	for ( Size n = 1; n<= num_cuts_closed; n++ )   obligate_cut_points.push_back( rna_params_.cutpoints_closed_[ n ] );
	for ( Size n = 1; n<= num_cuts_open  ; n++ )   obligate_cut_points.push_back( rna_params_.cutpoints_open_[n] );

	//////////////////////////////////////////////////////////////////////
	// base pair steps are a special kind of chunk, created "on-the-fly"
	// from a database. These stem base pairs will be obligate pairs
	//  (see above get_pairings_from_line ), and we can define cutpoints
	//  ahead of time.
	//////////////////////////////////////////////////////////////////////
	utility::vector1< Size > base_pair_step_starts;
	RNA_BasePairHandler const rna_base_pair_handler( rna_params_ ); // has get_base_pair_steps() function.
	if ( bps_moves_ ) {
		utility::vector1< BasePairStep > base_pair_steps = rna_base_pair_handler.get_base_pair_steps( false /* just canonical */ );

		for ( Size n = 1; n <= base_pair_steps.size(); n++ ) {
			BasePairStep const & base_pair_step = base_pair_steps[n];
			runtime_assert( rna_params_.check_in_pairing_sets( obligate_pairing_sets, BasePair( base_pair_step.i(),      base_pair_step.j_next() ) ) );
			runtime_assert( rna_params_.check_in_pairing_sets( obligate_pairing_sets, BasePair( base_pair_step.i_next(), base_pair_step.j()      ) ) );
			runtime_assert ( base_pair_step.i_next() == base_pair_step.i()+1 );
			if ( base_pair_step.j_next() != (base_pair_step.j() + 1) ) continue;

			// used below in cutpoint setting...
			base_pair_step_starts.push_back( base_pair_step.i() );
			base_pair_step_starts.push_back( base_pair_step.j() );

			// choose one side of the base pair step to place a cutpoint.
			if ( obligate_cut_points.has_value( base_pair_step.i() ) ) continue;
			if ( obligate_cut_points.has_value( base_pair_step.j() ) ) continue;
			// flip a coin
			if ( numeric::random::rg().random_range( 0, 1 ) ) {
				obligate_cut_points.push_back( base_pair_step.i() );
			} else {
				obligate_cut_points.push_back( base_pair_step.j() );
			}
		}
	}

	////////////////////////////////////////////////////////////////////////
	// Two possibilities for desired fold tree topology:
	//   Jumps dominate, or cuts dominate.
	Size const num_pairings_to_force = std::max( num_obligate_pairing_sets + num_chain_connections,
		num_cuts_total );


	////////////////////////////////////////////////////////////////////////
	ObjexxFCL::FArray2D <int> jump_points( 2, num_pairings_to_force );
	ObjexxFCL::FArray1D <int> cuts( num_pairings_to_force );

	//////////////////////////////////////////////////////////////////////
	// If a cut needs to be randomly chosen, will generally try to
	// place it in a loopy region.
	FArray1D_float cut_bias( nres, 1.0 );
	std::string const & rna_secstruct_legacy( get_rna_secstruct_legacy( pose ) );
	for ( Size i = 1; i < nres; i++ ) {
		// no cuts outside RNA
		if ( !pose.residue(i+1).is_RNA() ) {
			cut_bias(i) = 0.0;
			continue;
		}
		// Reduced probability of cuts inside helices
		if ( rna_secstruct_legacy[i] == 'H' && rna_secstruct_legacy[i+1] == 'H' ) {
			cut_bias( i )   = 0.1;
		}
		// No cuts inside user_input domains.
		Size domain( rna_jump_mover.atom_level_domain_map()->get_domain( core::id::AtomID( named_atom_id_to_atom_id( core::id::NamedAtomID( " P  ", i+1 ), pose )  ) ) );
		if ( domain > 0 ) {
			if ( bps_moves_ && base_pair_step_starts.has_value( i ) ) {
				cut_bias( i ) = 0.01;
			} else {
				cut_bias( i ) = 0.0;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////
	// Jump residues.
	//////////////////////////////////////////////////////////////////////
	Size ntries( 0 );
	Size const MAX_TRIES( 1000 );
	bool success( false );

	while ( !success && ++ntries < MAX_TRIES ) {

		// First, put in a jump from each obligate pairing set.
		Size count( 0 );

		for ( Size n = 1; n <= num_obligate_pairing_sets ; n++ ) {
			Size const pairing_index_in_stem( static_cast<Size>( numeric::random::rg().uniform() * obligate_pairing_sets[n].size() )  + 1 );
			Size const which_pairing = obligate_pairing_sets[n][pairing_index_in_stem];
			count++;
			jump_points(1, count) = rna_params_.rna_pairing_list_[which_pairing].res1();
			jump_points(2, count) = rna_params_.rna_pairing_list_[which_pairing].res2();
			//TR << "JUMPS1 " <<  jump_points(1,count) << ' ' << jump_points(2,count ) << std::endl;
		}

		// "Chain connections" provide less information about specific residues to pair --
		//  but they're assumed to be obligatory.
		for ( Size n = 1; n <= num_chain_connections ; n++ ) {
			utility::vector1 < Size > const & res_list1( rna_params_.chain_connections_[n].first );
			utility::vector1 < Size > const & res_list2( rna_params_.chain_connections_[n].second);
			Size jump_pos1 = numeric::random::rg().random_element( res_list1 );
			Size jump_pos2 = numeric::random::rg().random_element( res_list2 );
			count++;
			jump_points(1, count) =  std::min( jump_pos1, jump_pos2 );
			jump_points(2, count) =  std::max( jump_pos1, jump_pos2 );
			//   TR << "JUMPS2 " <<  jump_points(1,count) << ' ' << jump_points(2,count ) << std::endl;
		}
		//  TR << std::endl;

		// Then, to fill out pairings, look at remaining possible pairing sets (these
		// should typically be Watson-Crick stems, but this setup is general )
		// Note that there might be three stems defined, but we only want two --
		//  following picks a random set of two.
		//
		// Following has been ad hoc for a long time -- can be made more systematic
		//  based on fold_tree build up that has been worked out in stepwise monte carlo.
		//
		FArray1D < bool > used_set( num_stem_pairing_sets, false );
		Size num_sets_left( num_stem_pairing_sets );

		while ( count < num_pairings_to_force ) {

			// Find a random pairing among what's remaining.
			Size const set_index( static_cast<Size>( numeric::random::rg().uniform() * num_sets_left )  + 1 );
			Size m( 0 ), set_count( 0 );
			for ( m = 1; m <= num_stem_pairing_sets; m++ ) {
				if ( !used_set( m ) ) {
					set_count++;
					if ( set_count == set_index ) break;
				}
			}

			if ( m > num_stem_pairing_sets ) {
				utility_exit_with_message( "Could not find a stem_pairing. Number of stem_pairing_sets: "+I(3,num_stem_pairing_sets) );
			}

			Size const pairing_index_in_set( static_cast<Size>( numeric::random::rg().uniform() * stem_pairing_sets[m].size() )  + 1 );
			Size const which_pairing = stem_pairing_sets[m][pairing_index_in_set];

			//   std::cout  << "USING SET: " << m  << " ==> " << pairing_index_in_set << std::endl;

			count++;
			jump_points(1, count) = rna_params_.rna_pairing_list_[which_pairing].res1();
			jump_points(2, count) = rna_params_.rna_pairing_list_[which_pairing].res2();

			used_set( m ) = true;
			num_sets_left--;
		}

		////////////////////////////////////////////////////////////////////////////////
		// Do it, get the fold tree. and set up the pose.
		////////////////////////////////////////////////////////////////////////////////
		std::vector< int > obligate_cut_points_reformat;
		for ( Size q = 1; q <= obligate_cut_points.size(); q++ ) obligate_cut_points_reformat.push_back( obligate_cut_points[q] );

		// TR << TR.Cyan << "Making attempt " << ntries << std::endl;
		// TR << TR.Cyan << "obligate_cutpoints " << obligate_cut_points_reformat << std::endl;
		// for (Size n = 1; n <= num_pairings_to_force; n++ ){
		// 	TR << TR.Cyan << "JUMPS " << jump_points(1, n) <<
		//     " " <<  jump_points(2, n)  <<  std::endl;
		// }

		success = f.random_tree_from_jump_points( nres, num_pairings_to_force, jump_points, obligate_cut_points_reformat, cut_bias, 1, true /*enable 1 or NRES jumps*/ );
	}

	if ( !success )  utility_exit_with_message( "Couldn't find a freaking tree!" );

	// Hold on to torsion angles in case we need to set up chainbreak residues...
	pose::Pose pose_copy = pose;

	fill_in_default_jump_atoms( f, pose );

	if ( rna_params_.virtual_anchor_attachment_points_.size() > 0 ) {

		f.reorder( pose.total_residue() ); //reroot so that virtual residue is fixed.

		if ( root_at_first_rigid_body_ ) {
			utility::vector1< Size > rigid_body_jumps = get_rigid_body_jumps( pose );
			runtime_assert( rigid_body_jumps.size() > 0 );
			Size const anchor_rsd = pose.fold_tree().upstream_jump_residue( rigid_body_jumps[1] );
			TR << "ROOTING AT RSD" << anchor_rsd << std::endl;
			f.reorder( anchor_rsd ); //reroot so that partner of virtual residue is fixed
		}

	} else {

		// also useful -- if user has an input pdb, put the root in there, if possible.
		//  for (Size n = pose.total_residue(); n >= 1; n-- ){ // not sure why I did this backwards...
		for ( Size n = 1; n <= pose.total_residue(); n++ ) {
			if ( pose.residue(n).is_RNA() &&
					rna_jump_mover.atom_level_domain_map()->get_domain( named_atom_id_to_atom_id( id::NamedAtomID( " C1'", n ), pose ) ) == 1 /*1 means the first inputted pose*/ &&
					f.possible_root(n) ) {
				f.reorder( n );
				break;
			}
		}

	}

	// Make it so.
	pose.fold_tree( f );

	bool const random_jumps( true ); // For now this is true... perhaps should also have a more deterministic procedure.
	if ( random_jumps ) insert_base_pair_jumps( pose, rna_jump_mover, success );

	if ( !success ) {
		utility_exit_with_message( "Trouble inserting base pair jumps into pose -- check residues and edges." );
	}
}


///////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::setup_chainbreak_variants( pose::Pose & pose,
	protocols::toolbox::AtomLevelDomainMapOP atom_level_domain_map ) const
{

	pose::Pose pose_copy = pose;
	utility::vector1< Size > const & cutpoints_open( rna_params_.cutpoints_open_ );

	// Create cutpoint variants to force chainbreak score computation.
	for ( Size cutpos = 1; cutpos < pose.total_residue(); cutpos++ ) {

		if ( ! pose.fold_tree().is_cutpoint( cutpos ) ) continue;

		// Don't assign a chainbreak penalty if user said this was an "open" cutpoint.
		if ( std::find( cutpoints_open.begin(), cutpoints_open.end(), cutpos) != cutpoints_open.end() ) continue;

		core::pose::correctly_add_cutpoint_variants( pose, cutpos );

		for ( Size i = cutpos; i <= cutpos + 1; i++ ) {
			for ( Size j = 1; j <= pose.residue( i ).mainchain_torsions().size(); j++ ) {
				id::TorsionID torsion_id( i, id::BB, j );
				pose.set_torsion( torsion_id, pose_copy.torsion( torsion_id ) ) ;
			} // j
		} // i
	}

	atom_level_domain_map->renumber_after_variant_changes( pose );

}


////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::setup_virtual_phosphate_variants( pose::Pose & pose ) const {

	using namespace id;
	using namespace chemical;

	if ( pose.residue( 1 ).is_RNA() ) pose::add_variant_type_to_pose_residue( pose, VIRTUAL_PHOSPHATE, 1  );

	utility::vector1< Size > const & cutpoints_open( rna_params_.cutpoints_open_ );
	for ( Size i = 1; i <= cutpoints_open.size(); i++ ) {

		Size n = cutpoints_open[ i ];

		if ( n == pose.total_residue() ) {
			utility_exit_with_message( "Do not specify cutpoint_open at last residue of model" );
		}

		if ( pose.residue_type( n   ).has_variant_type( CUTPOINT_LOWER ) ||
				pose.residue_type( n+1 ).has_variant_type( CUTPOINT_UPPER ) ) {
			utility_exit_with_message( "conflicting cutpoint_open & cutpoint_closed" );
		}

		if ( pose.residue_type( n+1 ).is_RNA() ) {
			if ( ! pose.residue_type( n+1 ).has_variant_type( VIRTUAL_PHOSPHATE) ) {
				pose::add_variant_type_to_pose_residue( pose, VIRTUAL_PHOSPHATE, n+1  );
			}
		}

	}
}

} //setup
} //farna
} //protocols


