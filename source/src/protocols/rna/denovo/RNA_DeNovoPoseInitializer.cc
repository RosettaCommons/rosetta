// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RNA de novo fragment assembly Structure Parameters
/// @brief User input parameters for RNA structure inference
/// @details
/// @author Rhiju Das


// Unit Headers
#include <protocols/rna/denovo/RNA_DeNovoPoseInitializer.hh>
#include <core/pose/rna/secstruct_legacy/RNA_SecStructLegacyInfo.hh>
#include <core/import_pose/RNA_JumpMover.hh>
#include <core/import_pose/libraries/RNA_LibraryManager.hh>
#include <core/import_pose/RNA_DeNovoParameters.hh>
#include <core/import_pose/RNA_BasePairHandler.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.hh>
#include <core/import_pose/libraries/RNA_ChunkLibrary.hh>
#include <core/import_pose/libraries/ChunkSet.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/stepwise/monte_carlo/util.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/pose/copydofs/util.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/rna/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

#include <utility/exit.hh>
#include <utility/tools/make_vector1.hh>

#include <numeric/random/random.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <fstream>

#include <utility/vector1.hh>

#include <ObjexxFCL/FArray2D.hh> // AUTO IWYU For FArray2D, FArray2D<>::size_type, FArray...
#include <core/pose/rna/BasePairStep.hh> // AUTO IWYU For BasePairStep

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;
using namespace core::pose::rna::secstruct_legacy;
using namespace core;
using namespace core::pose::rna;
using namespace core::import_pose;
using namespace core::chemical::rna;
using namespace protocols::rna::denovo;

namespace protocols {
namespace rna {
namespace denovo {

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

static basic::Tracer TR( "protocols.rna.denovo.setup.RNA_DeNovoPoseInitializer" );

using namespace core;

RNA_DeNovoPoseInitializer::RNA_DeNovoPoseInitializer( RNA_DeNovoParameters const & rna_params ):
	rna_params_( rna_params ), // note that this object might be updated later...
	assume_non_stem_is_loop( false ),
	bps_moves_( false ),
	root_at_first_rigid_body_( false ),
	dock_each_chunk_( false ),
	dock_each_chunk_per_chain_( false ),
	center_jumps_in_single_stranded_( false ),
	new_fold_tree_initializer_( false ),
	model_with_density_( false )
{
}

RNA_DeNovoPoseInitializer::~RNA_DeNovoPoseInitializer() = default;

//////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::initialize_for_de_novo_protocol(
	core::pose::Pose & pose,
	bool const ignore_secstruct /* = false */ )
{
	initialize_secstruct( pose );
	if  ( ignore_secstruct ) override_secstruct( pose );
	if ( rna_params_.virtual_anchor_attachment_points().size() > 0 ) append_virtual_anchor( pose );
	setup_virtual_phosphate_variants( pose );
}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::override_secstruct( core::pose::Pose & pose ){
	rna_params_.set_rna_secstruct_legacy( std::string( pose.size(), 'X' ) );
	TR << "OVER-RIDING SECONDARY STRUCTURE WITH:   " << rna_params_.rna_secstruct_legacy() << std::endl;
	set_rna_secstruct_legacy( pose, rna_params_.rna_secstruct_legacy() );
}
/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::append_virtual_anchor( pose::Pose & pose )
{
	if ( rna_params_.virtual_anchor_attachment_points().size() == 0 ) return;

	TR.Debug << "Current last residue is type: " << pose.residue( pose.size() ).name3()  << std::endl;
	if ( pose.residue( pose.size() ).name3() == "XXX" ) return; //already did virtual residue attachment.

	// Fix up the pose.
	core::chemical::ResidueTypeCOP rsd_type( core::pose::virtual_type_for_pose( pose ) );
	core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *rsd_type ) );
	pose.append_residue_by_jump( *new_res, rna_params_.virtual_anchor_attachment_points()[1] );

	core::Size const virt_res = pose.size();

	// Info on pairings and cutpoints.
	rna_params_.add_cutpoint_open( virt_res - 1 );

	for ( core::Size n = 1; n <= rna_params_.virtual_anchor_attachment_points().size(); n++ ) {

		core::pose::rna::BasePair p;
		p.set_res1( rna_params_.virtual_anchor_attachment_points()[ n ] );
		p.set_res2( virt_res );
		rna_params_.add_rna_pairing( p );

		utility::vector1< core::Size > obligate_pairing_set;
		obligate_pairing_set.push_back( rna_params_.rna_pairing_list().size() );
		rna_params_.add_obligate_pairing_set( obligate_pairing_set );

	}

	if ( pose::full_model_info::full_model_info_defined( pose ) ) pose::full_model_info::append_virtual_residue_to_full_model_info( pose );
}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::initialize_secstruct( core::pose::Pose & pose  )
{
	std::string rna_secstruct_legacy = rna_params_.rna_secstruct_legacy();

	if ( !rna_params_.secstruct_defined() || rna_secstruct_legacy.size() == 0 ) {

		rna_secstruct_legacy = std::string( pose.size(), 'X' );

		if ( rna_params_.rna_pairing_list().size() > 0 && assume_non_stem_is_loop ) {
			rna_secstruct_legacy = std::string( pose.size(), 'L' );
		}

		for ( core::Size n = 1; n <= rna_params_.rna_pairing_list().size(); n++ ) {
			core::pose::rna::BasePair const & rna_pairing( rna_params_.rna_pairing_list()[ n ] );
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
	rna_params_.set_rna_secstruct_legacy( rna_secstruct_legacy );
	set_rna_secstruct_legacy( pose, rna_secstruct_legacy );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::insert_base_pair_jumps( pose::Pose & pose, RNA_JumpMover const & jump_mover,  bool & success ) const
{

	core::Size const num_jump( pose.num_jump() );

	for ( core::Size i = 1; i <= num_jump; i++ ) {

		// Check that we can actually insert here. At least one of the jump partners
		// should allow moves. (I guess the other one can stay fixed).
		core::Size const jump_pos1( pose.fold_tree().upstream_jump_residue( i ) );
		core::Size const jump_pos2( pose.fold_tree().downstream_jump_residue( i ) );

		// check that they're both RNA?
		if ( !( pose.residue( jump_pos1 ).is_RNA() && pose.residue( jump_pos2 ).is_RNA()) ) {
			success = true;
			continue;
		}

		if ( moveable_jump( jump_pos1, jump_pos2, *(jump_mover.atom_level_domain_map()) ) ) jump_mover.add_new_RNA_jump( pose, i, success );

		if ( !success && !new_fold_tree_initializer_ ) return;

	}
	if ( new_fold_tree_initializer_ ) {
		// sometimes sets up fold trees that the jump mover doesn't really understand
		// but it's not actually an issue, so let's just continue...
		success = true;
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////
// In principle, don't need to pass atom_level_domain_map into this function (its stored in rna_jump_mover),
// but I want to make sure that user known that  atom_level_domain_map will get updated -- rhiju, 2015
void
RNA_DeNovoPoseInitializer::setup_fold_tree_and_jumps_and_variants( pose::Pose & pose,
	RNA_JumpMover const & rna_jump_mover,
	core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map,
	libraries::RNA_ChunkLibrary const & rna_chunk_library,
	bool const & enumerate /*=false*/ ) const
{
	runtime_assert( rna_jump_mover.atom_level_domain_map() == atom_level_domain_map );
	setup_jumps( pose, rna_jump_mover, rna_chunk_library, enumerate );
	setup_chainbreak_variants(  pose, atom_level_domain_map );
	setup_block_stack_variants( pose, atom_level_domain_map );
	runtime_assert( rna_jump_mover.atom_level_domain_map() == atom_level_domain_map );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Only in use by cs_rosetta_rna.
// Creates a temporary JumpLibrary and AtomLevelDomainMap to do the jump setup.
void
RNA_DeNovoPoseInitializer::setup_fold_tree_and_jumps_and_variants( pose::Pose & pose ) const
{
	using namespace core::pose::toolbox;

	AtomLevelDomainMapOP atom_level_domain_map( new AtomLevelDomainMap( pose ) );
	RNA_JumpMover const rna_jump_mover( libraries::RNA_LibraryManager::get_instance()->rna_jump_library_cop(), atom_level_domain_map );
	libraries::RNA_ChunkLibrary rna_chunk_library; // blank rna_chunk_library
	setup_fold_tree_and_jumps_and_variants( pose, rna_jump_mover, atom_level_domain_map, rna_chunk_library );
}

///////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::setup_jumps( pose::Pose & pose,
	RNA_JumpMover const & rna_jump_mover,
	libraries::RNA_ChunkLibrary const & rna_chunk_library,
	bool const enumerate /*=false*/ ) const
{

	///// We need to do something differently here if we're going to use electron density
	// If we're using density, then there should be a virtual residue at the end of the sequence
	//
	//////
	using namespace core::pose::rna;

	core::Size const nres = pose.size();
	kinematics::FoldTree f( nres );


	/// NEW fold tree setup option
	if ( rna_params_.use_fold_tree_from_silent_file() ) {
		f = rna_params_.fold_tree_from_silent_file();
	} else if ( new_fold_tree_initializer_ ) { // default false
		setup_fold_tree_through_build_full_model_info( pose, rna_chunk_library, enumerate );
		f = pose.fold_tree();
	} else {
		f = setup_fold_tree_legacy( pose, rna_jump_mover );
	}

	///////////////////////////////////////////////////////////
	// Basic setup ==> How many jumps? cuts?
	///////////////////////////////////////////////////////////
	// Hold on to torsion angles in case we need to set up chainbreak residues...
	pose::Pose pose_copy = pose;

	fill_in_default_jump_atoms( f, pose );

	if ( rna_params_.virtual_anchor_attachment_points().size() > 0 ) {

		f.reorder( pose.size() ); //reroot so that virtual residue is fixed.

		if ( root_at_first_rigid_body_ ) {
			utility::vector1< core::Size > rigid_body_jumps = get_rigid_body_jumps( pose );
			runtime_assert( rigid_body_jumps.size() > 0 );
			core::Size const anchor_rsd = pose.fold_tree().upstream_jump_residue( rigid_body_jumps[1] );
			TR << "ROOTING AT RSD" << anchor_rsd << std::endl;
			f.reorder( anchor_rsd ); //reroot so that partner of virtual residue is fixed
		}

	} else if ( model_with_density_ ) {

		// if we're using density then we need to have the root at the virtual residue
		// reorder the fold tree so that it is rooted at the virtual residue
		f.reorder( pose.size() );
	} else { // I don't think we want to do this if we have protein

		// also useful -- if user has an input pdb, put the root in there, if possible.
		//  for (core::Size n = pose.size(); n >= 1; n-- ){ // not sure why I did this backwards...
		for ( core::Size n = 1; n <= pose.size(); n++ ) {
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
	bool success( true );

	// don't do this if we're using the new fold tree initializer?
	if ( !new_fold_tree_initializer_ ) {
		bool const random_jumps( true ); // For now this is true... perhaps should also have a more deterministic procedure.
		if ( random_jumps ) insert_base_pair_jumps( pose, rna_jump_mover, success );

		// argh: some of the jump atoms can get screwy after insert_base_pair_jumps.
		fill_in_default_jump_atoms( pose );

		if ( !success ) {
			utility_exit_with_message( "Trouble inserting base pair jumps into pose -- check residues and edges." );
		}
	}


}


///////////////////////////////////////////////////////////////
kinematics::FoldTree
RNA_DeNovoPoseInitializer::setup_fold_tree_legacy( pose::Pose & pose,
	RNA_JumpMover const & rna_jump_mover ) const
{

	using namespace core::pose::rna;

	core::Size const nres = pose.size();
	kinematics::FoldTree f( nres );

	core::Size const num_cuts_closed( rna_params_.cutpoints_closed().size() );
	core::Size const num_cuts_open  ( rna_params_.cutpoints_open().size() );
	core::Size const num_cuts_total ( num_cuts_closed + num_cuts_open );

	////////////////////////////////////////////////////////////////////////
	utility::vector1< utility::vector1< core::Size > > obligate_pairing_sets,  stem_pairing_sets;
	obligate_pairing_sets = rna_params_.obligate_pairing_sets();
	if ( bps_moves_ ) { // supplement obligate_pairing_sets with stems in freely moving regions.
		for ( auto const & stem_pairing_set : rna_params_.stem_pairing_sets() ) {
			for ( core::Size const idx : stem_pairing_set ) {
				BasePair const & base_pair = rna_params_.rna_pairing_list()[ idx ] ;
				if ( rna_params_.check_in_pairing_sets( obligate_pairing_sets, base_pair ) ) continue;
				if ( !base_pair_moving( base_pair, rna_jump_mover.atom_level_domain_map(), pose ) ) continue;
				obligate_pairing_sets.push_back( utility::tools::make_vector1( idx ) );
			}
		}
	} else {
		stem_pairing_sets = rna_params_.stem_pairing_sets();
	}

	core::Size const num_pairings( rna_params_.rna_pairing_list().size() );
	core::Size const num_obligate_pairing_sets( obligate_pairing_sets.size() );
	core::Size const num_stem_pairing_sets( stem_pairing_sets.size() );
	core::Size const num_chain_connections( rna_params_.chain_connections().size() );
	runtime_assert( num_stem_pairing_sets + num_obligate_pairing_sets <= num_pairings );

	//////////////////////////////////////////////////////////////////////
	// Cuts.
	//////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > obligate_cut_points;
	for ( core::Size n = 1; n<= num_cuts_closed; n++ )   obligate_cut_points.push_back( rna_params_.cutpoints_closed()[ n ] );
	for ( core::Size n = 1; n<= num_cuts_open  ; n++ )   obligate_cut_points.push_back( rna_params_.cutpoints_open()[n] );

	//////////////////////////////////////////////////////////////////////
	// base pair steps are a special kind of chunk, created "on-the-fly"
	// from a database. These stem base pairs will be obligate pairs
	//  (see above get_pairings_from_line ), and we can define cutpoints
	//  ahead of time.
	//////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > base_pair_step_starts;
	if ( bps_moves_ ) {
		RNA_BasePairHandler const rna_base_pair_handler( rna_params_ ); // has get_base_pair_steps() function.
		utility::vector1< BasePairStep > base_pair_steps = rna_base_pair_handler.get_base_pair_steps( false /* just canonical */ );

		for ( BasePairStep const & base_pair_step : base_pair_steps ) {
			if ( base_pair_step.j_next() != (base_pair_step.j() + 1) ) continue;

			// some base pair steps are actually not moveable.
			if ( !base_pair_step_moving( base_pair_step, rna_jump_mover.atom_level_domain_map(), pose ) ) continue;

			// following assert (jumps between top & bottom base pair of base pair step) do not hold if one of the base pair is inside a fixed input PDB.
			//   runtime_assert( rna_params_.check_in_pairing_sets( obligate_pairing_sets, BasePair( base_pair_step.i(),      base_pair_step.j_next() ) ) );
			//   runtime_assert( rna_params_.check_in_pairing_sets( obligate_pairing_sets, BasePair( base_pair_step.i_next(), base_pair_step.j()      ) ) );
			runtime_assert ( base_pair_step.i_next() == base_pair_step.i()+1 );

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
	// this sort of works with density -- but should really use -new_fold_tree_initializer if working with density!!
	core::Size num_density_cuts = 0;
	if ( model_with_density_ ) num_density_cuts = 1;
	core::Size const num_pairings_to_force = std::max( num_obligate_pairing_sets + num_chain_connections + num_density_cuts,
		num_cuts_total );

	////////////////////////////////////////////////////////////////////////
	ObjexxFCL::FArray2D< core::Size > jump_points( 2, num_pairings_to_force );
	ObjexxFCL::FArray1D< core::Size > cuts( num_pairings_to_force );

	//////////////////////////////////////////////////////////////////////
	// If a cut needs to be randomly chosen, will generally try to
	// place it in a loopy region.
	FArray1D_float cut_bias( nres, 1.0 );
	std::string const & rna_secstruct_legacy( get_rna_secstruct_legacy( pose ) );
	for ( core::Size i = 1; i < nres; i++ ) {
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
		core::Size domain( rna_jump_mover.atom_level_domain_map()->get_domain( core::id::AtomID( named_atom_id_to_atom_id( core::id::NamedAtomID( " P  ", i+1 ), pose )  ) ) );
		if ( domain > 0 ) {
			if ( bps_moves_ && base_pair_step_starts.has_value( i ) ) {
				cut_bias( i ) = 0.01;
			} else {
				cut_bias( i ) = 0.0;
			}
		}
	}

	// for ( core::Size i = 1; i < nres; i++ ) {
	//   TR  << TR.Blue << "CUT_BIAS " << i << " " << cut_bias( i ) << std::endl;
	// }

	//////////////////////////////////////////////////////////////////////
	// Jump residues.
	//////////////////////////////////////////////////////////////////////
	core::Size ntries( 0 );
	core::Size const MAX_TRIES( 1000 );
	bool success( false );

	while ( !success && ++ntries < MAX_TRIES ) {

		// First, put in a jump from each obligate pairing set.
		core::Size count( 0 );

		for ( core::Size n = 1; n <= num_obligate_pairing_sets ; n++ ) {
			core::Size const pairing_index_in_stem( static_cast<core::Size>( numeric::random::rg().uniform() * obligate_pairing_sets[n].size() )  + 1 );
			core::Size const which_pairing = obligate_pairing_sets[n][pairing_index_in_stem];
			count++;
			jump_points(1, count) = rna_params_.rna_pairing_list()[which_pairing].res1();
			jump_points(2, count) = rna_params_.rna_pairing_list()[which_pairing].res2();
			//    TR << "JUMPS1 " <<  jump_points(1,count) << ' ' << jump_points(2,count ) << std::endl;
		}

		// "Chain connections" provide less information about specific residues to pair --
		//  but they're assumed to be obligatory.
		for ( core::Size n = 1; n <= num_chain_connections ; n++ ) {
			utility::vector1 < core::Size > const & res_list1( rna_params_.chain_connections()[n].first );
			utility::vector1 < core::Size > const & res_list2( rna_params_.chain_connections()[n].second);
			core::Size jump_pos1 = numeric::random::rg().random_element( res_list1 );
			core::Size jump_pos2 = numeric::random::rg().random_element( res_list2 );
			count++;
			jump_points(1, count) =  std::min( jump_pos1, jump_pos2 );
			jump_points(2, count) =  std::max( jump_pos1, jump_pos2 );
			//  TR << "JUMPS2 " <<  jump_points(1,count) << ' ' << jump_points(2,count ) << std::endl;
		}
		//  TR << std::endl;

		// If we are using density, add the density jump and increase the count
		if ( model_with_density_ ) {
			// nres should be a virtual residue that was added in initial pose setup
			// Then we need to reorder the fold tree so that it is rooted at the virtual residue
			count++;
			jump_points(1, count) = nres;
			//jump_points(2, count) = nres - 1;
			jump_points(2, count) = 1;
		}

		// Then, to fill out pairings, look at remaining possible pairing sets (these
		// should typically be Watson-Crick stems, but this setup is general )
		// Note that there might be three stems defined, but we only want two --
		//  following picks a random set of two.
		//
		// Following has been ad hoc for a long time -- can be made more systematic
		//  based on fold_tree build up that has been worked out in stepwise monte carlo.
		//
		FArray1D < bool > used_set( num_stem_pairing_sets, false );
		core::Size num_sets_left( num_stem_pairing_sets );

		while ( count < num_pairings_to_force ) {

			// Find a random pairing among what's remaining.
			core::Size const set_index( static_cast<core::Size>( numeric::random::rg().uniform() * num_sets_left )  + 1 );
			core::Size m( 0 ), set_count( 0 );
			for ( m = 1; m <= num_stem_pairing_sets; m++ ) {
				if ( !used_set( m ) ) {
					set_count++;
					if ( set_count == set_index ) break;
				}
			}

			if ( m > num_stem_pairing_sets ) {
				utility_exit_with_message( "Could not find a stem_pairing. Number of stem_pairing_sets: "+I(3,num_stem_pairing_sets)+". Did you specify a secondary structure or obligate pair?" );
			}

			core::Size const pairing_index_in_set( static_cast<core::Size>( numeric::random::rg().uniform() * stem_pairing_sets[m].size() )  + 1 );
			core::Size const which_pairing = stem_pairing_sets[m][pairing_index_in_set];

			count++;
			jump_points(1, count) = rna_params_.rna_pairing_list()[which_pairing].res1();
			jump_points(2, count) = rna_params_.rna_pairing_list()[which_pairing].res2();

			used_set( m ) = true;
			num_sets_left--;
		}

		////////////////////////////////////////////////////////////////////////////////
		// Do it, get the fold tree. and set up the pose.
		////////////////////////////////////////////////////////////////////////////////
		std::vector< core::Size > obligate_cut_points_reformat;
		for ( core::Size q = 1; q <= obligate_cut_points.size(); q++ ) obligate_cut_points_reformat.push_back( obligate_cut_points[q] );

		// TR << TR.Cyan << "Making attempt " << ntries << std::endl;
		// TR << TR.Cyan << "obligate_cutpoints " << std::endl;
		// for ( core::Size q: obligate_cut_points_reformat ) TR << TR.Cyan << " " << pose.pdb_info()->chain( q ) << ":" << pose.pdb_info()->number( q ) << std::endl;
		// for (core::Size n = 1; n <= count; n++ ){
		// //for (core::Size n = 1; n <= num_pairings_to_force; n++ ){
		//   TR << TR.Cyan << "JUMPS" <<
		//   " " << pose.pdb_info()->chain( jump_points(1, n) ) << ":" << pose.pdb_info()->number( jump_points( 1, n ) ) <<
		//   " " << pose.pdb_info()->chain( jump_points(2, n) ) << ":" << pose.pdb_info()->number( jump_points( 2, n ) ) <<
		//   std::endl;
		//}

		success = f.random_tree_from_jump_points( nres, num_pairings_to_force, jump_points, obligate_cut_points_reformat, cut_bias, 1, true /*enable 1 or NRES jumps*/ );
	}


	if ( !success )  utility_exit_with_message( "Couldn't find a freaking tree!" );

	return f;
}

/////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::setup_fold_tree_through_build_full_model_info(
	pose::Pose & pose,
	libraries::RNA_ChunkLibrary const & chunk_library,
	bool const & enumerate /*=false*/ ) const
{

	using namespace core::pose::full_model_info;

	TR << "Setting up fold tree through build full model info" << std::endl;

	bool docking = false;
	if ( pose.residue( pose.size() ).name3() == "XXX" ) {
		docking = true;
	}

	pose::Pose pose_input = pose;
	pose::Pose start_pose;
	FullModelParametersCOP full_model_params_input_pose( FullModelParametersCOP( const_full_model_info( pose_input ).full_model_parameters()->clone()  ) );

	utility::vector1< pose::PoseOP > chunk_poses;
	chemical::ResidueTypeSetCOP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

	bool added_vrt_res = false; // for density stuff

	//////////////////////////////
	// Set up pose chunks (for input structures (-in:file:s, -in:file:silent))
	//////////////////////////////
	utility::vector1< utility::vector1< core::Size > > all_res_lists;

	// Create little pose chunks that correspond to each domain 1, 2, 3, ... N in the atom_level_domain_map as separate poses with reasonable fold trees.
	for ( auto chunk_set : chunk_library.chunk_sets() ) {
		utility::vector1< core::Size > res_list;
		for ( auto const & elem : chunk_set->res_map() ) res_list.push_back( elem.first );
		all_res_lists.push_back( res_list );
		// ChunkSet stores miniposes, not poses.
		pose::MiniPoseCOP mini_pose = chunk_set->mini_pose( 1 );

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// Ideally following should be a utility function that should be held as PoseOP minipose.create_pose();
		pose::PoseOP chunk_pose( new pose::Pose );
		pose::make_pose_from_sequence( *chunk_pose, mini_pose->sequence(), rsd_set );

		// add the proper variant types -- otherwise build_full_model does not work
		for ( core::Size i=1; i<= mini_pose->variant_types_list().size(); ++i ) {
			for ( core::Size j=1; j<=mini_pose->variant_types_list()[i].size(); ++j ) {
				chemical::VariantType vt = chemical::ResidueProperties::get_variant_from_string( mini_pose->variant_types_list()[i][j]);
				pose::add_variant_type_to_pose_residue( *chunk_pose, vt, i );
				//pose::add_variant_type_to_pose_residue( pose, VIRTUAL_PHOSPHATE, n+1  );
			}
		}


		chunk_pose->fold_tree( mini_pose->fold_tree() );
		// copy over coordinate from mini-pose to pose? that would allow a sanity check
		// make a resmap
		pose::ResMap resmap;
		for ( core::Size i = 1; i<=res_list.size(); ++i ) {
			resmap[i] = i;
		}
		core::pose::copydofs::copy_dofs( *chunk_pose, *mini_pose, resmap );
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// if this is the first chunk and we're modeling with density (or... other?)
		if ( (model_with_density_ || docking) && !added_vrt_res ) {
			// pose_input.residue( pose_input.size() ) better already be a virtual residue by here...
			//if the pose doesn't already have a virtual residue (I'm not really sure how it could, but just in case)
			if ( chunk_pose->residue( chunk_pose->size() ).name3() != "XXX" ) {
				pose::addVirtualResAsRoot( *chunk_pose ); // appends by jump to center of mass, I think it's OK
				// add it to the res list - it should be the last residue in the full pose
				res_list.push_back( pose_input.size() );
			}
			added_vrt_res = true;
		}

		// Then assign to each of them the same full_model_info object as in the pose_input,
		// but with res_list corresponding to the actual residues in each pose chunk.
		// (so the full_model_parameters will be the same in each full_model_info.)
		FullModelInfoOP full_model_info( new FullModelInfo( full_model_params_input_pose ) );
		full_model_info->set_res_list( res_list );
		set_full_model_info( *chunk_pose, full_model_info );
		chunk_poses.push_back( chunk_pose );
	}

	///// if we want to dock each chunk:
	///// check that we have a "chunk" from every chain, if not, make one
	if ( dock_each_chunk_ && (pose_input.residue( pose_input.size() ).name3() == "XXX") ) {
		// get a list of all the chains
		utility::vector1< std::string > chains = full_model_params_input_pose->conventional_chains();
		utility::vector1< std::string > chains_covered;
		for ( core::Size n=1; n<=chunk_poses.size(); ++n ) {
			for ( auto const i : const_full_model_info(*chunk_poses[n]).res_list() ) {
				std::string c = chains[i];
				chains_covered.push_back( c );
			}
		}

		// loop through the full list of chains, check whether they're represented in chains_covered
		// if not, make a little one residue pose for the first residue in that chain
		if ( !center_jumps_in_single_stranded_ ) {
			for ( core::Size i=1; i<=chains.size(); ++i ) {
				// If chain has already been covered, move on.
				if ( std::find( chains_covered.begin(), chains_covered.end(), chains[i] ) != chains_covered.end() ) continue;

				// it's not covered, make a 1 residue pose for the chain
				pose::PoseOP little_pose( new pose::Pose );
				std::string seq = full_model_params_input_pose->full_sequence().substr( i-1, 1 );

				pose::make_pose_from_sequence( *little_pose, seq, rsd_set );
				// add lower terminus variant type
				pose::add_variant_type_to_pose_residue( *little_pose, chemical::LOWER_TERMINUS_VARIANT, 1 );

				// add full model info
				utility::vector1< core::Size > res_list;
				res_list.push_back( i );
				FullModelInfoOP full_model_info( new FullModelInfo( full_model_params_input_pose ) );
				full_model_info->set_res_list( res_list );
				set_full_model_info( *little_pose, full_model_info );
				chunk_poses.push_back( little_pose );

				chains_covered.push_back( chains[i] ); // add the chain to the list of chains covered
			}
		} else { // center_jumps_in_single_stranded_

			// if the chain is all single stranded:
			// add the middle residue in the chain, rather than the first residue in the chain
			std::string const & rna_secstruct_legacy( rna_params_.rna_secstruct_legacy() );

			for ( core::Size i=1; i<=chains.size(); ++i ) {
				// If chain has already been covered, move on.
				if ( std::find( chains_covered.begin(), chains_covered.end(), chains[i] ) != chains_covered.end() ) continue;

				// it's not covered, make a 1 residue pose for the chain
				// this is the start of the chain then
				core::Size chain_start = i;
				core::Size chain_stop = i;
				// find the end of the chain
				for ( core::Size j=i; j<=chains.size(); ++j ) {
					if ( chains[j] != chains[i] ) {
						chain_stop = j-1;
						break;
					}
					if ( j == chains.size() ) chain_stop = j;
				}

				// figure out whether there's any secondary structure in this region
				// if not, add the middle residue
				std::string secstruct_chain = rna_secstruct_legacy.substr( chain_start -1, chain_stop-chain_start +1 );
				core::Size resid_to_add = i;
				if ( std::find( secstruct_chain.begin(), secstruct_chain.end(), 'H' ) == secstruct_chain.end() ) {
					resid_to_add = i + (( chain_stop - chain_start ) / 2);
				}

				pose::PoseOP little_pose( new pose::Pose );
				std::string seq = full_model_params_input_pose->full_sequence().substr( resid_to_add-1, 1 );

				pose::make_pose_from_sequence( *little_pose, seq, rsd_set );
				// add lower terminus variant type if we added the first residue in the chain
				if ( resid_to_add == chain_start ) {
					pose::add_variant_type_to_pose_residue( *little_pose, chemical::LOWER_TERMINUS_VARIANT, 1 );
				}
				if ( resid_to_add == chain_stop ) {
					pose::add_variant_type_to_pose_residue( *little_pose, chemical::UPPER_TERMINUS_VARIANT, 1 );
				}

				// add full model info
				utility::vector1< core::Size > res_list;
				res_list.push_back( resid_to_add );
				FullModelInfoOP full_model_info( new FullModelInfo( full_model_params_input_pose ) );
				full_model_info->set_res_list( res_list );
				set_full_model_info( *little_pose, full_model_info );
				chunk_poses.push_back( little_pose );

				chains_covered.push_back( chains[i] ); // add the chain to the list of chains covered
			}
		}

	}



	// /// dump out the chunk poses to check
	// for ( core::Size i=1; i<=chunk_poses.size(); ++i ) {
	//  std::ostringstream name;
	//  name << "chunk_pose_" << i << ".pdb";
	//  chunk_poses[i]->dump_pdb( name.str() );
	//  std::cout << "chunk_pose " << i << " " << chunk_poses[i]->annotated_sequence() << std::endl;
	//  chunk_poses[i]->fold_tree().show( std::cout );
	// }


	// Also create A-form helices for 'stems'. Use RNA_HelixAssembler.
	// Actually, punt on this -- we'll fill this in later.
	// yeah that's definitely going to be necessary though...
	// pose setup will fail if you have stems but don't input structures for them

	//////////////////////////////
	// Set up the starting pose
	//////////////////////////////

	if ( dock_each_chunk_ && (chunk_poses.size() > 0) && (pose_input.residue( pose_input.size() ).name3() == "XXX") ) {
		// merge all the chunks into 1 pose
		start_pose = *chunk_poses[ 1 ];

		// figure out which chunks we're going to merge (don't merge in multiple chunks from the same chain)
		utility::vector1< std::string > all_chains = full_model_params_input_pose->conventional_chains();
		utility::vector1< std::string > chains_in_start_pose;
		// a vector of booleans telling us whether to merge chunk poses or not, same length as chunk_poses
		utility::vector1< bool > chunks_to_merge;
		utility::vector1< pose::PoseOP > remaining_chunk_poses; // the chunk poses that aren't merged - will be included as "other poses"
		chunks_to_merge.push_back( true ); // include the first chunk

		for ( auto const i : const_full_model_info(*chunk_poses[1]).res_list() ) {
			std::string c = all_chains[i];
			chains_in_start_pose.push_back( c );
		}

		// loop through the remaining chunks
		for ( core::Size n = 2; n<= chunk_poses.size(); ++n ) {
			// if we already have a chunk from this chain, then we don't want to add this chunk
			// unless we specified dock_each_chunk_per_chain
			// add it to a new list of chunk poses
			utility::vector1< std::string > chains_in_chunk;
			bool merge_chunk = true;
			for ( auto const i : const_full_model_info( *chunk_poses[n] ).res_list() ) {
				std::string c = all_chains[ i ];
				chains_in_chunk.push_back( c );
				// check if this chain is already represented in the start pose, if so, we don't want to add this chunk as well
				if ( std::find( chains_in_start_pose.begin(), chains_in_start_pose.end(), c ) != chains_in_start_pose.end() ) { //HERE
					merge_chunk = false;
					if ( !dock_each_chunk_per_chain_ ) {
						remaining_chunk_poses.push_back( chunk_poses[n] );
					}
					break;
				}
			}
			// check whether this is a dock_chunks_res --> if yes, merge_chunk = true (but if not, merge_chunk could still be true)
			for ( auto const i : const_full_model_info( *chunk_poses[n] ).res_list() ) {
				if ( std::find( dock_chunks_res_.begin(), dock_chunks_res_.end(), i ) != dock_chunks_res_.end() ) {
					merge_chunk = true;
				}
			}
			// we're including the chunk, add the chains to the chains_in_start_pose
			if ( merge_chunk ) {
				for ( auto const & i : chains_in_chunk ) {
					chains_in_start_pose.push_back( i );
				}
			}
			chunks_to_merge.push_back( merge_chunk );
		}

		// for testing, go ack to merging all chunks
		if ( dock_each_chunk_per_chain_ ) {
			for ( core::Size n =1; n<=chunks_to_merge.size(); ++n ) {
				chunks_to_merge[n] = true;
			}
		}



		// chunk_poses[1] better have a virtual residue at the end
		core::Size vrt_resnum = start_pose.size();
		core::Size final_vrt_resnum = pose.size(); // full pose numbering!

		utility::vector1< core::Size > res_list_all_chunks_for_mapping;
		std::map< core::Size, std::pair< core::Size, core::Size> > full_number_to_chunk_num_and_resid; // full pose number: (chunk pose number, resid)
		for ( core::Size n = 1; n<= chunk_poses.size(); ++n ) {
			if ( !chunks_to_merge[n] ) continue;
			for ( core::Size j =1; j<=const_full_model_info( *chunk_poses[n] ).res_list().size(); ++j ) {
				res_list_all_chunks_for_mapping.push_back( const_full_model_info( *chunk_poses[n] ).res_list()[j] );
				full_number_to_chunk_num_and_resid[ const_full_model_info( *chunk_poses[n] ).res_list()[j] ] = std::make_pair( n, j );
			}
		}
		// reorder it
		std::sort( res_list_all_chunks_for_mapping.begin(), res_list_all_chunks_for_mapping.end() );
		pose::ResMap resmap_from_full_to_merged_chunks; /// full merged chunk number: (chunk_pose number, resid)
		std::map< core::Size, std::pair< core::Size, core::Size > > map_merged_chunk_res_to_chunk_poses;
		for ( core::Size i = 1; i<=res_list_all_chunks_for_mapping.size(); ++i ) {
			resmap_from_full_to_merged_chunks[ res_list_all_chunks_for_mapping[i] ] = i;
			map_merged_chunk_res_to_chunk_poses[ i ] = full_number_to_chunk_num_and_resid[ res_list_all_chunks_for_mapping[i] ];
		}

		// Set up the initial fold tree
		kinematics::FoldTree new_fold_tree;
		utility::vector1< core::Size > const res_list_start = const_full_model_info( start_pose ).res_list();
		// add the edges from the start_pose
		auto it( start_pose.fold_tree().begin() );
		for ( auto it_end = start_pose.fold_tree().end(); it != it_end; ++it ) {
			// just check that it doesn't include the virtual residue
			// (which should only be attached by a jump!!!)
			if ( (it->start() == vrt_resnum) || (it->stop() == vrt_resnum ) ) {
				// add the edge, but replace the vrt_resnum with the final vrt_resnum
				if ( it->start() == vrt_resnum ) {
					new_fold_tree.add_edge( resmap_from_full_to_merged_chunks[ final_vrt_resnum],
						resmap_from_full_to_merged_chunks[ res_list_start[it->stop()]], it->label() );
				} else {
					new_fold_tree.add_edge( resmap_from_full_to_merged_chunks[ res_list_start[it->start()]],
						resmap_from_full_to_merged_chunks[ final_vrt_resnum ], it->label() );
				}
			} else {
				// just add the edge as is
				new_fold_tree.add_edge( resmap_from_full_to_merged_chunks[ res_list_start[it->start()]],
					resmap_from_full_to_merged_chunks[ res_list_start[it->stop()]], it->label() );
				//new_fold_tree.add_edge( res_list_start[it->start()], res_list_start[it->stop()], it->label() );
			}
		}
		//// set the jump atoms
		//for ( core::Size n = 1; n<=start_pose.fold_tree().num_jump(); ++n ) {
		// new_fold_tree.set_jump_atoms( n, start_pose.fold_tree().upstream_atom( n ), start_pose.fold_tree().downstream_atom( n ), false );
		//}
		//fill_in_default_jump_atoms( new_fold_tree, start_pose );


		///// CHECK jump atoms
		//for ( core::Size n=1; n<=start_pose.fold_tree().num_jump(); ++n ) {
		// std::cout << "f.upstream_atom(" << n << ") " <<  start_pose.fold_tree().upstream_atom(n) << std::endl;
		// std::cout << "f.downstream_atom(" << n << ") "<< start_pose.fold_tree().downstream_atom(n) << std::endl;
		//}
		//std::cout << "And check the new fold tree as well" << std::endl;
		//new_fold_tree.show( std::cout );
		//for ( core::Size n=1; n<=new_fold_tree.num_jump(); ++n ) {
		// std::cout << "new_fold_tree.upstream_atom(" << n << ") " <<  new_fold_tree.upstream_atom(n) << std::endl;
		// std::cout << "new_fold_tree.downstream_atom(" << n << ") "<< new_fold_tree.downstream_atom(n) << std::endl;
		//}



		// Merge in the rest of the chunk poses
		for ( core::Size n = 2; n<= chunk_poses.size(); ++n ) {

			if ( !chunks_to_merge[n] ) {
				continue;
			}

			utility::vector1< core::Size > const & res_list_chunk = const_full_model_info( *chunk_poses[n] ).res_list();

			core::Size root_resnum = resmap_from_full_to_merged_chunks[ res_list_chunk[ chunk_poses[n]->fold_tree().root() ] ];
			//core::Size init_size = start_pose.size() -1; // don't include virtual residue

			// add everything to the new fold tree
			// figure out where in the pose this chunk actually belongs

			// add all the edges
			auto it ( chunk_poses[n]->fold_tree().begin() );
			for ( auto it_end = chunk_poses[n]->fold_tree().end(); it != it_end; ++it ) {
				// figure out where this edge belongs (i.e. get correct start and stop residues)
				core::Size edge_label = -1;
				if ( it->label() != -1 ) {
					edge_label = new_fold_tree.num_jump()+1;
				}
				new_fold_tree.add_edge( resmap_from_full_to_merged_chunks[ res_list_chunk[it->start()]],
					resmap_from_full_to_merged_chunks[ res_list_chunk[it->stop()] ], edge_label );
				//new_fold_tree.add_edge( it->start() + init_size, it->stop() + init_size, edge_label );
			}

			// add jump from the vrt res to the old ft root
			new_fold_tree.add_edge( resmap_from_full_to_merged_chunks[ final_vrt_resnum ], root_resnum, new_fold_tree.num_jump()+1 );

			// Done: new fold tree should be updated now with this additional chunk

		}
		// root the fold tree on the virtual residue

		new_fold_tree.reorder( res_list_all_chunks_for_mapping.size() );
		//new_fold_tree.reorder( start_pose.size() );

		// we need to preserve the jump atoms
		// f.set_jump_atoms( jump_num, atom_name, f.downstream_atom( jump_num ), false )
		// f.set_jump_atoms( jump_num, f.upstream_atom( jump_num ), atom_name, false )

		// let's actually make the pose
		// map_merged_chunk_res_to_chunk_poses
		pose::Pose new_start_pose;
		// let's add the first residue
		core::Size chunk_num_first = map_merged_chunk_res_to_chunk_poses[ 1 ].first;
		core::Size resid_first = map_merged_chunk_res_to_chunk_poses[ 1 ].second;
		new_start_pose.append_residue_by_jump( chunk_poses[ chunk_num_first ]->residue( resid_first ), 1 );
		// then add the virtual residue -- always on chunk # 1
		new_start_pose.append_residue_by_jump( chunk_poses[ 1 ]->residue( chunk_poses[1]->size() ), 1 );
		//new_start_pose.append_residue_by_jump( chunk_poses[ chunk_num_first ]->residue( chunk_poses[chunk_num_first]->size() ), 1 );
		//new_start_pose.conformation().append_residue_by_jump( pose.residue( pose.size() ), 1 );


		// now actually merge the poses
		// we really need the root_resnum to be 1
		utility::vector1< core::Size > residues_added; // in merged chunk numbering
		//core::Size merged_chunks_final_vrt_resnum = res_list_all_chunks_for_mapping.size();

		// insert_residue_by_jump( residue, desired_seq_pos, anchor in current numbering)
		utility::vector1< core::Size > jump_residues;
		kinematics::FoldTree const & nft = new_fold_tree;
		auto it_ft ( nft.begin() );
		for ( auto it_end_ft = nft.end(); it_ft != it_end_ft; ++it_ft ) {
			if ( it_ft->label() == -1 ) continue; // not a jump
			if ( !(std::find( jump_residues.begin(), jump_residues.end(), it_ft->start())!=jump_residues.end()) ) {
				jump_residues.push_back( it_ft->start() );
			}
			if ( !(std::find( jump_residues.begin(), jump_residues.end(), it_ft->stop())!=jump_residues.end()) ) {
				jump_residues.push_back( it_ft->stop() );
			}

		}

		// let's go through the residues in the full merged chunk numbering ( i.e. 1 to size of full merged chunk )
		for ( core::Size i=2; i<=map_merged_chunk_res_to_chunk_poses.size()-1; ++i ) { // already added virtual
			core::Size chunk_num = map_merged_chunk_res_to_chunk_poses[ i ].first;
			core::Size resid = map_merged_chunk_res_to_chunk_poses[ i ].second;

			// if it's a jump residue, insert it by a jump
			if ( std::find( jump_residues.begin(), jump_residues.end(), i)!=jump_residues.end() ) {
				// it's not really the correct jump, but should be OK b/c we reset with fold tree after
				new_start_pose.insert_residue_by_jump( chunk_poses[chunk_num]->residue( resid ), new_start_pose.size(), new_start_pose.size() );
			} else { // otherwise, insert it by a bond, just to previous residue
				new_start_pose.insert_residue_by_bond( chunk_poses[chunk_num]->residue( resid ), new_start_pose.size(), new_start_pose.size()-1);
			}
			// the fold tree is not going to be correct, but OK because we're resetting it right after this

		}

		start_pose = new_start_pose;
		// res_list_all_chunks_for_mapping is sorted

		// apply the fold tree to the pose
		fill_in_default_jump_atoms( new_fold_tree, start_pose );
		start_pose.fold_tree( new_fold_tree );

		// now update the full model info (the res list is messed up, just had info from first chunk)
		FullModelInfoOP full_model_info( new FullModelInfo( full_model_params_input_pose ) ); // keep the params from the input pose
		full_model_info->set_res_list( res_list_all_chunks_for_mapping ); // new res list with all chunk residues
		set_full_model_info( start_pose, full_model_info );

		// and now we need to add the remaining chunk poses as "other poses" (chunks that don't need to be docked individually into density)
		utility::vector1< pose::PoseOP > other_poses;
		for ( core::Size n = 1; n<=remaining_chunk_poses.size(); ++n ) {
			other_poses.push_back( remaining_chunk_poses[n] );
		}
		nonconst_full_model_info( start_pose ).set_other_pose_list( other_poses );

		// We need to set up the dock domain map
		//
		std::string const clean_desired_seq = core::pose::rna::remove_bracketed( const_full_model_info( start_pose ).full_sequence() );
		core::Size const desired_nres = clean_desired_seq.size();

		utility::vector1< core::Size > nonconst_cutpoints_who_fucking_cares = nonconst_full_model_info( start_pose ).cutpoint_open_in_full_model();
		utility::vector1< core::Size > const dock_domain_map = core::import_pose::figure_out_dock_domain_map( nonconst_cutpoints_who_fucking_cares,
			all_res_lists,
			const_full_model_info( start_pose ).working_res(),
			const_full_model_info( start_pose ).sample_res(),
			desired_nres );
		FullModelParametersOP fmp = nonconst_full_model_info( start_pose ).full_model_parameters()->clone();
		fmp->set_parameter( DOCK_DOMAIN,   dock_domain_map  );
		nonconst_full_model_info( start_pose ).set_full_model_parameters( fmp );

	} else if ( chunk_poses.size() > 0 ) {
		// make poses 2, 3, ... "other poses" in pose 1.
		start_pose = *chunk_poses[ 1 ];
		utility::vector1< pose::PoseOP > other_poses;
		for ( core::Size n = 2; n <= chunk_poses.size(); ++n ) {
			other_poses.push_back( chunk_poses[n] );
		}
		// and now make these the other poses
		nonconst_full_model_info( start_pose ).set_other_pose_list( other_poses );
		fill_in_default_jump_atoms( start_pose );


		std::string const clean_desired_seq = core::pose::rna::remove_bracketed( const_full_model_info( start_pose ).full_sequence() );
		core::Size const desired_nres = clean_desired_seq.size();

		utility::vector1< utility::vector1< core::Size > > all_res_lists;
		all_res_lists.push_back( const_full_model_info( start_pose ).res_list() );
		/*for ( core::Size n = 1; n<=remaining_chunk_poses.size(); ++n ) {
		all_res_lists.push_back( remaining_chunk_poses[n] );
		}*/
		utility::vector1< core::Size > nonconst_cutpoints_who_fucking_cares = nonconst_full_model_info( start_pose ).cutpoint_open_in_full_model();
		utility::vector1< core::Size > const dock_domain_map = core::import_pose::figure_out_dock_domain_map( nonconst_cutpoints_who_fucking_cares,
			all_res_lists,
			const_full_model_info( start_pose ).working_res(),
			const_full_model_info( start_pose ).sample_res(),
			desired_nres );
		FullModelParametersOP fmp = nonconst_full_model_info( start_pose ).full_model_parameters()->clone();
		fmp->set_parameter( DOCK_DOMAIN,   dock_domain_map  );
		nonconst_full_model_info( start_pose ).set_full_model_parameters( fmp );



	} else { // THIS DOESN'T REALLY WORK YET!
		// no chunk_poses -- so just create pose denovo.
		TR.Warning << "-new_fold_tree_initializer is not yet set up to work without any -s provided (please provide at least one -s structure or use legacy fold tree set up" << std::endl;
		pose::make_pose_from_sequence( start_pose, pose.sequence(), rsd_set );
		// and set the full_model_info
		// can't we just use the input_pose though?
		// start_pose = pose_input;
		//pose = pose_input;
	}

	// Do it: build the rest of the model, which will set up the full fold tree
	pose::PoseOP full_model_pose( new pose::Pose );
	stepwise::monte_carlo::build_full_model( start_pose, *full_model_pose, enumerate /*enumerate -- this will force stochastic choice*/, 0.0 /*bulge skip freq, don't skip bulges!*/ );
	pose = *full_model_pose;

}

///////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::setup_chainbreak_variants( pose::Pose & pose,
	core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map ) const
{
	pose::Pose pose_copy = pose;
	utility::vector1< core::Size > const & cutpoints_open( rna_params_.cutpoints_open() );

	// Create cutpoint variants to force chainbreak score computation.
	for ( core::Size cutpos = 1; cutpos < pose.size(); cutpos++ ) {

		if ( ! pose.fold_tree().is_cutpoint( cutpos ) ) continue;

		// Don't assign a chainbreak penalty if user said this was an "open" cutpoint.
		if ( cutpoints_open.has_value( cutpos ) ) continue;

		// Something happened and now this is sometimes called on a jump to VRT.
		// This shouldn't be necessary: maybe the FT setup is beign screwed up.
		if ( !pose.residue_type( cutpos ).is_polymer() || ( cutpos < pose.size() && !pose.residue_type( cutpos + 1 ).is_polymer() ) ) continue;

		core::pose::correctly_add_cutpoint_variants( pose, cutpos );

		for ( core::Size i = cutpos; i <= cutpos + 1; i++ ) {
			for ( core::Size j = 1; j <= pose.residue( i ).mainchain_torsions().size(); j++ ) {
				id::TorsionID torsion_id( i, id::BB, j );
				pose.set_torsion( torsion_id, pose_copy.torsion( torsion_id ) ) ;
			} // j
		} // i
	}

	for ( core::Size const n : rna_params_.cutpoints_cyclize() ) {
		// Need to track partner -- move back along chain until we see open cutpoint.
		core::Size m;
		for ( m = n; m > 1; m-- ) {
			if ( cutpoints_open.has_value( m - 1 ) ) break;
		}
		TR << TR.Green << "Cyclizing: " << n << " to " << m << std::endl;
		core::pose::correctly_add_cutpoint_variants( pose, n, true, m );
	}

	for ( core::Size const n : rna_params_.twoprime() ) {
		// Need to track partner -- move back along chain until we see open cutpoint.
		core::Size m = 1;
		for ( m = n; m > 1; m-- ) {
			if ( cutpoints_open.has_value( m - 1 ) ) break;
		}
		TR << TR.Green << "Twoprime-cyclizing: " << n << " to " << m << std::endl;
		core::pose::correctly_add_2prime_connection_variants( pose, n, m );
	}

	// Oh, this is a bummer. If this is empty, we have to separately check --
	// since uint(0-1) > 1
	for ( core::Size i = 1; !rna_params_.fiveprime_cap().empty() && i <= rna_params_.fiveprime_cap().size() - 1; i += 2 ) {

		// Add 5PrimeCap to first, assert second is 7MG, add cutpoint lower to it.
		core::Size const capped = rna_params_.fiveprime_cap()[i];
		core::Size const methylg = rna_params_.fiveprime_cap()[i+1];
		runtime_assert( pose.residue_type( methylg ).base_name() == "7MG" );

		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::THREE_PRIME_PHOSPHATE, capped );
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::UPPER_TERMINUS_VARIANT, capped );
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::LOWER_TERMINUS_VARIANT, capped );
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, capped );
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::FIVE_PRIME_PHOSPHATE, capped );
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::FIVEPRIME_CAP, capped );


		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::LOWER_TERMINUS_VARIANT, methylg );
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, methylg );
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::FIVE_PRIME_PHOSPHATE, methylg );
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, methylg );
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::UPPER_TERMINUS_VARIANT, methylg );

		pose.conformation().declare_chemical_bond( capped, "ZO3'", methylg, "P" );
		// Not obvious why I have to re-do this, but maybe it's because I nuke some polymeric variants?
		pose.conformation().declare_chemical_bond( capped, "O3'", capped+1, "P" );
	}

	atom_level_domain_map->renumber_after_variant_changes( pose );

}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::setup_block_stack_variants(
	pose::Pose & pose,
	core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map ) const
{
	using namespace chemical;
	for ( auto m : rna_params_.block_stack_above_res() ) add_variant_type_to_pose_residue( pose, BLOCK_STACK_ABOVE, m );
	for ( auto m : rna_params_.block_stack_below_res() ) add_variant_type_to_pose_residue( pose, BLOCK_STACK_BELOW, m );
	atom_level_domain_map->renumber_after_variant_changes( pose );
}


////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoPoseInitializer::setup_virtual_phosphate_variants( pose::Pose & pose ) const {

	using namespace id;
	using namespace chemical;

	// Special exception for FIVEPRIME_CAP use in rna_denovo
	if ( pose.residue( 1 ).is_RNA() && !pose.residue_type( 1 ).has_variant_type( FIVEPRIME_CAP ) ) pose::add_variant_type_to_pose_residue( pose, VIRTUAL_PHOSPHATE, 1  );

	utility::vector1< core::Size > const & cutpoints_open( rna_params_.cutpoints_open() );
	for ( core::Size i = 1; i <= cutpoints_open.size(); i++ ) {

		core::Size n = cutpoints_open[ i ];

		if ( n == pose.size() ) {
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

} //denovo
} //rna
} //protocols


